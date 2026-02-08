"""
Mortality model for FVS-Python.

Supports multiple FVS variants:
- SN (Southern): SDI-based mortality from Section 5.0 / EFVS 7.3.2
- LS (Lake States): 4-group background mortality from morts.f with VARADJ shade tolerance
- NE (Northeast): 4-group background mortality (same model as LS, different species mappings)

Implements FVS mortality equations:
- 5.0.1: Background mortality rate
- 5.0.2: Cycle adjustment
- 5.0.3: Mortality distribution
- 5.0.4: Final mortality calculation
"""
import math
import random
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple, TYPE_CHECKING

from .tree_utils import calculate_tree_basal_area, calculate_stand_basal_area

if TYPE_CHECKING:
    from .tree import Tree


@dataclass
class MortalityResult:
    """Result of mortality application.

    Attributes:
        survivors: List of trees that survived
        mortality_count: Number of trees that died
        trees_died: List of trees that died (for volume accounting)
    """
    survivors: List['Tree']
    mortality_count: int
    trees_died: List['Tree'] = field(default_factory=list)


class MortalityModel:
    """FVS SDI-based mortality model.

    Implements the full FVS mortality model from Section 5.0 and EFVS 7.3.2:
    1. Below 55% of SDImax: Use individual tree background mortality rates
    2. Above 55% of SDImax: Use stand-level density-related mortality

    Attributes:
        default_species: Default species code for coefficient lookup
        _mortality_coefficients: Cached mortality coefficients
    """

    # Class-level cache for mortality coefficients
    _coefficients_cache: Optional[Dict[str, Any]] = None
    _coefficients_loaded: bool = False

    # Default SDI threshold constants per EFVS 7.3.2
    LOWER_THRESHOLD = 0.55  # 55% of SDImax - density-related mortality begins
    UPPER_THRESHOLD = 0.85  # 85% of SDImax - asymptotic maximum density

    def __init__(self, default_species: str = 'LP', max_sdi: Optional[float] = None):
        """Initialize the mortality model.

        Args:
            default_species: Default species code for coefficient lookups
            max_sdi: Maximum SDI for the stand (if None, uses species default)
        """
        self.default_species = default_species
        self.max_sdi = max_sdi

        # Load coefficients if not already loaded
        if not MortalityModel._coefficients_loaded:
            self._load_mortality_coefficients()

    @classmethod
    def _load_mortality_coefficients(cls) -> None:
        """Load mortality coefficients from configuration."""
        try:
            from .config_loader import get_config_loader

            loader = get_config_loader()
            mortality_file = loader.cfg_dir / "sn_mortality_model.json"
            mortality_data = loader._load_config_file(mortality_file)

            background = {}
            mwt = {}

            # Extract background mortality coefficients (Table 5.0.1)
            if 'tables' in mortality_data and 'table_5_0_1' in mortality_data['tables']:
                coeffs = mortality_data['tables']['table_5_0_1'].get('coefficients', {})
                for species, values in coeffs.items():
                    background[species] = (values['p0'], values['p1'])

            # Extract MWT values (Table 5.0.2)
            if 'tables' in mortality_data and 'table_5_0_2' in mortality_data['tables']:
                mwt = mortality_data['tables']['table_5_0_2'].get('mwt_values', {})

            cls._coefficients_cache = {'background': background, 'mwt': mwt}
            cls._coefficients_loaded = True

        except Exception:
            # Return default values if config not available
            cls._coefficients_cache = cls._get_default_coefficients()
            cls._coefficients_loaded = True

    @staticmethod
    def _get_default_coefficients() -> Dict[str, Any]:
        """Get default mortality coefficients.

        Returns:
            Dictionary with default background and mwt values
        """
        return {
            'background': {
                'LP': (5.5876999, -0.0053480),
                'SP': (5.5876999, -0.0053480),
                'SA': (5.5876999, -0.0053480),
                'LL': (5.5876999, -0.0053480),
            },
            'mwt': {
                'LP': 0.7, 'SP': 0.7, 'SA': 0.7, 'LL': 0.7,
            }
        }

    def get_coefficients(self) -> Dict[str, Any]:
        """Get the mortality coefficients.

        Returns:
            Dictionary with 'background' and 'mwt' coefficients
        """
        if self._coefficients_cache is None:
            self._load_mortality_coefficients()
        return self._coefficients_cache

    def apply_mortality(
        self,
        trees: List['Tree'],
        cycle_length: int = 5,
        max_sdi: Optional[float] = None,
        random_seed: Optional[int] = None
    ) -> MortalityResult:
        """Apply FVS SDI-based mortality model.

        Implements the full FVS mortality model from Section 5.0 and EFVS 7.3.2.

        Args:
            trees: List of trees in the stand
            cycle_length: Length of projection cycle in years
            max_sdi: Maximum SDI (uses species default if None)
            random_seed: Optional seed for reproducibility

        Returns:
            MortalityResult with survivors and mortality count
        """
        if len(trees) <= 1:
            return MortalityResult(survivors=list(trees), mortality_count=0, trees_died=[])

        if random_seed is not None:
            random.seed(random_seed)

        # Use provided max_sdi or instance default
        max_sdi = max_sdi or self.max_sdi or 450

        # Calculate stand SDI using Reineke's equation
        current_sdi = self._calculate_stand_sdi(trees)
        relative_sdi = current_sdi / max_sdi

        # Calculate basal area percentile ranking for each tree
        tree_data = self._calculate_tree_percentiles(trees)

        # Get mortality coefficients
        coefficients = self.get_coefficients()

        # Apply appropriate mortality model based on density
        if relative_sdi <= self.LOWER_THRESHOLD:
            result = self._apply_background_mortality(
                tree_data, coefficients, cycle_length
            )
        else:
            result = self._apply_density_mortality(
                tree_data, coefficients, cycle_length,
                current_sdi, max_sdi, relative_sdi
            )

        return result

    def _calculate_stand_sdi(self, trees: List['Tree']) -> float:
        """Calculate stand SDI using Reineke's equation.

        Args:
            trees: List of trees

        Returns:
            Stand Density Index
        """
        if not trees:
            return 0.0

        tpa = len(trees)
        qmd_squared = sum(tree.dbh ** 2 for tree in trees) / tpa
        qmd = math.sqrt(qmd_squared)
        return tpa * (qmd / 10.0) ** 1.605

    def _calculate_tree_percentiles(
        self,
        trees: List['Tree']
    ) -> List[Tuple['Tree', float]]:
        """Calculate basal area percentile for each tree.

        Args:
            trees: List of trees

        Returns:
            List of (tree, percentile) tuples sorted by DBH
        """
        total_ba = calculate_stand_basal_area(trees)
        tree_data = []
        cumulative_ba = 0.0
        sorted_trees = sorted(trees, key=lambda t: t.dbh)

        for tree in sorted_trees:
            tree_ba = calculate_tree_basal_area(tree.dbh)
            cumulative_ba += tree_ba
            pct = (cumulative_ba / total_ba) * 100.0 if total_ba > 0 else 50.0
            tree_data.append((tree, pct))

        return tree_data

    def _apply_background_mortality(
        self,
        tree_data: List[Tuple['Tree', float]],
        coefficients: Dict[str, Any],
        cycle_length: int
    ) -> MortalityResult:
        """Apply background mortality only (below SDI threshold).

        Uses equations 5.0.1 and 5.0.2.

        Args:
            tree_data: List of (tree, percentile) tuples
            coefficients: Mortality coefficients
            cycle_length: Cycle length in years

        Returns:
            MortalityResult
        """
        survivors = []
        trees_died = []

        for tree, pct in tree_data:
            # Get species-specific coefficients (Equation 5.0.1)
            p0, p1 = coefficients['background'].get(
                tree.species, (5.5876999, -0.0053480)
            )

            # Individual tree background mortality rate (Equation 5.0.1)
            ri = 1.0 / (1.0 + math.exp(p0 + p1 * tree.dbh))

            # Adjust for cycle length (Equation 5.0.2)
            rip = 1.0 - ((1.0 - ri) ** cycle_length)

            # Apply mortality stochastically
            if random.random() > rip:
                survivors.append(tree)
            else:
                trees_died.append(tree)

        return MortalityResult(
            survivors=survivors,
            mortality_count=len(trees_died),
            trees_died=trees_died
        )

    def _apply_density_mortality(
        self,
        tree_data: List[Tuple['Tree', float]],
        coefficients: Dict[str, Any],
        cycle_length: int,
        current_sdi: float,
        max_sdi: float,
        relative_sdi: float
    ) -> MortalityResult:
        """Apply density-related mortality (above SDI threshold).

        Uses equations 5.0.1-5.0.4.

        Args:
            tree_data: List of (tree, percentile) tuples
            coefficients: Mortality coefficients
            cycle_length: Cycle length in years
            current_sdi: Current stand SDI
            max_sdi: Maximum SDI
            relative_sdi: Current SDI / max SDI

        Returns:
            MortalityResult
        """
        survivors = []
        trees_died = []

        # Calculate target SDI reduction
        if relative_sdi > self.UPPER_THRESHOLD:
            # Need to reduce density to asymptotic level
            target_sdi = self.UPPER_THRESHOLD * max_sdi
            excess_sdi = current_sdi - target_sdi
            sdi_to_remove = excess_sdi * 0.5  # Remove up to half per cycle
        else:
            # Between lower and upper thresholds - gradual transition
            target_sdi = current_sdi * (
                1.0 - 0.05 * (relative_sdi - self.LOWER_THRESHOLD) /
                (self.UPPER_THRESHOLD - self.LOWER_THRESHOLD)
            )
            sdi_to_remove = current_sdi - target_sdi

        # Calculate target removal fraction
        target_removal_fraction = sdi_to_remove / current_sdi if current_sdi > 0 else 0

        for tree, pct in tree_data:
            # Mortality distribution by relative height/size (Equation 5.0.3)
            mr = 0.84525 - (0.01074 * pct) + (0.0000002 * (pct ** 3))
            mr = max(0.01, min(1.0, mr))

            # Get species-specific mortality weight (MWT from Table 5.0.2)
            mwt = coefficients['mwt'].get(tree.species, 0.7)

            # Final mortality rate (Equation 5.0.4)
            mort = mr * mwt * 0.1 * target_removal_fraction * (cycle_length / 5.0)

            # Also add background mortality component
            p0, p1 = coefficients['background'].get(
                tree.species, (5.5876999, -0.0053480)
            )
            ri = 1.0 / (1.0 + math.exp(p0 + p1 * tree.dbh))
            rip = 1.0 - ((1.0 - ri) ** cycle_length)

            # Combined mortality probability (independent events)
            total_mort_prob = 1.0 - (1.0 - mort) * (1.0 - rip)

            # Apply mortality stochastically
            if random.random() > total_mort_prob:
                survivors.append(tree)
            else:
                trees_died.append(tree)

        return MortalityResult(
            survivors=survivors,
            mortality_count=len(trees_died),
            trees_died=trees_died
        )

    def calculate_background_mortality_rate(
        self,
        tree: 'Tree',
        cycle_length: int = 5
    ) -> float:
        """Calculate background mortality rate for a single tree.

        Implements equations 5.0.1 and 5.0.2.

        Args:
            tree: Tree to calculate mortality for
            cycle_length: Cycle length in years

        Returns:
            Mortality probability (0-1)
        """
        coefficients = self.get_coefficients()
        p0, p1 = coefficients['background'].get(
            tree.species, (5.5876999, -0.0053480)
        )

        # Individual tree background mortality rate (Equation 5.0.1)
        ri = 1.0 / (1.0 + math.exp(p0 + p1 * tree.dbh))

        # Adjust for cycle length (Equation 5.0.2)
        rip = 1.0 - ((1.0 - ri) ** cycle_length)

        return rip

    def calculate_mortality_distribution(self, basal_area_percentile: float) -> float:
        """Calculate mortality distribution factor (MR).

        Implements equation 5.0.3:
        MR = 0.84525 - 0.01074*PCT + 0.0000002*PCT³

        Args:
            basal_area_percentile: Tree's position in BA distribution (0-100)

        Returns:
            Mortality distribution factor (bounded 0.01-1.0)
        """
        pct = basal_area_percentile
        mr = 0.84525 - (0.01074 * pct) + (0.0000002 * (pct ** 3))
        return max(0.01, min(1.0, mr))


class LSMortalityModel:
    """FVS Lake States mortality model.

    Implements the LS variant mortality model from morts.f:
    1. 4-group background mortality with halved logistic rates
    2. SDI-based density mortality with VARADJ shade tolerance adjustment

    The 4 mortality groups have different logistic coefficients (PMSC, PMD)
    from morts.f:
    - Group 1: PMSC=5.1677, PMD=-0.00777
    - Group 2: PMSC=9.6943, PMD=-0.01273
    - Group 3: PMSC=5.5877, PMD=-0.00535
    - Group 4: PMSC=5.9617, PMD=-0.03401

    Background rates are halved compared to the logistic.
    Density mortality selects EITHER background OR density-based (never both),
    matching Fortran varmrt.f behavior.

    Subclasses (e.g., NEMortalityModel) override class attributes to use
    variant-specific coefficient files, species-group mappings, and SDI maximums.
    """

    # Coefficient file path (relative to cfg/), overridden by subclasses
    _COEFFICIENT_FILE = 'ls/ls_mortality_coefficients.json'

    # SDI threshold constants (same as SN)
    LOWER_THRESHOLD = 0.55
    UPPER_THRESHOLD = 0.85

    # 4-group background mortality coefficients (PMSC, PMD) from morts.f lines 99-100
    # PMSC = [5.1677, 9.6943, 5.5877, 5.9617]
    # PMD  = [-0.00777, -0.01273, -0.00535, -0.03401]
    MORTALITY_GROUPS = {
        1: {'P0': 5.1677, 'P1': -0.00777},
        2: {'P0': 9.6943, 'P1': -0.01273},
        3: {'P0': 5.5877, 'P1': -0.00535},
        4: {'P0': 5.9617, 'P1': -0.03401},
    }

    # Species-to-mortality-group mapping from IMAPLS in morts.f (68 species)
    SPECIES_MORTALITY_GROUP = {
        # Group 1 (28 species)
        'JP': 1, 'RN': 1, 'RP': 1, 'WS': 1, 'NS': 1, 'BF': 1,
        'BS': 1, 'TA': 1, 'WC': 1, 'EH': 1, 'GA': 1, 'SV': 1,
        'RM': 1, 'AE': 1, 'RL': 1, 'RE': 1, 'BW': 1, 'SM': 1,
        'BM': 1, 'AB': 1, 'HH': 1, 'BK': 1, 'AH': 1, 'DW': 1,
        'BG': 1, 'WI': 1, 'BL': 1, 'SS': 1,
        # Group 3 (4 species)
        'SC': 3, 'WP': 3, 'OS': 3, 'RC': 3,
        # Group 4 (35 species)
        'BA': 4, 'EC': 4, 'BC': 4, 'YB': 4, 'WA': 4, 'WO': 4,
        'SW': 4, 'BR': 4, 'CK': 4, 'RO': 4, 'BO': 4, 'NP': 4,
        'BH': 4, 'PH': 4, 'SH': 4, 'BT': 4, 'QA': 4, 'BP': 4,
        'PB': 4, 'BN': 4, 'WN': 4, 'OH': 4, 'BE': 4, 'ST': 4,
        'MM': 4, 'AC': 4, 'HK': 4, 'HT': 4, 'AP': 4, 'SY': 4,
        'PR': 4, 'CC': 4, 'PL': 4, 'DM': 4, 'MA': 4,
    }

    # Default SDI maximums (overridden by subclasses)
    _SDI_MAXIMUMS = {
        'JP': 400, 'SC': 400, 'RN': 500, 'RP': 500, 'WP': 450,
        'WS': 500, 'NS': 500, 'BF': 400, 'BS': 400, 'TA': 350,
        'WC': 400, 'EH': 500, 'SM': 450, 'RM': 400, 'QA': 350,
        'PB': 350, 'RO': 400, 'WO': 400, 'YB': 400, 'AB': 450,
    }
    DEFAULT_SDI_MAX = 400

    # Fallback shade tolerances (overridden by subclasses)
    _FALLBACK_SHADE_TOLERANCE = {
        'JP': 0.30, 'RN': 0.30, 'WP': 0.50, 'BF': 0.90,
        'SM': 0.90, 'RM': 0.85, 'QA': 0.10, 'PB': 0.30,
        'RO': 0.50, 'WO': 0.50,
    }

    def __init__(self, default_species: str = 'RN', max_sdi: Optional[float] = None,
                 variant: str = 'LS'):
        """Initialize the LS mortality model.

        Args:
            default_species: Default species code for coefficient lookups
            max_sdi: Maximum SDI for the stand (if None, uses species default)
            variant: Variant code
        """
        self.default_species = default_species
        self.max_sdi = max_sdi
        self.variant = variant
        self._shade_tolerance = {}
        self._load_shade_tolerance()

    def _load_shade_tolerance(self):
        """Load shade tolerance (VARADJ) from mortality coefficients."""
        try:
            from .config_loader import load_coefficient_file
            data = load_coefficient_file(self._COEFFICIENT_FILE)
            coeffs = data.get('coefficients', {})
            for species, values in coeffs.items():
                self._shade_tolerance[species] = values.get('shade_tolerance', 0.30)
        except (FileNotFoundError, KeyError):
            self._shade_tolerance = dict(self._FALLBACK_SHADE_TOLERANCE)

    def _get_background_rate(self, tree: 'Tree', cycle_length: int) -> float:
        """Calculate LS 4-group background mortality rate.

        Rate is halved compared to the base logistic equation per morts.f.

        Args:
            tree: Tree object
            cycle_length: Cycle length in years

        Returns:
            Mortality probability (0-1)
        """
        group = self.SPECIES_MORTALITY_GROUP.get(tree.species, 3)
        group_coeffs = self.MORTALITY_GROUPS[group]
        p0 = group_coeffs['P0']
        p1 = group_coeffs['P1']

        # Base annual mortality rate (logistic)
        ri = 1.0 / (1.0 + math.exp(p0 + p1 * tree.dbh))

        # LS halves the background rate (slower mortality than SN)
        ri *= 0.5

        # Adjust for cycle length
        rip = 1.0 - ((1.0 - ri) ** cycle_length)

        return rip

    def apply_mortality(
        self,
        trees: List['Tree'],
        cycle_length: int = 5,
        max_sdi: Optional[float] = None,
        random_seed: Optional[int] = None
    ) -> MortalityResult:
        """Apply LS mortality model with 4-group background and SDI density.

        Args:
            trees: List of trees in the stand
            cycle_length: Length of projection cycle in years
            max_sdi: Maximum SDI (uses species default if None)
            random_seed: Optional seed for reproducibility

        Returns:
            MortalityResult with survivors and mortality count
        """
        if len(trees) <= 1:
            return MortalityResult(survivors=list(trees), mortality_count=0, trees_died=[])

        if random_seed is not None:
            random.seed(random_seed)

        max_sdi = max_sdi or self.max_sdi or self._SDI_MAXIMUMS.get(
            self.default_species, self.DEFAULT_SDI_MAX
        )

        # Calculate stand SDI
        current_sdi = self._calculate_stand_sdi(trees)
        relative_sdi = current_sdi / max_sdi

        # Calculate basal area percentiles
        tree_data = self._calculate_tree_percentiles(trees)

        survivors = []
        trees_died = []

        # Calculate density mortality fraction if above threshold
        density_removal_fraction = 0.0
        if relative_sdi > self.LOWER_THRESHOLD:
            if relative_sdi > self.UPPER_THRESHOLD:
                target_sdi = self.UPPER_THRESHOLD * max_sdi
                excess_sdi = current_sdi - target_sdi
                sdi_to_remove = excess_sdi * 0.5
            else:
                target_sdi = current_sdi * (
                    1.0 - 0.05 * (relative_sdi - self.LOWER_THRESHOLD) /
                    (self.UPPER_THRESHOLD - self.LOWER_THRESHOLD)
                )
                sdi_to_remove = current_sdi - target_sdi
            density_removal_fraction = sdi_to_remove / current_sdi if current_sdi > 0 else 0

        for tree, pct in tree_data:
            # Background mortality (LS 4-group, halved)
            rip = self._get_background_rate(tree, cycle_length)

            # Density-dependent mortality component
            if density_removal_fraction > 0:
                # Mortality distribution by size (same equation as SN)
                mr = 0.84525 - (0.01074 * pct) + (0.0000002 * (pct ** 3))
                mr = max(0.01, min(1.0, mr))

                # VARADJ shade tolerance adjustment from varmrt.f:
                # EFFTR(I) = PEFF * (1.0 - VARADJ(JSPC)) * 0.1
                varadj = self._shade_tolerance.get(tree.species, 0.30)
                shade_modifier = (1.0 - varadj) * 0.1

                mort = mr * shade_modifier * density_removal_fraction * (cycle_length / 5.0)
            else:
                mort = 0.0

            # Select mortality rate (Fortran uses one or the other, never both)
            # varmrt.f: IF(T .LE. TEM .OR. RN .LE. 0.0) RIP = RI
            total_mort_prob = mort if density_removal_fraction > 0 else rip

            if random.random() > total_mort_prob:
                survivors.append(tree)
            else:
                trees_died.append(tree)

        return MortalityResult(
            survivors=survivors,
            mortality_count=len(trees_died),
            trees_died=trees_died
        )

    def _calculate_stand_sdi(self, trees: List['Tree']) -> float:
        """Calculate stand SDI using Reineke's equation."""
        if not trees:
            return 0.0
        tpa = len(trees)
        qmd_squared = sum(tree.dbh ** 2 for tree in trees) / tpa
        qmd = math.sqrt(qmd_squared)
        return tpa * (qmd / 10.0) ** 1.605

    def _calculate_tree_percentiles(
        self, trees: List['Tree']
    ) -> List[Tuple['Tree', float]]:
        """Calculate basal area percentile for each tree."""
        total_ba = calculate_stand_basal_area(trees)
        tree_data = []
        cumulative_ba = 0.0
        sorted_trees = sorted(trees, key=lambda t: t.dbh)
        for tree in sorted_trees:
            tree_ba = calculate_tree_basal_area(tree.dbh)
            cumulative_ba += tree_ba
            pct = (cumulative_ba / total_ba) * 100.0 if total_ba > 0 else 50.0
            tree_data.append((tree, pct))
        return tree_data

    def calculate_background_mortality_rate(
        self, tree: 'Tree', cycle_length: int = 5
    ) -> float:
        """Calculate background mortality rate for a single tree."""
        return self._get_background_rate(tree, cycle_length)


class NEMortalityModel(LSMortalityModel):
    """FVS Northeast mortality model.

    NE uses the same 4-group background + SDI density model as LS (TWIGS family),
    with NE-specific species-group mappings, shade tolerances, and SDI maximums.
    """

    _COEFFICIENT_FILE = 'ne/ne_mortality_coefficients.json'

    # NE species-to-mortality-group mapping from IMAPNE in morts.f (108 species)
    SPECIES_MORTALITY_GROUP = {
        # Group 1 (38 species)
        'BF': 1, 'TA': 1, 'WS': 1, 'RS': 1, 'NS': 1, 'BS': 1,
        'PI': 1, 'RN': 1, 'WC': 1, 'AW': 1, 'EH': 1, 'HM': 1,
        'JP': 1, 'RM': 1, 'SM': 1, 'BM': 1, 'SV': 1, 'SB': 1,
        'AB': 1, 'AS': 1, 'GA': 1, 'PA': 1, 'BU': 1, 'YY': 1,
        'MG': 1, 'BG': 1, 'SD': 1, 'BK': 1, 'BL': 1, 'SS': 1,
        'BW': 1, 'WB': 1, 'EL': 1, 'AE': 1, 'RL': 1, 'AH': 1,
        'DW': 1, 'HH': 1,
        # Group 2 (1 species)
        'OS': 2,
        # Group 3 (11 species)
        'WP': 3, 'LP': 3, 'VP': 3, 'RC': 3, 'JU': 3, 'OP': 3,
        'SP': 3, 'TM': 3, 'PP': 3, 'PD': 3, 'SC': 3,
        # Group 4 (57 species)
        'YB': 4, 'RB': 4, 'PB': 4, 'GB': 4, 'HI': 4, 'PH': 4,
        'SL': 4, 'SH': 4, 'MH': 4, 'WA': 4, 'BA': 4, 'YP': 4,
        'SU': 4, 'CT': 4, 'QA': 4, 'BP': 4, 'EC': 4, 'BT': 4,
        'PY': 4, 'BC': 4, 'WO': 4, 'BR': 4, 'CK': 4, 'PO': 4,
        'OK': 4, 'SO': 4, 'QI': 4, 'WK': 4, 'PN': 4, 'CO': 4,
        'SW': 4, 'SN': 4, 'RO': 4, 'SK': 4, 'BO': 4, 'CB': 4,
        'WR': 4, 'HK': 4, 'PS': 4, 'HY': 4, 'BN': 4, 'WN': 4,
        'OO': 4, 'MV': 4, 'AP': 4, 'WT': 4, 'PW': 4, 'SY': 4,
        'WL': 4, 'OH': 4, 'BE': 4, 'ST': 4, 'AI': 4, 'SE': 4,
        'HT': 4, 'PL': 4, 'PR': 4,
    }

    # NE SDI maximums (references StandMetricsCalculator.NE_SDI_MAXIMUMS)
    _SDI_MAXIMUMS = {
        'BF': 400, 'TA': 350, 'WS': 500, 'RS': 450, 'NS': 450,
        'BS': 400, 'PI': 400, 'RN': 500, 'WP': 450, 'LP': 450,
        'VP': 350, 'OP': 400, 'JP': 400, 'SP': 400, 'TM': 350,
        'PP': 350, 'PD': 350, 'SC': 400, 'WC': 400, 'AW': 350,
        'RC': 350, 'JU': 350, 'EH': 500, 'HM': 450, 'OS': 400,
        'RM': 400, 'SM': 450, 'BM': 450, 'SV': 400, 'WO': 400,
        'RO': 400, 'YB': 400, 'AB': 450, 'WA': 350, 'BC': 400,
        'YP': 450, 'WN': 400, 'QA': 350, 'PB': 350,
    }

    _FALLBACK_SHADE_TOLERANCE = {
        'RM': 0.85, 'SM': 0.90, 'WP': 0.50, 'RO': 0.50,
        'WO': 0.50, 'BF': 0.90, 'RS': 0.70, 'EH': 0.90,
        'YB': 0.50, 'AB': 0.90, 'WA': 0.30, 'BC': 0.40,
    }

    def __init__(self, default_species: str = 'RM', max_sdi: Optional[float] = None,
                 variant: str = 'NE'):
        """Initialize the NE mortality model.

        Args:
            default_species: Default species code for coefficient lookups
            max_sdi: Maximum SDI for the stand (if None, uses species default)
            variant: Variant code
        """
        super().__init__(default_species=default_species, max_sdi=max_sdi, variant=variant)


class CSMortalityModel(LSMortalityModel):
    """FVS Central States mortality model.

    CS uses the same 4-group background + SDI density model as LS/NE (TWIGS family),
    with CS-specific species-group mappings, shade tolerances, and SDI maximums.
    """

    _COEFFICIENT_FILE = 'cs/cs_mortality_coefficients.json'

    # CS species-to-mortality-group mapping from IMAPCS in morts.f (96 species)
    SPECIES_MORTALITY_GROUP = {
        # Group 1 (30 species)
        'TL': 1, 'BG': 1, 'AB': 1, 'PA': 1, 'UA': 1, 'RM': 1,
        'SV': 1, 'AE': 1, 'WE': 1, 'EL': 1, 'SI': 1, 'RL': 1,
        'RE': 1, 'BW': 1, 'SM': 1, 'AS': 1, 'GA': 1, 'SS': 1,
        'OB': 1, 'BK': 1, 'WI': 1, 'BL': 1, 'AH': 1, 'RD': 1,
        'DW': 1, 'KC': 1, 'OO': 1, 'MB': 1, 'HH': 1, 'SD': 1,
        # Group 2 (2 species)
        'JU': 2, 'OS': 2,
        # Group 3 (6 species)
        'RC': 3, 'SP': 3, 'VP': 3, 'LP': 3, 'WP': 3, 'BY': 3,
        # Group 4 (55 species)
        'WN': 4, 'BN': 4, 'TS': 4, 'WT': 4, 'SH': 4, 'SL': 4,
        'MH': 4, 'PH': 4, 'HI': 4, 'WH': 4, 'BH': 4, 'PE': 4,
        'BI': 4, 'BA': 4, 'EC': 4, 'BE': 4, 'BC': 4, 'SG': 4,
        'HK': 4, 'YP': 4, 'WA': 4, 'WO': 4, 'RO': 4, 'SK': 4,
        'BO': 4, 'SO': 4, 'BJ': 4, 'CK': 4, 'SW': 4, 'BR': 4,
        'SN': 4, 'PO': 4, 'DO': 4, 'CO': 4, 'PN': 4, 'CB': 4,
        'QI': 4, 'OV': 4, 'WK': 4, 'NK': 4, 'WL': 4, 'QS': 4,
        'CA': 4, 'PS': 4, 'HL': 4, 'BP': 4, 'BT': 4, 'QA': 4,
        'SY': 4, 'RB': 4, 'SU': 4, 'OH': 4, 'CT': 4, 'MV': 4,
        'HT': 4,
    }

    # CS SDI maximums (references StandMetricsCalculator.CS_SDI_MAXIMUMS)
    _SDI_MAXIMUMS = {
        'RC': 400, 'JU': 350, 'SP': 450, 'VP': 400, 'LP': 450,
        'OS': 400, 'WP': 450, 'BY': 400,
        'WN': 400, 'BN': 350, 'TL': 400, 'TS': 400,
        'HS': 380, 'SH': 380, 'SL': 380, 'MH': 380, 'PH': 380,
        'HI': 380, 'WH': 380, 'BH': 380, 'PE': 380, 'BI': 380,
        'AB': 450, 'BA': 350, 'PA': 350, 'UA': 350,
        'RM': 400, 'SM': 450, 'BE': 400, 'SV': 400,
        'BC': 400, 'AE': 400, 'SG': 400, 'WE': 350, 'EL': 350,
        'SI': 350, 'RL': 350, 'RE': 350,
        'YP': 450, 'BW': 350, 'WA': 400, 'GA': 400, 'AS': 400,
        'WO': 420, 'RO': 420, 'SK': 400, 'BO': 420, 'SO': 400,
        'BJ': 350, 'CK': 420, 'SW': 420, 'BR': 420, 'SN': 400,
        'PO': 380, 'DO': 380, 'CO': 420, 'PN': 400, 'CB': 420,
        'QI': 380, 'OV': 380, 'WK': 380, 'NK': 400, 'WL': 400,
        'QS': 400,
    }

    _FALLBACK_SHADE_TOLERANCE = {
        'WO': 0.50, 'RO': 0.50, 'SM': 0.90, 'WN': 0.30,
        'YP': 0.30, 'WA': 0.30, 'BC': 0.40, 'SP': 0.30,
        'RM': 0.85, 'AB': 0.90, 'BO': 0.50, 'SH': 0.50,
    }

    def __init__(self, default_species: str = 'WO', max_sdi: Optional[float] = None,
                 variant: str = 'CS'):
        """Initialize the CS mortality model.

        Args:
            default_species: Default species code for coefficient lookups
            max_sdi: Maximum SDI for the stand (if None, uses species default)
            variant: Variant code
        """
        super().__init__(default_species=default_species, max_sdi=max_sdi, variant=variant)


# Module-level convenience functions
_default_model: Optional[MortalityModel] = None


def create_mortality_model(
    default_species: str = 'LP',
    max_sdi: Optional[float] = None,
    variant: Optional[str] = None
):
    """Factory function to create a variant-appropriate mortality model.

    Args:
        default_species: Default species code
        max_sdi: Maximum SDI for the stand
        variant: FVS variant code (e.g., 'SN', 'LS')

    Returns:
        MortalityModel or LSMortalityModel instance
    """
    if variant is None:
        from .config_loader import get_default_variant
        variant = get_default_variant()

    if variant == 'LS':
        return LSMortalityModel(default_species=default_species, max_sdi=max_sdi)
    elif variant == 'NE':
        return NEMortalityModel(default_species=default_species, max_sdi=max_sdi)
    elif variant == 'CS':
        return CSMortalityModel(default_species=default_species, max_sdi=max_sdi)
    elif variant == 'PN':
        # PN uses the same base FVS SDI mortality model as SN
        # but with PN-specific SDI maximums (loaded via StandMetricsCalculator)
        from .stand_metrics import StandMetricsCalculator
        if max_sdi is None:
            pn_sdi_maxs = StandMetricsCalculator.PN_SDI_MAXIMUMS
            max_sdi = pn_sdi_maxs.get(default_species, 850)
        return MortalityModel(default_species=default_species, max_sdi=max_sdi)
    elif variant == 'WC':
        # WC uses the same SDI mortality model with WC-specific SDI maximums
        from .stand_metrics import StandMetricsCalculator
        if max_sdi is None:
            wc_sdi_maxs = StandMetricsCalculator.WC_SDI_MAXIMUMS
            max_sdi = wc_sdi_maxs.get(default_species, 800)
        return MortalityModel(default_species=default_species, max_sdi=max_sdi)
    elif variant == 'OP':
        # OP uses the same SDI mortality model with OP-specific SDI maximums
        from .stand_metrics import StandMetricsCalculator
        if max_sdi is None:
            op_sdi_maxs = StandMetricsCalculator.OP_SDI_MAXIMUMS
            max_sdi = op_sdi_maxs.get(default_species, 850)
        return MortalityModel(default_species=default_species, max_sdi=max_sdi)
    else:
        return MortalityModel(default_species=default_species, max_sdi=max_sdi)


def get_mortality_model(species: str = 'LP') -> MortalityModel:
    """Get or create a mortality model instance.

    Args:
        species: Default species code

    Returns:
        MortalityModel instance
    """
    global _default_model
    if _default_model is None:
        _default_model = MortalityModel(species)
    return _default_model


def apply_stand_mortality(
    trees: List['Tree'],
    cycle_length: int = 5,
    max_sdi: Optional[float] = None
) -> MortalityResult:
    """Apply mortality to a list of trees.

    Convenience function using the default model.

    Args:
        trees: List of trees
        cycle_length: Cycle length in years
        max_sdi: Maximum SDI (uses default if None)

    Returns:
        MortalityResult
    """
    return get_mortality_model().apply_mortality(
        trees, cycle_length=cycle_length, max_sdi=max_sdi
    )
