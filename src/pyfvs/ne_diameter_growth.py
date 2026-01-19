"""
Northeast (NE) variant diameter growth model.

Implements the NE-TWIGS diameter growth equation for the FVS Northeast variant.
The NE variant uses a basal area increment model:

NE-TWIGS equation:
    Potential BA Growth = B1 * SI * (1 - exp(-B2 * DBH))
    Adjusted Growth = POTBAG * 0.7 * BAGMOD

where:
    - B1, B2 are species-specific coefficients
    - SI = site index
    - DBH = diameter at breast height
    - BAGMOD = competition modifier from basal area in larger trees (BAL)

The growth is then converted from basal area to diameter increment.
This model iterates 10 times annually per the dgf.f source code.

Key differences from other variants:
    - Uses a potential basal area growth approach
    - Simple 2-coefficient species equations (B1, B2)
    - 0.7 modifier applied to potential growth
    - Competition through BAL-based modifier
    - 10-year base cycle

Source: FVS Northeast Variant dgf.f, USDA Forest Service
"""
import math
from typing import Dict, Any, Optional

from .model_base import ParameterizedModel
from .config_loader import load_coefficient_file


class NEDiameterGrowthModel(ParameterizedModel):
    """Northeast variant diameter growth model.

    Calculates diameter increment using the NE-TWIGS approach.

    The NE variant uses potential basal area growth:
        POTBAG = B1 * SI * (1 - exp(-B2 * DBH))
        Adjusted = POTBAG * 0.7 * BAGMOD

    Then converts to diameter increment.

    Attributes:
        species_code: Species code (e.g., 'RM', 'SM', 'WP')
        coefficients: Species-specific coefficients (B1, B2)
    """

    COEFFICIENT_FILE = 'ne/ne_diameter_growth_coefficients.json'
    COEFFICIENT_KEY = 'coefficients'
    DEFAULT_SPECIES = 'RM'  # Red Maple is the default site species for NE

    # Fallback parameters for key NE species (from dgf.f)
    FALLBACK_PARAMETERS = {
        'RM': {  # Red Maple
            'B1': 0.0007906,
            'B2': 0.0651982,
        },
        'SM': {  # Sugar Maple
            'B1': 0.0007439,
            'B2': 0.0706905,
        },
        'WP': {  # Eastern White Pine
            'B1': 0.0011303,
            'B2': 0.0934796,
        },
        'RO': {  # Northern Red Oak
            'B1': 0.0008920,
            'B2': 0.0979702,
        },
        'YB': {  # Yellow Birch
            'B1': 0.0006668,
            'B2': 0.0768212,
        },
        'AB': {  # American Beech
            'B1': 0.0006911,
            'B2': 0.0730441,
        },
        'EH': {  # Eastern Hemlock
            'B1': 0.0008737,
            'B2': 0.0940538,
        },
        'BF': {  # Balsam Fir
            'B1': 0.0008829,
            'B2': 0.0602785,
        },
        'RS': {  # Red Spruce
            'B1': 0.0008236,
            'B2': 0.0549439,
        },
        'WA': {  # White Ash
            'B1': 0.0008992,
            'B2': 0.0925395,
        },
    }

    def __init__(self, species_code: str = None, variant: str = 'NE'):
        """Initialize the NE diameter growth model.

        Args:
            species_code: Species code (e.g., 'RM', 'SM', 'WP').
                         Defaults to DEFAULT_SPECIES (RM).
            variant: FVS variant (should be 'NE' for this model)
        """
        self.variant = variant
        super().__init__(species_code)

    def _get_coefficient_data(self) -> Dict[str, Any]:
        """Load coefficient data from JSON file using ConfigLoader with caching.

        Returns:
            Dictionary containing the full coefficient file data.
        """
        try:
            return load_coefficient_file(self.COEFFICIENT_FILE, variant='NE')
        except FileNotFoundError:
            return {}

    def _calculate_bagmod(self, bal: float, ba: float) -> float:
        """Calculate basal area growth modifier based on competition.

        The BAGMOD reduces growth based on basal area in larger trees (BAL).
        This is a simplified version of the BALMOD subroutine in dgf.f.

        Args:
            bal: Basal area in larger trees (sq ft/acre)
            ba: Total stand basal area (sq ft/acre)

        Returns:
            BAGMOD: Growth modifier (0-1)
        """
        # Simple competition modifier based on relative position
        # Trees with more BAL (more overtopped) grow slower
        if ba <= 0:
            return 1.0

        # Competition ratio: 0 = dominant, 1 = fully suppressed
        competition_ratio = min(1.0, bal / max(1.0, ba))

        # Modifier decreases as competition increases
        # This follows the general FVS approach of reducing growth
        # for suppressed trees
        bagmod = 1.0 - (0.5 * competition_ratio)

        return max(0.1, min(1.0, bagmod))

    def calculate_dds(
        self,
        dbh: float,
        crown_ratio: float,
        site_index: float,
        ba: float,
        bal: float,
        time_step: float = 10.0
    ) -> float:
        """Calculate diameter squared increment (DDS).

        The NE variant uses an iterative basal area growth approach:
            POTBAG = B1 * SI * (1 - exp(-B2 * DBH))
            Growth = POTBAG * 0.7 * BAGMOD

        Per the FVS dgf.f source, this iterates once per year in the cycle,
        with diameter updating each iteration. This accumulates growth over
        the time step.

        Args:
            dbh: Diameter at breast height (inches)
            crown_ratio: Crown ratio as proportion (0-1)
            site_index: Site index (base age 50 for NE) in feet
            ba: Stand basal area (sq ft/acre)
            bal: Basal area in larger trees (sq ft/acre)
            time_step: Growth period in years (default 10 for NE)

        Returns:
            DDS: Change in diameter squared (sq inches inside bark)
        """
        p = self.coefficients

        # Get coefficients
        b1 = p.get('B1', 0.0008)
        b2 = p.get('B2', 0.08)

        # Ensure valid inputs
        current_dbh = max(0.1, dbh)
        si_safe = max(20.0, site_index)
        cr_safe = max(0.1, min(1.0, crown_ratio))

        # Calculate competition modifier (constant for the period)
        bagmod = self._calculate_bagmod(bal, ba)

        # Crown ratio effect (trees with lower CR grow slower)
        cr_modifier = 0.5 + (0.5 * cr_safe)  # Range 0.55 to 1.0

        # Number of annual iterations
        num_years = int(time_step)

        # Iterate once per year, updating diameter each year
        # This matches the FVS dgf.f DO 1000 ILOOP=1,10 loop
        for _ in range(num_years):
            # Calculate potential basal area growth for this year
            # POTBAG = B1 * SI * (1 - exp(-B2 * DBH))
            potbag = b1 * si_safe * (1.0 - math.exp(-b2 * current_dbh))

            # Apply 0.7 modifier (as in dgf.f)
            potbag = potbag * 0.7

            # Calculate adjusted annual basal area growth
            annual_ba_growth = potbag * bagmod * cr_modifier

            # Current tree basal area (sq ft)
            # BA = pi/4 * D^2 / 144
            current_ba = (math.pi / 4.0) * (current_dbh ** 2) / 144.0

            # New basal area
            new_ba = current_ba + annual_ba_growth

            # Convert back to diameter
            # D = sqrt(BA * 144 * 4 / pi)
            current_dbh = math.sqrt(new_ba * 144.0 * 4.0 / math.pi)

        # DDS is the change in diameter squared over the whole period
        dds = max(0.0, current_dbh ** 2 - dbh ** 2)

        # Apply reasonable bounds
        # Max DDS ~ 50 sq in corresponds to ~3.5" growth per decade for 10" tree
        dds = min(50.0, dds)

        return dds

    def calculate_diameter_growth(
        self,
        dbh: float,
        crown_ratio: float,
        site_index: float,
        ba: float,
        bal: float,
        bark_ratio: float = 0.9,
        time_step: float = 10.0
    ) -> float:
        """Calculate diameter growth in inches.

        Converts DDS to actual diameter increment, applying bark ratio
        to work with inside-bark measurements like FVS.

        Args:
            dbh: Current DBH (outside bark) in inches
            crown_ratio: Crown ratio as proportion (0-1)
            site_index: Site index (base age 50) in feet
            ba: Stand basal area (sq ft/acre)
            bal: Basal area in larger trees (sq ft/acre)
            bark_ratio: DIB/DOB ratio (default 0.9)
            time_step: Growth period in years (default 10)

        Returns:
            Diameter increment in inches (outside bark)
        """
        # Calculate DDS
        dds = self.calculate_dds(
            dbh=dbh,
            crown_ratio=crown_ratio,
            site_index=site_index,
            ba=ba,
            bal=bal,
            time_step=time_step
        )

        # Convert to inside-bark diameter
        dib_old = dbh * bark_ratio
        dib_old_sq = dib_old * dib_old

        # Apply DDS to inside-bark diameter
        dib_new = math.sqrt(dib_old_sq + dds)

        # Convert back to outside-bark diameter
        dbh_new = dib_new / bark_ratio

        # Return the increment
        return max(0.0, dbh_new - dbh)


# Module-level cache for model instances
_model_cache: Dict[str, NEDiameterGrowthModel] = {}


def create_ne_diameter_growth_model(species_code: str = 'RM') -> NEDiameterGrowthModel:
    """Factory function to create a cached NE diameter growth model.

    Args:
        species_code: Species code (e.g., 'RM', 'SM', 'WP')

    Returns:
        Cached NEDiameterGrowthModel instance
    """
    species_upper = species_code.upper()
    if species_upper not in _model_cache:
        _model_cache[species_upper] = NEDiameterGrowthModel(species_upper)
    return _model_cache[species_upper]


# Backwards compatibility alias
get_ne_diameter_growth_model = create_ne_diameter_growth_model


def calculate_ne_diameter_growth(
    dbh: float,
    crown_ratio: float,
    site_index: float,
    ba: float,
    bal: float,
    species_code: str = 'RM',
    bark_ratio: float = 0.9,
    time_step: float = 10.0
) -> float:
    """Convenience function to calculate NE diameter growth.

    Args:
        dbh: Current DBH (outside bark) in inches
        crown_ratio: Crown ratio as proportion (0-1)
        site_index: Site index (base age 50) in feet
        ba: Stand basal area (sq ft/acre)
        bal: Basal area in larger trees (sq ft/acre)
        species_code: Species code (default 'RM' - Red Maple)
        bark_ratio: DIB/DOB ratio (default 0.9)
        time_step: Growth period in years (default 10)

    Returns:
        Diameter increment in inches (outside bark)
    """
    model = get_ne_diameter_growth_model(species_code)
    return model.calculate_diameter_growth(
        dbh=dbh,
        crown_ratio=crown_ratio,
        site_index=site_index,
        ba=ba,
        bal=bal,
        bark_ratio=bark_ratio,
        time_step=time_step
    )
