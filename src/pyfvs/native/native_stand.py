"""
High-level NativeStand API that mirrors the pyfvs Stand interface.

NativeStand provides a familiar API for running simulations through the
native FVS Fortran library, producing results in the same format as
pyfvs.Stand for direct comparison.

Usage:
    >>> from pyfvs.native import NativeStand
    >>> with NativeStand(variant='SN') as ns:
    ...     ns.initialize_planted(500, 70, 'LP')
    ...     ns.grow(50)
    ...     metrics = ns.get_metrics()
    ...     print(f"TPA: {metrics['tpa']:.0f}, BA: {metrics['basal_area']:.1f}")

Important limitations:
    - FVS uses Fortran COMMON blocks for global state. Only one NativeStand
      instance per variant can be active at a time within a process.
    - Requires the FVS shared library to be built and installed.
      See native/BUILD.md for instructions.
"""

import logging
import os
import shutil
import tempfile
from pathlib import Path
from typing import Dict, List, Optional

from pyfvs.exceptions import FVSNativeError
from pyfvs.native.fvs_bindings import FVSBindings
from pyfvs.native.species_map import get_species_index

logger = logging.getLogger(__name__)


class NativeStand:
    """High-level interface to the native FVS shared library.

    Mirrors the pyfvs Stand API for direct comparison of simulation results.
    Generates FVS keyword files and manages temporary working directories.

    Args:
        variant: FVS variant code (e.g., 'SN', 'PN', 'LS').
        work_dir: Optional working directory. If None, a temporary
            directory is created and cleaned up on close().
    """

    def __init__(self, variant: str = "SN", work_dir: Optional[str] = None):
        self.variant = variant.upper()
        self._bindings = FVSBindings(self.variant)
        self._initialized = False
        self._years_grown = 0
        self._cycle_length = 5 if self.variant in ("SN", "OP") else 10

        # Working directory for keyword files
        if work_dir:
            self._work_dir = Path(work_dir)
            self._work_dir.mkdir(parents=True, exist_ok=True)
            self._owns_work_dir = False
        else:
            self._work_dir = Path(tempfile.mkdtemp(prefix="fvs_native_"))
            self._owns_work_dir = True

        self._keyword_file: Optional[Path] = None
        self._species: Optional[str] = None
        self._site_index: Optional[float] = None
        self._trees_per_acre: Optional[float] = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    def close(self):
        """Clean up temporary working directory."""
        if self._owns_work_dir and self._work_dir.exists():
            shutil.rmtree(self._work_dir, ignore_errors=True)

    # =========================================================================
    # Initialization
    # =========================================================================

    def initialize_planted(
        self,
        trees_per_acre: float,
        site_index: float,
        species: str,
        age: int = 1,
        dbh: float = 0.1,
        height: float = 1.0,
        ecounit: str = "",
    ) -> None:
        """Initialize a planted stand, matching the pyfvs Stand API.

        Generates an FVS keyword file with STDIDENT, SITECODE, DESIGN,
        PLANT/TREEFMT, NUMCYCLE, and PROCESS keywords, then initializes
        the native FVS library.

        Args:
            trees_per_acre: Planting density (trees per acre).
            site_index: Site index (base age 25 or 50 depending on variant).
            species: Two-letter species code (e.g., 'LP', 'DF').
            age: Initial stand age in years.
            dbh: Initial DBH in inches.
            height: Initial height in feet.
            ecounit: Ecological unit code (if applicable).
        """
        self._species = species.upper()
        self._site_index = site_index
        self._trees_per_acre = trees_per_acre

        # Map species to FVS integer index
        species_index = get_species_index(self._species, self.variant)

        # Calculate number of cycles needed (will be set during grow())
        num_cycles = 10  # default, updated in grow()

        # Generate keyword file
        self._keyword_file = self._generate_keyword_file(
            species_index=species_index,
            trees_per_acre=trees_per_acre,
            site_index=site_index,
            age=age,
            dbh=dbh,
            height=height,
            num_cycles=num_cycles,
            ecounit=ecounit,
        )

        logger.info(
            "NativeStand initialized: %s variant, %s species, "
            "SI=%.0f, TPA=%.0f, keyword=%s",
            self.variant,
            self._species,
            site_index,
            trees_per_acre,
            self._keyword_file,
        )

        self._initialized = True

    def _generate_keyword_file(
        self,
        species_index: int,
        trees_per_acre: float,
        site_index: float,
        age: int,
        dbh: float,
        height: float,
        num_cycles: int,
        ecounit: str = "",
    ) -> Path:
        """Generate an FVS keyword file for the simulation.

        Returns:
            Path to the generated .key file.
        """
        key_path = self._work_dir / "stand.key"
        tree_path = self._work_dir / "stand.tre"

        # Write tree data file (FVS free-format tree list)
        # Format: Plot Tree Species DBH Height TPA CrownRatio
        with open(tree_path, "w") as tree_file:
            tree_file.write(
                f"   1    1  {species_index:3d}  {dbh:6.1f}  {height:6.1f}"
                f"  {trees_per_acre:8.1f}    0\n"
            )

        # Build keyword file
        lines = []
        lines.append(f"STDIDENT")
        lines.append(f"PYFVS_NATIVE_{self.variant}_{self._species}")

        # Site index
        lines.append(f"SITECODE")
        lines.append(f"         {species_index:3d}  {site_index:6.1f}")

        # Design (1 plot, no sampling weight adjustments)
        lines.append(f"DESIGN")
        lines.append(f"       0       1       1     1.0     1.0     1.0")

        # Stand age
        lines.append(f"STDINFO")
        lines.append(f"         0  {age:4d}       0       0       0")

        # Number of cycles
        lines.append(f"NUMCYCLE")
        lines.append(f"  {num_cycles:8d}")

        # Time interval (cycle length)
        lines.append(f"TIMEINT")
        lines.append(f"       0  {self._cycle_length:8d}")

        # Tree data input
        lines.append(f"TREEFMT")
        lines.append(f"(I4,T1,I7,F3.0,I1,A3,F4.1,F3.1,2F3.0,F4.1,I1,3F2.0,2I1,I2,2I3,2I1,F3.0)")
        lines.append(f"TREEIN")
        lines.append(f"{tree_path}")

        # FVS PROCESS keyword (starts simulation)
        lines.append(f"PROCESS")
        lines.append(f"STOP")

        with open(key_path, "w") as key_file:
            key_file.write("\n".join(lines) + "\n")

        return key_path

    # =========================================================================
    # Simulation
    # =========================================================================

    def grow(self, years: int = 50) -> None:
        """Run the FVS simulation for the specified number of years.

        The simulation runs all cycles at once through the native FVS
        library, matching the keyword file configuration.

        Args:
            years: Number of years to simulate.

        Raises:
            FVSNativeError: If the simulation has not been initialized
                or if the native library returns an error.
        """
        if not self._initialized:
            raise FVSNativeError(
                "Stand not initialized. Call initialize_planted() first."
            )

        # Recalculate number of cycles and regenerate keyword file
        num_cycles = max(1, years // self._cycle_length)
        self._keyword_file = self._generate_keyword_file(
            species_index=get_species_index(self._species, self.variant),
            trees_per_acre=self._trees_per_acre,
            site_index=self._site_index,
            age=1,
            dbh=0.1,
            height=1.0,
            num_cycles=num_cycles,
        )

        # Initialize FVS with keyword file
        self._bindings.set_cmd_line(str(self._keyword_file))

        # Run simulation - FVS returns 2 when simulation is complete
        return_code = self._bindings.run()
        while return_code == 0:
            return_code = self._bindings.run()

        if return_code not in (0, 2):
            raise FVSNativeError(
                f"FVS simulation failed with return code {return_code}"
            )

        self._years_grown = years
        logger.info(
            "NativeStand simulation complete: %d years, %d cycles",
            years,
            num_cycles,
        )

    # =========================================================================
    # Output retrieval
    # =========================================================================

    def get_metrics(self) -> Dict[str, float]:
        """Get current stand metrics in the same format as Stand.get_metrics().

        Uses fvsEvmonAttr() for stand-level metrics — no reimplementation
        of SDI, BA, or other formulas.

        Returns:
            Dictionary with keys: tpa, basal_area, qmd, top_height,
            sdi, volume, age.
        """
        try:
            tpa = self._bindings.get_evmon_attr("atpa")
            basal_area = self._bindings.get_evmon_attr("aba")
            top_height = self._bindings.get_evmon_attr("atopht")
            sdi = self._bindings.get_evmon_attr("asdi")
            qmd = self._bindings.get_evmon_attr("aqmd")
        except FVSNativeError:
            # Fall back to before-growth values if after-growth not available
            tpa = self._bindings.get_evmon_attr("btpa")
            basal_area = self._bindings.get_evmon_attr("bba")
            top_height = self._bindings.get_evmon_attr("btopht")
            sdi = self._bindings.get_evmon_attr("bsdi")
            qmd = self._bindings.get_evmon_attr("bqmd")

        # Get volume from tree attributes
        try:
            tcuft_array = self._bindings.get_tree_attr("tcuft")
            tpa_array = self._bindings.get_tree_attr("tpa")
            volume = float(sum(tcuft_array * tpa_array)) if len(tcuft_array) > 0 else 0.0
        except FVSNativeError:
            volume = 0.0

        return {
            "tpa": tpa,
            "basal_area": basal_area,
            "qmd": qmd,
            "top_height": top_height,
            "sdi": sdi,
            "volume": volume,
            "age": self._years_grown,
        }

    def get_tree_list(self) -> List[Dict[str, float]]:
        """Get individual tree data in the same format as Stand tree lists.

        Returns:
            List of dictionaries, one per tree, with keys: species, dbh,
            height, tpa, crown_ratio, tcuft.
        """
        dims = self._bindings.get_dim_sizes()
        ntrees = dims["ntrees"]

        if ntrees <= 0:
            return []

        dbh = self._bindings.get_tree_attr("dbh", ntrees)
        ht = self._bindings.get_tree_attr("ht", ntrees)
        tpa = self._bindings.get_tree_attr("tpa", ntrees)
        cratio = self._bindings.get_tree_attr("cratio", ntrees)

        try:
            tcuft = self._bindings.get_tree_attr("tcuft", ntrees)
        except FVSNativeError:
            tcuft = [0.0] * ntrees

        try:
            species = self._bindings.get_tree_attr("species", ntrees)
        except FVSNativeError:
            species = [0] * ntrees

        tree_list = []
        for i in range(ntrees):
            tree_list.append({
                "species_index": int(species[i]) if i < len(species) else 0,
                "dbh": float(dbh[i]),
                "height": float(ht[i]),
                "tpa": float(tpa[i]),
                "crown_ratio": float(cratio[i]),
                "tcuft": float(tcuft[i]) if i < len(tcuft) else 0.0,
            })

        return tree_list

    def get_yield_table(self) -> List[Dict[str, float]]:
        """Get the complete yield table from the FVS simulation.

        Returns cycle-by-cycle summary data in a format compatible with
        the pyfvs validation system.

        Returns:
            List of dictionaries, one per cycle, with keys matching
            FVS summary output fields.
        """
        dims = self._bindings.get_dim_sizes()
        ncycles = dims["ncycles"]

        yield_table = []
        for cycle in range(ncycles + 1):  # Include initial conditions (cycle 0)
            try:
                summary = self._bindings.get_summary(cycle)
                yield_table.append(summary)
            except FVSNativeError:
                logger.warning("Could not retrieve summary for cycle %d", cycle)

        return yield_table
