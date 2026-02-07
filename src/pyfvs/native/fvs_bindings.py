"""
Low-level ctypes bindings to FVS API subroutines.

Wraps the Fortran API functions defined in apisubs.f from the USDA Forest
Service FVS source code. All functions use GCC Fortran calling conventions
by default (lowercase with trailing underscore, hidden string length args).

FVS API Reference (from apisubs.f):
    fvsSetCmdLine  - Initialize FVS with keyword file path
    fvs            - Run the full simulation
    fvsDimSizes    - Get dimension sizes (max trees, cycles, species)
    fvsTreeAttr    - Get/set individual tree attributes
    fvsSpeciesAttr - Get species-level attributes
    fvsEvmonAttr   - Get event monitor / stand-level attributes
    fvsSummary     - Get cycle summary output data
    fvsAddTrees    - Add trees programmatically
    fvsCutTrees    - Remove trees from the simulation

Important notes:
    - FVS uses Fortran COMMON blocks for global state. Only ONE simulation
      can be active per variant at a time within a process.
    - GCC Fortran passes string lengths as hidden trailing arguments.
    - Tree attribute arrays are Fortran-ordered (1-based indexing).
    - Return code 0 = success, non-zero = error.
"""

import ctypes
import logging
from typing import Optional

import numpy as np

from pyfvs.exceptions import FVSNativeError
from pyfvs.native.library_loader import load_fvs_library, get_symbol_name

logger = logging.getLogger(__name__)


class FVSBindings:
    """Typed ctypes bindings to native FVS API functions.

    Wraps the raw Fortran subroutines with Python-friendly signatures,
    handling string encoding, array allocation, and error checking.

    Args:
        variant: FVS variant code (e.g., 'SN', 'PN', 'LS').
    """

    def __init__(self, variant: str = "SN"):
        self.variant = variant.upper()
        self._lib = load_fvs_library(self.variant)
        self._convention = self._lib._fvs_symbol_convention

    def _get_func(self, base_name: str):
        """Get a function from the library using the detected naming convention."""
        symbol = get_symbol_name(base_name, self._convention)
        try:
            return getattr(self._lib, symbol)
        except AttributeError:
            raise FVSNativeError(
                f"Symbol '{symbol}' not found in FVS {self.variant} library. "
                f"The library may be an incompatible version."
            )

    def _check_return_code(self, return_code: int, function_name: str):
        """Raise FVSNativeError if return code indicates failure."""
        if return_code != 0:
            raise FVSNativeError(
                f"FVS {function_name}() returned error code {return_code}"
            )

    # =========================================================================
    # Initialization
    # =========================================================================

    def set_cmd_line(self, keyword_file: str) -> int:
        """Initialize FVS with a keyword file path.

        Corresponds to: SUBROUTINE FVSSETCMDLINE(CMD, CMDLEN, IRTNCD)

        Args:
            keyword_file: Path to the FVS keyword (.key) file.

        Returns:
            Return code (0 = success).
        """
        func = self._get_func("fvssetcmdline")

        cmd_bytes = keyword_file.encode("ascii")
        cmd_len = len(cmd_bytes)
        cmd_buf = ctypes.create_string_buffer(cmd_bytes, cmd_len)
        cmd_len_c = ctypes.c_int(cmd_len)
        return_code = ctypes.c_int(0)

        # GCC Fortran: hidden string length at end of argument list
        func(cmd_buf, ctypes.byref(cmd_len_c), ctypes.byref(return_code),
             ctypes.c_int(cmd_len))

        self._check_return_code(return_code.value, "fvsSetCmdLine")
        return return_code.value

    def run(self) -> int:
        """Run the FVS simulation to completion.

        Corresponds to: SUBROUTINE FVS(IRTNCD)

        Returns:
            Return code (0 = success, 2 = simulation complete).
        """
        func = self._get_func("fvs")
        return_code = ctypes.c_int(0)
        func(ctypes.byref(return_code))
        return return_code.value

    # =========================================================================
    # Dimension queries
    # =========================================================================

    def get_dim_sizes(self) -> dict:
        """Get FVS dimension sizes (max trees, cycles, species, etc.).

        Corresponds to: SUBROUTINE FVSDIMSIZES(NTREES, NCYCLES, NPLOTS,
                                                MAXTREES, MAXSPECIES, MAXPLOTS,
                                                MAXCYCLES, IRTNCD)

        Returns:
            Dictionary with keys: ntrees, ncycles, nplots, maxtrees,
            maxspecies, maxplots, maxcycles.
        """
        func = self._get_func("fvsdimsizes")

        ntrees = ctypes.c_int(0)
        ncycles = ctypes.c_int(0)
        nplots = ctypes.c_int(0)
        maxtrees = ctypes.c_int(0)
        maxspecies = ctypes.c_int(0)
        maxplots = ctypes.c_int(0)
        maxcycles = ctypes.c_int(0)
        return_code = ctypes.c_int(0)

        func(
            ctypes.byref(ntrees),
            ctypes.byref(ncycles),
            ctypes.byref(nplots),
            ctypes.byref(maxtrees),
            ctypes.byref(maxspecies),
            ctypes.byref(maxplots),
            ctypes.byref(maxcycles),
            ctypes.byref(return_code),
        )

        self._check_return_code(return_code.value, "fvsDimSizes")

        return {
            "ntrees": ntrees.value,
            "ncycles": ncycles.value,
            "nplots": nplots.value,
            "maxtrees": maxtrees.value,
            "maxspecies": maxspecies.value,
            "maxplots": maxplots.value,
            "maxcycles": maxcycles.value,
        }

    # =========================================================================
    # Tree attributes
    # =========================================================================

    def get_tree_attr(self, attr_name: str, ntrees: Optional[int] = None) -> np.ndarray:
        """Get a tree-level attribute array from FVS.

        Corresponds to: SUBROUTINE FVSTREEATTR(ESSION, ATTRNAME, IATTRL,
                                                NVALS, VALS, IRTNCD)

        Common attribute names:
            - 'dbh'   : Diameter at breast height (inches)
            - 'ht'    : Total height (feet)
            - 'tpa'   : Trees per acre
            - 'cratio': Crown ratio (0-1)
            - 'tcuft' : Total cubic foot volume
            - 'mcuft' : Merchantable cubic foot volume
            - 'bdft'  : Board foot volume
            - 'species': Species index (integer)
            - 'age'   : Tree age (years)

        Args:
            attr_name: Attribute name string.
            ntrees: Number of trees (if None, queries dim sizes first).

        Returns:
            numpy array of attribute values (length = ntrees).
        """
        if ntrees is None:
            dims = self.get_dim_sizes()
            ntrees = dims["ntrees"]

        if ntrees <= 0:
            return np.array([], dtype=np.float64)

        func = self._get_func("fvstreeattr")

        attr_bytes = attr_name.encode("ascii")
        attr_len = len(attr_bytes)
        attr_buf = ctypes.create_string_buffer(attr_bytes, attr_len)
        attr_len_c = ctypes.c_int(attr_len)
        nvals = ctypes.c_int(ntrees)
        vals = (ctypes.c_double * ntrees)()
        return_code = ctypes.c_int(0)

        func(
            attr_buf,
            ctypes.byref(attr_len_c),
            ctypes.byref(nvals),
            vals,
            ctypes.byref(return_code),
            ctypes.c_int(attr_len),  # hidden string length
        )

        self._check_return_code(return_code.value, f"fvsTreeAttr('{attr_name}')")
        return np.array(vals[:nvals.value], dtype=np.float64)

    def set_tree_attr(self, attr_name: str, values: np.ndarray) -> None:
        """Set a tree-level attribute array in FVS.

        Corresponds to: SUBROUTINE FVSSETATTRDB(ATTRNAME, IATTRL,
                                                 NVALS, VALS, IRTNCD)

        Args:
            attr_name: Attribute name string (e.g., 'dbh', 'ht', 'tpa').
            values: numpy array of values to set.
        """
        func = self._get_func("fvssetattrdb")

        attr_bytes = attr_name.encode("ascii")
        attr_len = len(attr_bytes)
        attr_buf = ctypes.create_string_buffer(attr_bytes, attr_len)
        attr_len_c = ctypes.c_int(attr_len)
        nvals = ctypes.c_int(len(values))
        vals = (ctypes.c_double * len(values))(*values)
        return_code = ctypes.c_int(0)

        func(
            attr_buf,
            ctypes.byref(attr_len_c),
            ctypes.byref(nvals),
            vals,
            ctypes.byref(return_code),
            ctypes.c_int(attr_len),
        )

        self._check_return_code(return_code.value, f"fvsSetAttrDb('{attr_name}')")

    # =========================================================================
    # Species attributes
    # =========================================================================

    def get_species_attr(self, attr_name: str) -> np.ndarray:
        """Get a species-level attribute array from FVS.

        Corresponds to: SUBROUTINE FVSSPECIESATTR(ATTRNAME, IATTRL,
                                                   NVALS, VALS, IRTNCD)

        Common attribute names:
            - 'spmcdbh' : Species minimum merchantable DBH
            - 'sptpa'   : Species trees per acre

        Args:
            attr_name: Attribute name string.

        Returns:
            numpy array of species attribute values.
        """
        dims = self.get_dim_sizes()
        nspecies = dims["maxspecies"]

        func = self._get_func("fvsspeciesattr")

        attr_bytes = attr_name.encode("ascii")
        attr_len = len(attr_bytes)
        attr_buf = ctypes.create_string_buffer(attr_bytes, attr_len)
        attr_len_c = ctypes.c_int(attr_len)
        nvals = ctypes.c_int(nspecies)
        vals = (ctypes.c_double * nspecies)()
        return_code = ctypes.c_int(0)

        func(
            attr_buf,
            ctypes.byref(attr_len_c),
            ctypes.byref(nvals),
            vals,
            ctypes.byref(return_code),
            ctypes.c_int(attr_len),
        )

        self._check_return_code(return_code.value, f"fvsSpeciesAttr('{attr_name}')")
        return np.array(vals[:nvals.value], dtype=np.float64)

    # =========================================================================
    # Stand-level / Event Monitor attributes
    # =========================================================================

    def get_evmon_attr(self, attr_name: str) -> float:
        """Get a stand-level attribute from the FVS Event Monitor.

        Corresponds to: SUBROUTINE FVSEVMONATTR(ATTRNAME, IATTRL,
                                                  NVALS, VALS, IRTNCD)

        Common attribute names:
            - 'bba'    : Before-growth basal area (sq ft/acre)
            - 'btpa'   : Before-growth trees per acre
            - 'bsdi'   : Before-growth Stand Density Index
            - 'btopht' : Before-growth top height (feet)
            - 'bqmd'   : Before-growth QMD (inches)
            - 'aba'    : After-growth basal area
            - 'atpa'   : After-growth trees per acre
            - 'asdi'   : After-growth SDI
            - 'atopht' : After-growth top height
            - 'aqmd'   : After-growth QMD

        Args:
            attr_name: Attribute name string.

        Returns:
            Scalar float value.
        """
        func = self._get_func("fvsevmonattr")

        attr_bytes = attr_name.encode("ascii")
        attr_len = len(attr_bytes)
        attr_buf = ctypes.create_string_buffer(attr_bytes, attr_len)
        attr_len_c = ctypes.c_int(attr_len)
        nvals = ctypes.c_int(1)
        val = ctypes.c_double(0.0)
        return_code = ctypes.c_int(0)

        func(
            attr_buf,
            ctypes.byref(attr_len_c),
            ctypes.byref(nvals),
            ctypes.byref(val),
            ctypes.byref(return_code),
            ctypes.c_int(attr_len),
        )

        self._check_return_code(return_code.value, f"fvsEvmonAttr('{attr_name}')")
        return val.value

    # =========================================================================
    # Summary output
    # =========================================================================

    def get_summary(self, cycle: int) -> dict:
        """Get summary output data for a given cycle.

        Corresponds to: SUBROUTINE FVSSUMMARY(ICYC, ISUMARY, OSUMARY, IRTNCD)

        The summary array contains standard FVS output table values
        for the specified cycle.

        Args:
            cycle: Cycle number (0 = initial conditions, 1+ = growth cycles).

        Returns:
            Dictionary with summary metrics for the cycle.
        """
        func = self._get_func("fvssummary")

        icyc = ctypes.c_int(cycle)
        # FVS summary has ~20 output values
        nsumary = 20
        isumary = ctypes.c_int(nsumary)
        osumary = (ctypes.c_double * nsumary)()
        return_code = ctypes.c_int(0)

        func(
            ctypes.byref(icyc),
            ctypes.byref(isumary),
            osumary,
            ctypes.byref(return_code),
        )

        self._check_return_code(return_code.value, f"fvsSummary(cycle={cycle})")

        # Map summary array positions to named fields
        # Positions based on FVS summary table output (sumary.f)
        vals = list(osumary)
        return {
            "cycle": cycle,
            "year": int(vals[0]) if len(vals) > 0 else 0,
            "age": int(vals[1]) if len(vals) > 1 else 0,
            "tpa": vals[2] if len(vals) > 2 else 0.0,
            "ba": vals[3] if len(vals) > 3 else 0.0,
            "sdi": vals[4] if len(vals) > 4 else 0.0,
            "ccf": vals[5] if len(vals) > 5 else 0.0,
            "top_height": vals[6] if len(vals) > 6 else 0.0,
            "qmd": vals[7] if len(vals) > 7 else 0.0,
            "total_cuft": vals[8] if len(vals) > 8 else 0.0,
            "merch_cuft": vals[9] if len(vals) > 9 else 0.0,
            "merch_bdft": vals[10] if len(vals) > 10 else 0.0,
            "removed_tpa": vals[11] if len(vals) > 11 else 0.0,
            "removed_cuft": vals[12] if len(vals) > 12 else 0.0,
            "removed_bdft": vals[13] if len(vals) > 13 else 0.0,
            "mortality_tpa": vals[14] if len(vals) > 14 else 0.0,
            "mortality_cuft": vals[15] if len(vals) > 15 else 0.0,
            "mai_cuft": vals[16] if len(vals) > 16 else 0.0,
            "forest_type": int(vals[17]) if len(vals) > 17 else 0,
            "size_class": int(vals[18]) if len(vals) > 18 else 0,
            "stocking_class": int(vals[19]) if len(vals) > 19 else 0,
        }

    # =========================================================================
    # Tree manipulation
    # =========================================================================

    def add_trees(
        self,
        species_indices: np.ndarray,
        dbh_values: np.ndarray,
        height_values: np.ndarray,
        tpa_values: np.ndarray,
        crown_ratios: Optional[np.ndarray] = None,
    ) -> None:
        """Add trees to the FVS simulation programmatically.

        Corresponds to: SUBROUTINE FVSADDTREES(NTREES, ISPN, DBH, HT,
                                                CRATIO, PLOT, IRTNCD)

        Args:
            species_indices: Array of integer species indices.
            dbh_values: Array of DBH values (inches).
            height_values: Array of total heights (feet).
            tpa_values: Array of trees per acre values.
            crown_ratios: Optional array of crown ratios (0-100 scale).
        """
        func = self._get_func("fvsaddtrees")

        ntrees = len(species_indices)
        if not (len(dbh_values) == len(height_values) == len(tpa_values) == ntrees):
            raise ValueError("All input arrays must have the same length")

        ntrees_c = ctypes.c_int(ntrees)
        ispn = (ctypes.c_int * ntrees)(*species_indices.astype(int))
        dbh = (ctypes.c_double * ntrees)(*dbh_values)
        ht = (ctypes.c_double * ntrees)(*height_values)

        if crown_ratios is not None:
            cr = (ctypes.c_double * ntrees)(*crown_ratios)
        else:
            cr = (ctypes.c_double * ntrees)(*([0.0] * ntrees))

        # Plot number (1 for single-plot stands)
        plot = (ctypes.c_int * ntrees)(*([1] * ntrees))
        return_code = ctypes.c_int(0)

        func(
            ctypes.byref(ntrees_c),
            ispn,
            dbh,
            ht,
            cr,
            plot,
            ctypes.byref(return_code),
        )

        self._check_return_code(return_code.value, "fvsAddTrees")

    def cut_trees(self, tree_indices: np.ndarray) -> None:
        """Remove trees from the FVS simulation.

        Corresponds to: SUBROUTINE FVSCUTTREES(NTREES, ITREELIST, IRTNCD)

        Args:
            tree_indices: Array of 1-based tree indices to remove.
        """
        func = self._get_func("fvscuttrees")

        ntrees = len(tree_indices)
        ntrees_c = ctypes.c_int(ntrees)
        indices = (ctypes.c_int * ntrees)(*tree_indices.astype(int))
        return_code = ctypes.c_int(0)

        func(
            ctypes.byref(ntrees_c),
            indices,
            ctypes.byref(return_code),
        )

        self._check_return_code(return_code.value, "fvsCutTrees")
