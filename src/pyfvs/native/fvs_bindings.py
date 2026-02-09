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
                                                MAXTREES, MAXSPECIES,
                                                MAXPLOTS, MAXCYCLES)

        Note: fvsDimSizes has no return code parameter in the Fortran source.

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

        func(
            ctypes.byref(ntrees),
            ctypes.byref(ncycles),
            ctypes.byref(nplots),
            ctypes.byref(maxtrees),
            ctypes.byref(maxspecies),
            ctypes.byref(maxplots),
            ctypes.byref(maxcycles),
        )

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

        Corresponds to: SUBROUTINE FVSTREEATTR(NAME, NCH, ACTION, NTREES, ATTR, RTNCODE)

        Common attribute names (case sensitive):
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
            attr_name: Attribute name string (case sensitive).
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

        # Name parameter
        name_bytes = attr_name.encode("ascii")
        name_len = len(name_bytes)
        name_buf = ctypes.create_string_buffer(name_bytes, name_len)
        nch = ctypes.c_int(name_len)

        # Action parameter - "get"
        action_bytes = b"get"
        action_len = len(action_bytes)
        action_buf = ctypes.create_string_buffer(action_bytes, action_len)

        nvals = ctypes.c_int(ntrees)
        vals = (ctypes.c_double * ntrees)()
        return_code = ctypes.c_int(0)

        # GCC Fortran: hidden string lengths at end of arg list
        func(
            name_buf,
            ctypes.byref(nch),
            action_buf,
            ctypes.byref(nvals),
            vals,
            ctypes.byref(return_code),
            ctypes.c_int(name_len),    # hidden length for name
            ctypes.c_int(action_len),  # hidden length for action
        )

        self._check_return_code(return_code.value, f"fvsTreeAttr('{attr_name}')")
        return np.array(vals[:nvals.value], dtype=np.float64)

    def set_tree_attr(self, attr_name: str, values: np.ndarray) -> None:
        """Set a tree-level attribute array in FVS.

        Uses fvsTreeAttr with action="set" to modify tree attributes.

        Args:
            attr_name: Attribute name string (e.g., 'dbh', 'ht', 'tpa').
            values: numpy array of values to set.
        """
        func = self._get_func("fvstreeattr")

        # Name parameter
        name_bytes = attr_name.encode("ascii")
        name_len = len(name_bytes)
        name_buf = ctypes.create_string_buffer(name_bytes, name_len)
        nch = ctypes.c_int(name_len)

        # Action parameter - "set"
        action_bytes = b"set"
        action_len = len(action_bytes)
        action_buf = ctypes.create_string_buffer(action_bytes, action_len)

        nvals = ctypes.c_int(len(values))
        vals = (ctypes.c_double * len(values))(*values)
        return_code = ctypes.c_int(0)

        # GCC Fortran: hidden string lengths at end of arg list
        func(
            name_buf,
            ctypes.byref(nch),
            action_buf,
            ctypes.byref(nvals),
            vals,
            ctypes.byref(return_code),
            ctypes.c_int(name_len),    # hidden length for name
            ctypes.c_int(action_len),  # hidden length for action
        )

        self._check_return_code(return_code.value, f"fvsTreeAttr('{attr_name}', set)")

    # =========================================================================
    # Species attributes
    # =========================================================================

    def get_species_attr(self, attr_name: str) -> np.ndarray:
        """Get a species-level attribute array from FVS.

        Corresponds to: SUBROUTINE FVSSPECIESATTR(NAME, NCH, ACTION, ATTR, RTNCODE)

        Common attribute names (case sensitive):
            - 'spmcdbh' : Species minimum merchantable DBH
            - 'sptpa'   : Species trees per acre

        Args:
            attr_name: Attribute name string (case sensitive).

        Returns:
            numpy array of species attribute values.
        """
        dims = self.get_dim_sizes()
        nspecies = dims["maxspecies"]

        func = self._get_func("fvsspeciesattr")

        # Name parameter
        name_bytes = attr_name.encode("ascii")
        name_len = len(name_bytes)
        name_buf = ctypes.create_string_buffer(name_bytes, name_len)
        nch = ctypes.c_int(name_len)

        # Action parameter - "get"
        action_bytes = b"get"
        action_len = len(action_bytes)
        action_buf = ctypes.create_string_buffer(action_bytes, action_len)

        vals = (ctypes.c_double * nspecies)()
        return_code = ctypes.c_int(0)

        # GCC Fortran: hidden string lengths at end of arg list
        func(
            name_buf,
            ctypes.byref(nch),
            action_buf,
            vals,
            ctypes.byref(return_code),
            ctypes.c_int(name_len),    # hidden length for name
            ctypes.c_int(action_len),  # hidden length for action
        )

        self._check_return_code(return_code.value, f"fvsSpeciesAttr('{attr_name}')")
        return np.array(vals[:nspecies], dtype=np.float64)

    # =========================================================================
    # Stand-level / Event Monitor attributes
    # =========================================================================

    def get_evmon_attr(self, attr_name: str) -> float:
        """Get a stand-level attribute from the FVS Event Monitor.

        Corresponds to: SUBROUTINE FVSEVMONATTR(NAME, NCH, ACTION, ATTR, RTNCODE)

        Common attribute names (case sensitive):
            - 'bba'    : Before-growth basal area (sq ft/acre)
            - 'btpa'   : Before-growth trees per acre
            - 'bsdi'   : Before-growth Stand Density Index
            - 'btopht' : Before-growth top height (feet)
            - 'badbh'  : Before-growth average DBH (inches)
            - 'aba'    : After-growth basal area
            - 'atpa'   : After-growth trees per acre
            - 'asdi'   : After-growth SDI
            - 'atopht' : After-growth top height
            - 'aadbh'  : After-growth average DBH
            - 'age'    : Stand age
            - 'year'   : Calendar year
            - 'cycle'  : Current cycle number
            - 'site'   : Site index

        Args:
            attr_name: Attribute name string (case sensitive).

        Returns:
            Scalar float value.
        """
        func = self._get_func("fvsevmonattr")

        # Name parameter
        name_bytes = attr_name.encode("ascii")
        name_len = len(name_bytes)
        name_buf = ctypes.create_string_buffer(name_bytes, name_len)
        nch = ctypes.c_int(name_len)

        # Action parameter - "get"
        action_bytes = b"get"
        action_len = len(action_bytes)
        action_buf = ctypes.create_string_buffer(action_bytes, action_len)

        val = ctypes.c_double(0.0)
        return_code = ctypes.c_int(0)

        # GCC Fortran hidden string lengths go at the end of the arg list
        func(
            name_buf,
            ctypes.byref(nch),
            action_buf,
            ctypes.byref(val),
            ctypes.byref(return_code),
            ctypes.c_int(name_len),    # hidden length for name
            ctypes.c_int(action_len),  # hidden length for action
        )

        self._check_return_code(return_code.value, f"fvsEvmonAttr('{attr_name}')")
        return val.value

    # =========================================================================
    # Summary output
    # =========================================================================

    def get_summary(self, cycle: int) -> dict:
        """Get summary output data for a given cycle.

        Corresponds to: SUBROUTINE FVSSUMMARY(SUMMARY, ICYCLE, NCYCLES,
                                               MAXROW, MAXCOL, RTNCODE)

        The summary array is INTEGER*4(22) from IOSUM in OUTCOM.F77.
        Values are integer-encoded (e.g., BA in sq ft, TPA as integer).

        IOSUM columns (from FVS source vbase/sumout.f):
             1: Year
             2: Age
             3: Trees/acre (start of period)
             4: Total cu ft (start of period)
             5: Merch cu ft (start of period)
             6: Merch bd ft (start of period)
             7: Removed trees/acre
             8: Removed total cu ft
             9: Removed merch cu ft
            10: Removed merch bd ft
            11: Basal area/acre (after treatment)
            12: CCF (after treatment)
            13: Avg dominant height (after treatment)
            14: Period length (years)
            15: Accretion (annual cu ft/acre)
            16: Mortality (annual cu ft/acre)
            17: Sample weight
            18: Forest cover type code
            19: Size class
            20: Stocking class
            21: Cubic saw volume (start of period)
            22: Removed cubic saw volume

        Note: Begin-of-period BA, SDI, CCF, TopHt, QMD and after-treatment
        SDI, QMD, MAI are stored in the SUMTAB common block (IOLDBA, ISDI,
        IBTCCF, IBTAVH, QSDBT, ISDIAT, QDBHAT, BCYMAI), NOT in IOSUM.
        They are not accessible through fvsSummary.

        Args:
            cycle: Cycle number (1 = initial conditions, 2+ = growth cycles).
                   Range: 1 to ncycles+1.

        Returns:
            Dictionary with summary metrics for the cycle.
        """
        func = self._get_func("fvssummary")

        # summary is INTEGER*4(22) output array
        summary = (ctypes.c_int * 22)()
        icycle = ctypes.c_int(cycle)
        ncycles = ctypes.c_int(0)
        maxrow = ctypes.c_int(0)
        maxcol = ctypes.c_int(0)
        return_code = ctypes.c_int(0)

        func(
            summary,
            ctypes.byref(icycle),
            ctypes.byref(ncycles),
            ctypes.byref(maxrow),
            ctypes.byref(maxcol),
            ctypes.byref(return_code),
        )

        self._check_return_code(return_code.value, f"fvsSummary(cycle={cycle})")

        vals = [summary[i] for i in range(22)]
        return {
            "cycle": cycle,
            "ncycles": ncycles.value,
            "year": vals[0],
            "age": vals[1],
            "begin_tpa": vals[2],
            "begin_tcuft": vals[3],
            "begin_mcuft": vals[4],
            "begin_bdft": vals[5],
            "removed_tpa": vals[6],
            "removed_tcuft": vals[7],
            "removed_mcuft": vals[8],
            "removed_bdft": vals[9],
            "after_ba": vals[10],
            "after_ccf": vals[11],
            "after_top_ht": vals[12],
            "period_length": vals[13],
            "accretion_cuft": vals[14],
            "mortality_cuft": vals[15],
            "sample_weight": vals[16],
            "forest_type": vals[17],
            "size_class": vals[18],
            "stocking_class": vals[19],
            "begin_scuft": vals[20],
            "removed_scuft": vals[21],
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
        plot_indices: Optional[np.ndarray] = None,
    ) -> None:
        """Add trees to the FVS simulation programmatically.

        Corresponds to: SUBROUTINE FVSADDTREES(IN_DBH, IN_SPECIES, IN_HT,
                                                IN_CRATIO, IN_PLOT, IN_TPA,
                                                NTREES, RTNCODE)

        All array parameters are real(kind=8) (double precision) in Fortran,
        including species indices and plot numbers which are converted to
        integer internally by FVS.

        Args:
            species_indices: Array of integer species indices.
            dbh_values: Array of DBH values (inches).
            height_values: Array of total heights (feet).
            tpa_values: Array of trees per acre values.
            crown_ratios: Optional array of crown ratios (0-100 scale).
            plot_indices: Optional array of plot numbers (default 1).
        """
        func = self._get_func("fvsaddtrees")

        ntrees = len(species_indices)
        if not (len(dbh_values) == len(height_values) == len(tpa_values) == ntrees):
            raise ValueError("All input arrays must have the same length")

        # All arrays are real(kind=8) in Fortran — use c_double for everything
        in_dbh = (ctypes.c_double * ntrees)(*dbh_values)
        in_species = (ctypes.c_double * ntrees)(*species_indices.astype(float))
        in_ht = (ctypes.c_double * ntrees)(*height_values)

        if crown_ratios is not None:
            in_cratio = (ctypes.c_double * ntrees)(*crown_ratios)
        else:
            in_cratio = (ctypes.c_double * ntrees)(*([0.0] * ntrees))

        if plot_indices is not None:
            in_plot = (ctypes.c_double * ntrees)(*plot_indices.astype(float))
        else:
            in_plot = (ctypes.c_double * ntrees)(*([1.0] * ntrees))

        in_tpa = (ctypes.c_double * ntrees)(*tpa_values)
        ntrees_c = ctypes.c_int(ntrees)
        return_code = ctypes.c_int(0)

        # Fortran order: in_dbh, in_species, in_ht, in_cratio,
        #                in_plot, in_tpa, ntrees, rtnCode
        func(
            in_dbh,
            in_species,
            in_ht,
            in_cratio,
            in_plot,
            in_tpa,
            ctypes.byref(ntrees_c),
            ctypes.byref(return_code),
        )

        self._check_return_code(return_code.value, "fvsAddTrees")

    def cut_trees(self, cut_flags: np.ndarray) -> None:
        """Remove trees from the FVS simulation.

        Corresponds to: SUBROUTINE FVSCUTTREES(PTOCUT, NTREES, RTNCODE)

        Note: This subroutine is marked "Not yet implemented" in the FVS
        Fortran source (apisubs.f) and always returns rtnCode=1.

        Args:
            cut_flags: Double precision array of length ntrees. Non-zero
                values mark trees for cutting.
        """
        func = self._get_func("fvscuttrees")

        ntrees = len(cut_flags)
        # pToCut is double precision in Fortran
        p_to_cut = (ctypes.c_double * ntrees)(*cut_flags.astype(float))
        ntrees_c = ctypes.c_int(ntrees)
        return_code = ctypes.c_int(0)

        # Fortran order: pToCut, ntrees, rtnCode
        func(
            p_to_cut,
            ctypes.byref(ntrees_c),
            ctypes.byref(return_code),
        )

        self._check_return_code(return_code.value, "fvsCutTrees")
