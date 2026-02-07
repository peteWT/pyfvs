"""
Native FVS bindings for validating pyfvs against the USDA Fortran implementation.

This subpackage provides ctypes bindings to the official FVS shared libraries.
All imports are lazy — importing pyfvs or pyfvs.native will never fail even
when the FVS shared library is not installed.

Quick check:
    >>> from pyfvs.native import fvs_library_available
    >>> fvs_library_available('SN')
    False  # Library not installed

When the library IS available:
    >>> from pyfvs.native import NativeStand
    >>> with NativeStand(variant='SN') as ns:
    ...     ns.initialize_planted(500, 70, 'LP')
    ...     ns.grow(50)
    ...     print(ns.get_metrics())

See native/BUILD.md for instructions on building the FVS shared libraries.
"""

# Species map is pure data — always importable
from pyfvs.native.species_map import (
    get_species_index,
    get_species_code,
    get_supported_variants,
    get_species_count,
    VARIANT_SPECIES_MAPS,
)

# Library availability check — always importable (doesn't load the library)
from pyfvs.native.library_loader import fvs_library_available, get_library_info


def __getattr__(name):
    """Lazy imports for classes that require the FVS shared library."""
    if name == "NativeStand":
        from pyfvs.native.native_stand import NativeStand
        return NativeStand
    if name == "FVSBindings":
        from pyfvs.native.fvs_bindings import FVSBindings
        return FVSBindings
    if name == "load_fvs_library":
        from pyfvs.native.library_loader import load_fvs_library
        return load_fvs_library
    raise AttributeError(f"module 'pyfvs.native' has no attribute {name!r}")


__all__ = [
    # Always available (no library needed)
    "fvs_library_available",
    "get_library_info",
    "get_species_index",
    "get_species_code",
    "get_supported_variants",
    "get_species_count",
    "VARIANT_SPECIES_MAPS",
    # Lazy-loaded (require FVS library)
    "NativeStand",
    "FVSBindings",
    "load_fvs_library",
]
