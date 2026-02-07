"""
Platform-aware FVS shared library loader.

Discovers and loads native FVS shared libraries built from the USDA Forest
Service Fortran source code (https://github.com/USDAForestService/ForestVegetationSimulator).

Search order for library files:
    1. FVS_LIB_PATH environment variable (explicit path to directory)
    2. ~/.fvs/lib/
    3. /usr/local/lib/
    4. ./lib/ (relative to working directory)

Library naming convention: FVS{variant}.{ext}
    - macOS: FVSsn.dylib
    - Linux: FVSsn.so
    - Windows: FVSsn.dll

Usage:
    >>> from pyfvs.native.library_loader import load_fvs_library, fvs_library_available
    >>> if fvs_library_available('SN'):
    ...     lib = load_fvs_library('SN')
"""

import ctypes
import logging
import os
import platform
from pathlib import Path
from typing import Optional

from pyfvs.exceptions import FVSNativeError

logger = logging.getLogger(__name__)

# Singleton cache: variant -> loaded ctypes.CDLL
_library_cache: dict = {}


def _get_platform_extension() -> str:
    """Return the shared library extension for the current platform."""
    system = platform.system()
    if system == "Darwin":
        return ".dylib"
    elif system == "Windows":
        return ".dll"
    else:
        return ".so"


def _get_search_paths() -> list:
    """Return ordered list of directories to search for FVS libraries."""
    paths = []

    # 1. FVS_LIB_PATH environment variable
    env_path = os.environ.get("FVS_LIB_PATH")
    if env_path:
        paths.append(Path(env_path))

    # 2. ~/.fvs/lib/
    home_fvs = Path.home() / ".fvs" / "lib"
    paths.append(home_fvs)

    # 3. /usr/local/lib/
    paths.append(Path("/usr/local/lib"))

    # 4. ./lib/ relative to working directory
    paths.append(Path.cwd() / "lib")

    return paths


def _find_library_path(variant: str) -> Optional[Path]:
    """Search for the FVS shared library file for the given variant.

    Args:
        variant: FVS variant code (e.g., 'SN', 'PN', 'LS').

    Returns:
        Path to the library file, or None if not found.
    """
    ext = _get_platform_extension()
    # Try both cases: FVSsn and FVSSn
    lib_names = [
        f"FVS{variant.lower()}{ext}",
        f"FVS{variant.upper()}{ext}",
        f"fvs{variant.lower()}{ext}",
    ]

    for search_dir in _get_search_paths():
        if not search_dir.is_dir():
            continue
        for lib_name in lib_names:
            candidate = search_dir / lib_name
            if candidate.is_file():
                logger.debug("Found FVS library: %s", candidate)
                return candidate

    return None


def _detect_symbol_convention(lib: ctypes.CDLL) -> str:
    """Detect whether the library uses GCC or Intel Fortran symbol naming.

    GCC convention: lowercase with trailing underscore (e.g., fvs_)
    Intel convention: uppercase without underscore (e.g., FVS)

    Args:
        lib: Loaded ctypes library.

    Returns:
        'gcc' or 'intel' indicating the symbol convention.

    Raises:
        FVSNativeError: If neither convention is detected.
    """
    # Try GCC convention first (most common on macOS/Linux with gfortran)
    try:
        _ = lib.fvssetcmdline_
        return "gcc"
    except AttributeError:
        pass

    # Try Intel convention
    try:
        _ = lib.FVSSETCMDLINE
        return "intel"
    except AttributeError:
        pass

    raise FVSNativeError(
        "Cannot detect Fortran symbol convention in FVS library. "
        "Expected 'fvssetcmdline_' (GCC) or 'FVSSETCMDLINE' (Intel)."
    )


def get_symbol_name(base_name: str, convention: str) -> str:
    """Convert a base function name to the platform-specific symbol name.

    Args:
        base_name: Base function name (e.g., 'fvssetcmdline').
        convention: Either 'gcc' or 'intel'.

    Returns:
        Symbol name as it appears in the shared library.
    """
    if convention == "gcc":
        return f"{base_name.lower()}_"
    else:
        return base_name.upper()


def load_fvs_library(variant: str = "SN") -> ctypes.CDLL:
    """Load the native FVS shared library for the specified variant.

    The library is cached per variant, so subsequent calls return the
    same ctypes.CDLL instance.

    Args:
        variant: FVS variant code (e.g., 'SN', 'PN', 'LS').

    Returns:
        Loaded ctypes.CDLL instance.

    Raises:
        FVSNativeError: If the library cannot be found or loaded.
    """
    variant_upper = variant.upper()

    if variant_upper in _library_cache:
        return _library_cache[variant_upper]

    lib_path = _find_library_path(variant_upper)
    if lib_path is None:
        search_dirs = [str(p) for p in _get_search_paths()]
        raise FVSNativeError(
            f"FVS {variant_upper} shared library not found. "
            f"Searched: {search_dirs}. "
            f"Set FVS_LIB_PATH environment variable or install to ~/.fvs/lib/. "
            f"See pyfvs/native/BUILD.md for build instructions."
        )

    try:
        lib = ctypes.CDLL(str(lib_path))
    except OSError as exc:
        raise FVSNativeError(
            f"Failed to load FVS library at {lib_path}: {exc}"
        ) from exc

    # Detect and store the symbol convention
    convention = _detect_symbol_convention(lib)
    lib._fvs_symbol_convention = convention
    lib._fvs_variant = variant_upper
    logger.info(
        "Loaded FVS %s library from %s (convention: %s)",
        variant_upper,
        lib_path,
        convention,
    )

    _library_cache[variant_upper] = lib
    return lib


def fvs_library_available(variant: str = "SN") -> bool:
    """Check if the native FVS library is available for a given variant.

    This does NOT load the library — it only checks if the file exists.

    Args:
        variant: FVS variant code (e.g., 'SN', 'PN', 'LS').

    Returns:
        True if the library file exists in any search path.
    """
    return _find_library_path(variant.upper()) is not None


def get_library_info(variant: str = "SN") -> dict:
    """Return information about the FVS library for a given variant.

    Args:
        variant: FVS variant code.

    Returns:
        Dictionary with library location, platform, and availability info.
    """
    lib_path = _find_library_path(variant.upper())
    return {
        "variant": variant.upper(),
        "available": lib_path is not None,
        "path": str(lib_path) if lib_path else None,
        "platform": platform.system(),
        "extension": _get_platform_extension(),
        "search_paths": [str(p) for p in _get_search_paths()],
    }


def clear_library_cache():
    """Clear the library cache. Primarily useful for testing."""
    _library_cache.clear()
