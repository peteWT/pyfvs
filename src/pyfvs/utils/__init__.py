"""
Utility functions for FVS-Python.

This module provides common utilities used throughout the codebase.
"""

from .string_utils import normalize_code, normalize_species_code, normalize_ecounit
from .age_lookup import (
    generate_dbh_age_table,
    get_age_at_dbh,
    get_dbh_at_age,
    estimate_age_from_dbh,
    get_years_to_merchantable,
)

__all__ = [
    "normalize_code",
    "normalize_species_code",
    "normalize_ecounit",
    "generate_dbh_age_table",
    "get_age_at_dbh",
    "get_dbh_at_age",
    "estimate_age_from_dbh",
    "get_years_to_merchantable",
]
