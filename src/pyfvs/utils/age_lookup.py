"""
DBH-to-age lookup utilities for FVS-Python.

This module provides functions to generate and query DBH-age lookup tables
from FVS growth simulations. These are used for the Bruck et al. discount
timber price method, which requires estimating tree age from diameter.
"""

from functools import lru_cache
from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np


def generate_dbh_age_table(
    species: str,
    site_index: float,
    trees_per_acre: int = 500,
    max_age: int = 50,
    time_step: int = 1,
) -> pd.DataFrame:
    """
    Generate a lookup table mapping age to DBH for a species.

    Uses the SimulationEngine to run a growth simulation and extract
    the DBH-age trajectory. The result is a table where each row
    represents stand metrics at a specific age.

    Args:
        species: Species code (e.g., 'LP', 'SP', 'SA', 'LL')
        site_index: Site index (base age 25) in feet
        trees_per_acre: Initial planting density (default 500)
        max_age: Maximum age to simulate (default 50 years)
        time_step: Years between growth periods (default 1 for finer resolution)

    Returns:
        DataFrame with columns: age, mean_dbh, qmd, tpa, volume
    """
    # Import here to avoid circular imports
    from pyfvs.stand import Stand

    # Initialize stand
    stand = Stand.initialize_planted(
        trees_per_acre=trees_per_acre,
        site_index=site_index,
        species=species
    )

    # Collect metrics at each age
    records = []

    # Initial metrics
    metrics = stand.get_metrics()
    records.append({
        'age': metrics['age'],
        'mean_dbh': metrics['mean_dbh'],
        'qmd': metrics.get('qmd', metrics['mean_dbh']),
        'tpa': metrics['tpa'],
        'volume': metrics['volume'],
    })

    # Grow and collect metrics
    for age in range(time_step, max_age + 1, time_step):
        stand.grow(years=time_step)
        metrics = stand.get_metrics()
        records.append({
            'age': metrics['age'],
            'mean_dbh': metrics['mean_dbh'],
            'qmd': metrics.get('qmd', metrics['mean_dbh']),
            'tpa': metrics['tpa'],
            'volume': metrics['volume'],
        })

    return pd.DataFrame(records)


def get_age_at_dbh(
    dbh_age_table: pd.DataFrame,
    target_dbh: float,
    dbh_column: str = 'mean_dbh',
) -> Optional[int]:
    """
    Interpolate age from DBH using a pre-computed growth curve.

    Uses linear interpolation between growth increments. If the target DBH
    is below the minimum in the table, returns 0. If above the maximum,
    returns the maximum age.

    Args:
        dbh_age_table: DataFrame with 'age' and DBH columns
        target_dbh: Target diameter at breast height (inches)
        dbh_column: Column name for DBH values (default 'mean_dbh')

    Returns:
        Estimated age in years, or None if interpolation fails
    """
    if dbh_age_table.empty:
        return None

    # Get age and DBH arrays
    ages = dbh_age_table['age'].values
    dbhs = dbh_age_table[dbh_column].values

    # Handle edge cases
    if target_dbh <= dbhs[0]:
        return int(ages[0])
    if target_dbh >= dbhs[-1]:
        return int(ages[-1])

    # Linear interpolation
    for i in range(len(dbhs) - 1):
        if dbhs[i] <= target_dbh <= dbhs[i + 1]:
            # Linear interpolation between points
            dbh_range = dbhs[i + 1] - dbhs[i]
            if dbh_range == 0:
                return int(ages[i])

            fraction = (target_dbh - dbhs[i]) / dbh_range
            age_range = ages[i + 1] - ages[i]
            interpolated_age = ages[i] + fraction * age_range
            return int(round(interpolated_age))

    return None


def get_dbh_at_age(
    dbh_age_table: pd.DataFrame,
    target_age: int,
    dbh_column: str = 'mean_dbh',
) -> Optional[float]:
    """
    Interpolate DBH from age using a pre-computed growth curve.

    Args:
        dbh_age_table: DataFrame with 'age' and DBH columns
        target_age: Target age in years
        dbh_column: Column name for DBH values

    Returns:
        Estimated DBH in inches, or None if interpolation fails
    """
    if dbh_age_table.empty:
        return None

    ages = dbh_age_table['age'].values
    dbhs = dbh_age_table[dbh_column].values

    # Handle edge cases
    if target_age <= ages[0]:
        return float(dbhs[0])
    if target_age >= ages[-1]:
        return float(dbhs[-1])

    # Linear interpolation
    return float(np.interp(target_age, ages, dbhs))


# Species code to name mapping for caching
SPECIES_CODES = {
    'LP': 'loblolly',
    'SP': 'shortleaf',
    'SA': 'slash',
    'LL': 'longleaf',
}


@lru_cache(maxsize=32)
def get_cached_dbh_age_table(
    species: str,
    site_index: float,
    trees_per_acre: int = 500,
    max_age: int = 50,
) -> Tuple[Tuple[int, ...], Tuple[float, ...]]:
    """
    Get cached DBH-age table as tuples for use in cached functions.

    Args:
        species: Species code
        site_index: Site index
        trees_per_acre: Initial TPA
        max_age: Maximum age

    Returns:
        Tuple of (ages, dbhs) as tuples for caching
    """
    df = generate_dbh_age_table(species, site_index, trees_per_acre, max_age)
    return (
        tuple(df['age'].values),
        tuple(df['mean_dbh'].values),
    )


def estimate_age_from_dbh(
    dbh: float,
    species: str,
    site_index: float = 65.0,
    trees_per_acre: int = 500,
) -> int:
    """
    Estimate tree age from DBH using FVS growth curves.

    This is the main entry point for age estimation. It generates a
    growth curve (cached) and interpolates the age for the given DBH.

    Args:
        dbh: Diameter at breast height (inches)
        species: Species code (LP, SP, SA, LL)
        site_index: Site index (base age 25), default 65 for southern region
        trees_per_acre: Planting density assumption (default 500)

    Returns:
        Estimated age in years
    """
    ages, dbhs = get_cached_dbh_age_table(species, site_index, trees_per_acre)

    # Handle edge cases
    if dbh <= dbhs[0]:
        return ages[0]
    if dbh >= dbhs[-1]:
        return ages[-1]

    # Linear interpolation
    for i in range(len(dbhs) - 1):
        if dbhs[i] <= dbh <= dbhs[i + 1]:
            dbh_range = dbhs[i + 1] - dbhs[i]
            if dbh_range == 0:
                return ages[i]

            fraction = (dbh - dbhs[i]) / dbh_range
            age_range = ages[i + 1] - ages[i]
            return int(round(ages[i] + fraction * age_range))

    return ages[-1]


def get_years_to_merchantable(
    current_dbh: float,
    species: str,
    merch_dbh: float = 5.0,
    site_index: float = 65.0,
) -> int:
    """
    Calculate years until a tree reaches merchantable size.

    Args:
        current_dbh: Current diameter at breast height (inches)
        species: Species code
        merch_dbh: Merchantable DBH threshold (default 5.0 inches)
        site_index: Site index (default 65)

    Returns:
        Years until merchantable, or 0 if already merchantable
    """
    if current_dbh >= merch_dbh:
        return 0

    current_age = estimate_age_from_dbh(current_dbh, species, site_index)
    merch_age = estimate_age_from_dbh(merch_dbh, species, site_index)

    return max(0, merch_age - current_age)


__all__ = [
    "generate_dbh_age_table",
    "get_age_at_dbh",
    "get_dbh_at_age",
    "estimate_age_from_dbh",
    "get_years_to_merchantable",
    "get_cached_dbh_age_table",
]
