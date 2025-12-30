# FVS-Python Manuscript Validation Report

## Overview

This report validates fvs-python yield predictions against the
timber asset account manuscript data.

### Source
- **Manuscript**: 'Toward a timber asset account for the United States'
- **Authors**: Bruck, Mihiar, Mei, Brandeis, Chambers, Hass, Wentland, Warziniack
- **FVS Version**: FS2025.1

## Species Simulated

- **LP** (Loblolly Pine): SI=55 (North), SI=65 (South)
- **SA** (Slash Pine): SI=55 (North), SI=65 (South)
- **SP** (Shortleaf Pine): SI=55 (North), SI=65 (South)
- **LL** (Longleaf Pine): SI=55 (North), SI=65 (South)

## Summary Statistics

### Yields at Age 25

| Species   | Region   |   Site_Index |   Volume_Tons |   Mean_DBH |   TPA |
|:----------|:---------|-------------:|--------------:|-----------:|------:|
| LP        | North    |           55 |          73.9 |       7.81 |   451 |
| LP        | South    |           65 |          98.4 |       8.23 |   464 |
| SA        | North    |           55 |          62.9 |       6.73 |   462 |
| SA        | South    |           65 |          85.4 |       7.28 |   462 |
| SP        | North    |           55 |          71.3 |       7.18 |   447 |
| SP        | South    |           65 |          95.7 |       7.64 |   456 |
| LL        | North    |           55 |          54.8 |       6.36 |   458 |
| LL        | South    |           65 |          75.9 |       6.9  |   466 |

### Files Generated

- `full_yield_simulation.csv`: Complete simulation results
- `yield_summary_by_species_age.csv`: Summary by species and age
- `table1_full_comparison.csv`: Table 1 validation details
- `lev_mai_comparison.csv`: LEV vs MAI rotation ages
