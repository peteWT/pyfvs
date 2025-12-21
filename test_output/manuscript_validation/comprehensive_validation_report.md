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
| LP        | North    |           55 |          40   |       6.49 |   455 |
| LP        | South    |           65 |          49.9 |       6.9  |   447 |
| SA        | North    |           55 |          31.8 |       5.51 |   445 |
| SA        | South    |           65 |          42.8 |       5.96 |   455 |
| SP        | North    |           55 |          39.6 |       6.03 |   461 |
| SP        | South    |           65 |          51.6 |       6.5  |   458 |
| LL        | North    |           55 |          29.6 |       5.3  |   461 |
| LL        | South    |           65 |          39   |       5.77 |   459 |

### Files Generated

- `full_yield_simulation.csv`: Complete simulation results
- `yield_summary_by_species_age.csv`: Summary by species and age
- `table1_full_comparison.csv`: Table 1 validation details
- `lev_mai_comparison.csv`: LEV vs MAI rotation ages
