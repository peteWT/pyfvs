# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Forest Vegetation Simulator (FVS) - Python implementation of the Southern variant for simulating growth and yield of southern yellow pine species (Loblolly, Shortleaf, Longleaf, and Slash pine). Uses object-oriented architecture with modular growth models.

## Key Development Commands

```bash
# Install in development mode with uv
uv pip install -e .

# Run all tests with coverage
uv run pytest

# Run specific test file
uv run pytest tests/test_tree.py

# Run with verbose output
uv run pytest -v --tb=short

# Format code (required before commits)
uv run black src/fvs_python tests

# Lint and type check
uv run flake8 src/fvs_python tests
uv run mypy src/fvs_python

# Run simulation via CLI
uv run fvs-simulate run --years 50 --species LP --site-index 70

# Run main module directly
uv run python -m fvs_python.main
```

## Architecture Overview

### Core Simulation Flow
```
Stand.initialize_planted() → Tree objects → grow() each cycle → Stand.apply_mortality()
```

1. **Tree** (`tree.py`) - Individual tree attributes and growth methods. Uses small-tree model (height-driven) or large-tree model (diameter-driven) based on DBH threshold.
2. **Stand** (`stand.py`) - Manages tree collections using composition pattern. Delegates to specialized component classes.
3. **Growth models** - Modular components in separate files implementing specific FVS equations.
4. **ConfigLoader** (`config_loader.py`) - Unified configuration loading for YAML, TOML, and JSON files with caching.

### Stand Class Composition (Refactored)
The Stand class uses composition with these specialized components:
- **StandMetricsCalculator** (`stand_metrics.py`) - CCF, QMD, SDI, basal area, top height calculations
- **MortalityModel** (`mortality.py`) - Background and density-dependent mortality
- **HarvestManager** (`harvest.py`) - Thinning operations (from below/above, selection, clearcut)
- **CompetitionCalculator** (`competition.py`) - Tree-level competition metrics (PBAL, rank, relative height)
- **StandOutputGenerator** (`stand_output.py`) - Yield tables, tree lists, stock tables, exports

### Growth Model Transition (Critical Logic)
- `DBH < Xmin`: Small-tree model only (height growth drives DBH)
- `Xmin <= DBH <= Xmax`: Weighted blend of both models
- `DBH >= Xmax`: Large-tree model only (DBH growth drives height)
- Default transition: Xmin=1.0", Xmax=3.0" (species-specific in config)

### Key Module Dependencies
```
tree.py
  ├── large_tree_height_growth.py (FVS Section 4.7.2 equations)
  ├── height_diameter.py (Curtis-Arney/Wykoff equations)
  ├── crown_ratio.py (Weibull-based crown ratio)
  ├── bark_ratio.py (Clark 1991 DIB/DOB)
  ├── crown_width.py (forest/open grown equations)
  └── validation.py (parameter bounds checking)

stand.py (composition)
  ├── stand_metrics.py (StandMetricsCalculator)
  ├── mortality.py (MortalityModel)
  ├── harvest.py (HarvestManager)
  ├── competition.py (CompetitionCalculator)
  ├── stand_output.py (StandOutputGenerator)
  └── tree.py

config_loader.py
  └── All coefficient modules use load_coefficient_file() with caching
```

### Configuration System
- Species configs: `/cfg/species/*.yaml` (~80+ species)
- Model coefficients: `/cfg/sn_*.json` (height-diameter, crown width, bark ratio, CCF, etc.)
- Functional forms: `/cfg/functional_forms.yaml` (equation specifications)
- All JSON loading uses `ConfigLoader.load_coefficient_file()` with centralized caching

### Testing
- Unit tests in `/tests/` for individual modules
- Comprehensive integration tests in `test_tree_comprehensive.py`
- Test outputs and plots saved to `/test_output/` for manual verification

## Known Issues

1. **Competition Factor Not Used in Small Tree Model** (`tree.py:_grow_small_tree()`)
   - Small trees don't respond to competition until transitioning to large-tree model

2. **Volume Calculations** (`tree.py:get_volume()`)
   - May not match official FVS exactly; verify against NVEL library if available

3. **Fort Bragg Special Equations** (`cfg/sn_diameter_growth_coefficients.json`)
   - Coefficients for IFOR=20 exist in config but forest ID system not yet implemented

## Recently Fixed

1. **HGMDCR Crown Ratio Modifier** - Fixed to include 100x multiplier per FVS Fortran source (htgf.f)
2. **SDI-Mortality Model** - Implemented full FVS equations 5.0.1-5.0.4 with 55%/85% SDImax thresholds
3. **Chapman-Richards Time Step** - Clarified that equation naturally handles any time step (no fix needed)
4. **Topographic Effects** - Verified slope/aspect correctly applied in CONSPP term
5. **Crown Ratio Time Step Scaling** - Fixed `_update_crown_ratio_weibull()` to scale changes by time_step
6. **Stand.grow() Time Step Handling** - Fixed to respect years parameter instead of forcing 5-year cycles
7. **Mortality Time Step Scaling** - Fixed to pass cycle_length parameter to mortality calculation
8. **Stand Class Decomposition** - Refactored 2000-line Stand class into 5 focused component modules using composition
9. **Tree Height Growth Duplication** - Eliminated duplicate code; `tree.py` now calls `large_tree_height_growth.py` module
10. **Configuration Unification** - All modules now use `ConfigLoader.load_coefficient_file()` with centralized caching
11. **Bark Ratio Path Bug** - Fixed critical bug where `bark_ratio.py` looked in `/docs/` instead of `/cfg/`
12. **Ecological Unit Propagation** - Fixed `tree.grow()` to receive ecounit/forest_type from Stand; previously trees always used species config base_ecounit (232=0.0 for LP), ignoring Stand's ecounit setting
13. **Volume Equations** - Replaced simple form-factor calculation with combined-variable equations (V = a + b × D²H) from Amateis & Burkhart (1987), matching published research with R² > 0.97
14. **Height Growth Cap** - Removed artificial 4.0 ft/5yr cap on POTHTG that was limiting young tree height growth; now uses site-index-based maximum (SI × 0.20)
15. **Relative Height Default** - Fixed relative height (RELHT) to default to 1.0 for codominant trees instead of incorrectly comparing tree height to site index, which was suppressing height growth

## Ecological Unit Effects on Growth

The FVS diameter growth model includes ecological unit adjustments. **Province 232 (Georgia)** is the BASE for loblolly pine (LP), with coefficient 0.0. Other provinces have significant effects:

| Province | Effect on ln(DDS) | Growth Impact |
|----------|-------------------|---------------|
| 232 (Georgia) | 0.0 (base) | 0.24 in/year |
| 231L (Lowland) | +0.256 | 0.30 in/year |
| 255 (Prairie) | +0.275 | 0.31 in/year |
| M222 | +0.582 | 0.36 in/year |
| M231 (Mountain) | +0.790 | **0.43 in/year** |

To achieve growth rates matching typical yield tables (0.3-0.5 in/year), use appropriate ecounit:
```python
stand = Stand(site_index=70, species='LP', ecounit='M231')  # Mountain province
stand = Stand(site_index=70, species='LP', ecounit='231L')  # Section 231 lowland
```

## Validation Status (Manuscript Comparison)

Validation against timber asset account manuscript ("Toward a timber asset account for the United States"):
- **16 of 25 tests pass** - Core mechanics work correctly
- **Current yield ratios: 15-21% of manuscript expectations** (improved from initial 5-10%)
- **Improvements made**:
  - Combined-variable volume equations (Amateis & Burkhart 1987) - validated against published research
  - Fixed height growth cap (was limiting POTHTG to 4 ft/5yr)
  - Fixed relative height calculation (was suppressing codominant tree growth)
  - Ecounit propagation from Stand to Trees (M231 adds +0.790 to growth)
- **Remaining gap**: Tree dimensions at age 25 (DBH=9.4", Ht=40 ft) are below manuscript expectations (DBH~15", Ht~55 ft)
- See `test_output/manuscript_validation/` for validation reports

## Development Priorities

### Code Quality
1. **Consolidate Simulation Functions**: Three overlapping functions in `main.py` (run_simulation, simulate_stand_growth, generate_yield_table)

### Testing & Validation
1. Re-run manuscript validation tests with appropriate ecounit settings
2. Add regression tests with known good outputs
3. Test with large stands (1000+ trees) for performance

### Growth Model Calibration
1. ~~Investigate diameter growth rates~~ **RESOLVED** - Use appropriate ecounit for region
2. Implement NVEL volume equations or find alternative for macOS
3. Consider adding user-settable calibration factor (COR) for fine-tuning
