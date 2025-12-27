# API Reference

This section provides detailed documentation for the PyFVS Python API.

## Core Classes

The primary interface for PyFVS simulations:

| Class | Description |
|-------|-------------|
| [`Stand`](stand.md) | Forest stand management - initialization, growth, harvest operations |
| [`Tree`](tree.md) | Individual tree with growth models and attributes |
| [`SimulationEngine`](simulation-engine.md) | High-level simulation orchestration and batch processing |

## Quick Reference

### Stand Initialization

```python
from pyfvs import Stand

# Planted stand (most common)
stand = Stand.initialize_planted(
    trees_per_acre=500,
    site_index=70,
    species='LP',
    ecounit='M231'
)

# Direct constructor
stand = Stand(site_index=70, species='LP')
stand.add_tree(dbh=6.0, height=45.0, age=10)
```

### Growth Simulation

```python
# Grow for 25 years
stand.grow(years=25)

# Grow with custom time step
stand.grow(years=30, time_step=10)
```

### Harvest Operations

```python
# Thin from below to target TPA
stand.thin_from_below(target_tpa=200)

# Thin from above
stand.thin_from_above(target_tpa=300)

# Selection harvest to target basal area
stand.selection_thin(target_basal_area=80)
```

### Output Methods

```python
# Get current metrics
metrics = stand.get_metrics()

# Get yield table (growth history)
yield_table = stand.get_yield_table()

# Get individual tree data
tree_list = stand.get_tree_list()
```

## Return Value Schema

### Stand Metrics Dictionary

The `get_metrics()` method returns:

| Key | Type | Description |
|-----|------|-------------|
| `tpa` | float | Trees per acre |
| `basal_area` | float | Basal area (ft²/acre) |
| `volume` | float | Volume (ft³/acre) |
| `qmd` | float | Quadratic mean diameter (inches) |
| `top_height` | float | Top height / dominant height (feet) |
| `ccf` | float | Crown competition factor |
| `sdi` | float | Stand density index |

### Yield Table DataFrame

The `get_yield_table()` method returns a pandas DataFrame with columns:

| Column | Description |
|--------|-------------|
| `age` | Stand age (years) |
| `tpa` | Trees per acre |
| `basal_area` | Basal area (ft²/acre) |
| `volume` | Volume (ft³/acre) |
| `mean_dbh` | Mean DBH (inches) |
| `mean_height` | Mean height (feet) |
| `qmd` | Quadratic mean diameter (inches) |
| `top_height` | Dominant height (feet) |
| `ccf` | Crown competition factor |
| `sdi` | Stand density index |
| `mortality` | Trees that died this period |
