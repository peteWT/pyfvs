# PyFVS Documentation

[![PyPI version](https://badge.fury.io/py/pyfvs-fia.svg)](https://badge.fury.io/py/pyfvs-fia)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**PyFVS** is a Python implementation of the Forest Vegetation Simulator (FVS) Southern variant for simulating growth and yield of southern yellow pine species.

Part of the **FIA Python Ecosystem**:

- [PyFIA](https://github.com/mihiarc/pyfia) - Survey/plot data analysis
- [GridFIA](https://github.com/mihiarc/gridfia) - Spatial raster analysis
- **PyFVS** - Growth/yield simulation (this package)
- [AskFIA](https://github.com/mihiarc/askfia) - AI conversational interface

## Quick Example

```python
from pyfvs import Stand

# Create a planted loblolly pine stand
stand = Stand.initialize_planted(
    trees_per_acre=500,
    site_index=70,
    species='LP',
    ecounit='M231'  # Mountain province for realistic growth
)

# Simulate 25 years of growth
stand.grow(years=25)

# Get stand metrics
metrics = stand.get_metrics()
print(f"Volume: {metrics['volume']:.0f} ft³/acre")
print(f"Trees/acre: {metrics['tpa']:.0f}")
print(f"Mean DBH: {metrics['qmd']:.1f} inches")
```

## Documentation Overview

<div class="grid cards" markdown>

-   :material-rocket-launch: **[Getting Started](getting-started.md)**

    Installation and your first simulation

-   :material-tree: **[Stand Class](api/stand.md)**

    Manage forest stands and growth cycles

-   :material-pine-tree: **[Tree Class](api/tree.md)**

    Individual tree growth models

-   :material-cog: **[SimulationEngine](api/simulation-engine.md)**

    High-level simulation orchestration

</div>

## Features

| Feature | Description |
|---------|-------------|
| **Individual Tree Models** | Small-tree (height-driven) and large-tree (DBH-driven) growth |
| **Stand Dynamics** | Mortality, competition, crown competition factor |
| **Harvest Operations** | Thin from below/above, selection harvest, clearcut |
| **Volume Calculations** | Combined-variable equations (Amateis & Burkhart 1987) |
| **Ecological Effects** | Regional growth modifiers by ecological province |
| **Multi-format Output** | CSV, JSON, Excel exports with yield tables |

## Supported Species

| Code | Species | Scientific Name |
|------|---------|-----------------|
| LP | Loblolly Pine | *Pinus taeda* |
| SP | Shortleaf Pine | *Pinus echinata* |
| SA | Slash Pine | *Pinus elliottii* |
| LL | Longleaf Pine | *Pinus palustris* |

## Installation

```bash
pip install pyfvs-fia
```

Or with uv:

```bash
uv add pyfvs-fia
```
