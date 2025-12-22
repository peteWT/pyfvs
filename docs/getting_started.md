# Getting Started with FVS-Python

This guide walks you through installing FVS-Python and running your first forest growth simulation.

## Prerequisites

- Python 3.8 or higher
- [uv](https://docs.astral.sh/uv/) (recommended) or pip
- Git (for cloning the repository)

## Installation

### Option 1: Install with uv (Recommended)

```bash
# Clone the repository
git clone https://github.com/mihiarc/fvs-python.git
cd fvs-python

# Create virtual environment and install in development mode
uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
uv pip install -e .
```

### Option 2: Install with pip

```bash
git clone https://github.com/mihiarc/fvs-python.git
cd fvs-python

python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
pip install -e .
```

### Verify Installation

```bash
fvs-python --version
```

You should see: `FVS-Python 1.0.0`

## Quick Start

### Running Your First Simulation

The simplest way to run a simulation is via the command line:

```bash
fvs-python simulate --species LP --tpa 500 --site-index 70 --years 30
```

This simulates a Loblolly Pine plantation with:
- 500 trees per acre initial density
- Site index of 70 feet (base age 25)
- 30 years of growth

### Understanding the Output

The simulation produces a yield table showing stand development over time:

```
Age  TPA   Mean DBH  Mean Height  Volume (cuft/acre)
---  ---   --------  -----------  ------------------
  0  500      0.5        1.0              0
  5  498      2.1        8.4            142
 10  495      4.8       22.6            892
 15  488      7.2       38.4          2,156
 20  476      9.1       51.2          3,842
 25  461     10.8       62.4          5,621
 30  442     12.3       71.8          7,312
```

## Using the Python API

For more control, use the Python API directly:

```python
from fvs_python import Stand

# Create a planted stand
stand = Stand.initialize_planted(
    trees_per_acre=500,
    site_index=70,
    species='LP',
    ecounit='M231'  # Mountain province for higher growth rates
)

# Grow for 25 years
stand.grow(years=25)

# Get stand metrics
metrics = stand.get_metrics()
print(f"Trees per acre: {metrics['tpa']:.0f}")
print(f"Mean DBH: {metrics['mean_dbh']:.1f} inches")
print(f"Mean height: {metrics['mean_height']:.1f} feet")
print(f"Volume: {metrics['volume']:.0f} cubic feet/acre")
```

### Working with Individual Trees

```python
from fvs_python import Tree

# Create an individual tree
tree = Tree(
    species='LP',
    dbh=6.0,
    height=45.0,
    age=15,
    crown_ratio=0.45
)

# Grow the tree for 5 years
tree.grow(
    site_index=70,
    basal_area=120.0,
    time_step=5
)

print(f"New DBH: {tree.dbh:.1f} inches")
print(f"New height: {tree.height:.1f} feet")
print(f"Volume: {tree.get_volume():.1f} cubic feet")
```

## Supported Species

FVS-Python supports four southern yellow pine species:

| Code | Common Name | Scientific Name |
|------|-------------|-----------------|
| LP | Loblolly Pine | *Pinus taeda* |
| SP | Shortleaf Pine | *Pinus echinata* |
| LL | Longleaf Pine | *Pinus palustris* |
| SA | Slash Pine | *Pinus elliottii* |

List available species with:

```bash
fvs-python list-species --detailed
```

## Key Parameters

### Site Index

Site index (SI) measures site productivity - the expected height of dominant trees at a base age (25 years for southern pines). Common ranges:

| Site Quality | Site Index (ft) |
|--------------|-----------------|
| Low | 50-60 |
| Medium | 60-75 |
| High | 75-90 |

### Ecological Units (Ecounit)

Ecounits affect growth rates significantly. Province 232 (Georgia) is the base for Loblolly Pine:

| Ecounit | Region | Growth Effect |
|---------|--------|---------------|
| 232 | Georgia (base) | No modifier |
| 231L | Lowland | +6% diameter growth |
| 255 | Prairie | +6% diameter growth |
| M231 | Mountain | **+54% diameter growth** |

For realistic growth rates, use an appropriate ecounit:

```python
stand = Stand.initialize_planted(
    trees_per_acre=500,
    site_index=70,
    species='LP',
    ecounit='M231'  # Highest growth rates
)
```

## Common Simulation Scenarios

### Unthinned Plantation

```python
from fvs_python import Stand

# High-density unthinned for maximum volume
stand = Stand.initialize_planted(
    trees_per_acre=800,
    site_index=55,
    species='LP',
    ecounit='M231'
)
stand.grow(years=25)
metrics = stand.get_metrics()
print(f"Final volume: {metrics['volume']:.0f} cuft/acre")
```

### Thinned Plantation

```python
from fvs_python import Stand

# Start dense, thin to grow larger trees
stand = Stand.initialize_planted(
    trees_per_acre=500,
    site_index=70,
    species='LP',
    ecounit='M231'
)

# Grow 10 years
stand.grow(years=10)

# Thin from below to 75 trees per acre
stand.thin_from_below(target_tpa=75)

# Grow 15 more years
stand.grow(years=15)

metrics = stand.get_metrics()
print(f"Final TPA: {metrics['tpa']:.0f}")
print(f"Mean DBH: {metrics['mean_dbh']:.1f} inches")  # Larger trees!
```

### Generating Yield Tables

```bash
fvs-python yield-table \
    --species LP SP \
    --site-indices 60 70 80 \
    --densities 300 500 700 \
    --years 50
```

Or with Python:

```python
from fvs_python import SimulationEngine
from pathlib import Path

engine = SimulationEngine(Path('./output'))

yield_table = engine.simulate_yield_table(
    species=['LP', 'SP'],
    site_indices=[60.0, 70.0, 80.0],
    planting_densities=[300, 500, 700],
    years=50
)

yield_table.to_csv('yield_table.csv', index=False)
```

## Configuration

### Viewing Species Configuration

```bash
fvs-python show-config --species LP --format yaml
```

### Configuration Files

Configuration files are in the `cfg/` directory:

```
cfg/
├── species/           # Species-specific parameters (90 files)
│   ├── lp_loblolly_pine.yaml
│   ├── sp_shortleaf_pine.yaml
│   └── ...
├── functional_forms.yaml  # Growth equation specifications
├── sn_diameter_growth_coefficients.json
├── sn_height_diameter_coefficients.json
└── ...
```

### Validating Configuration

```bash
fvs-python validate-config
```

## Output Files

By default, simulations save results to `./output/`:

```
output/
├── simulation_results.csv    # Yield table data
├── fvs-python.log           # Log file
└── plots/                   # Generated charts (if enabled)
```

Disable file output with `--no-save`:

```bash
fvs-python simulate --species LP --tpa 500 --site-index 70 --no-save
```

## Troubleshooting

### Common Issues

**Import Error: Module not found**

Make sure you installed in development mode:
```bash
uv pip install -e .
```

**Low growth rates**

Check your ecounit setting. Province 232 (base) has lower growth than other provinces:
```python
stand = Stand.initialize_planted(..., ecounit='M231')  # Higher growth
```

**Volume seems wrong**

FVS-Python uses combined-variable volume equations (Amateis & Burkhart 1987). Values are for merchantable stem volume, not total tree volume.

### Getting Help

- Check the README.md in the project root for architecture overview
- See CLAUDE.md for known issues and recent fixes
- Review the [Validation Specification](FVS_PYTHON_VALIDATION_SPEC.md) for model accuracy

## Next Steps

- Read the Architecture Overview in README.md to understand the codebase
- Explore the SINGLE_TREE_GROWTH_GUIDE.md for equation details
- Check the [API Reference](api/index.rst) for complete module documentation
- Review test outputs in `test_output/` for validation results
