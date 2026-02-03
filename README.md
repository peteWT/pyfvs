<div align="center">
  <h1>pyFVS</h1>

  <p><strong>Forest growth modeling in Python</strong></p>

  <p>
    <a href="https://pypi.org/project/pyfvs/"><img src="https://img.shields.io/pypi/v/pyfvs?color=006D6D&label=PyPI" alt="PyPI"></a>
    <a href="https://mihiarc.github.io/pyfvs/"><img src="https://img.shields.io/badge/docs-GitHub%20Pages-006D6D" alt="Documentation"></a>
    <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-006D6D" alt="License: MIT"></a>
    <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/python-3.9+-006D6D" alt="Python 3.9+"></a>
  </p>
</div>

---

A Python implementation of the Forest Vegetation Simulator (FVS) supporting multiple regional variants. Simulate individual tree growth and yield across the United States.

## Supported Variants

| Variant | Region | Species Count | Key Species |
|---------|--------|---------------|-------------|
| **SN** | Southern US | 90 | Loblolly Pine, Shortleaf Pine, Longleaf Pine, Slash Pine |
| **LS** | Lake States (MI, WI, MN) | 67 | Red Pine, Jack Pine, Sugar Maple |
| **PN** | Pacific Northwest Coast | 39 | Douglas-fir, Western Hemlock, Sitka Spruce |
| **WC** | West Cascades (OR, WA) | 37 | Douglas-fir, Western Hemlock, Western Red Cedar |
| **NE** | Northeast (13 states) | 108 | Red Maple, Sugar Maple, Northern Red Oak |
| **CS** | Central States (IL, IN, IA, MO) | 96 | White Oak, Black Walnut, Yellow-Poplar |
| **OP** | ORGANON Pacific Northwest | 18 | Douglas-fir, Red Alder |

## Southern Yellow Pines (SN Variant)

| Code | Species | Scientific Name |
|------|---------|-----------------|
| LP | Loblolly Pine | *Pinus taeda* |
| SP | Shortleaf Pine | *Pinus echinata* |
| LL | Longleaf Pine | *Pinus palustris* |
| SA | Slash Pine | *Pinus elliottii* |

## Quick Start

```bash
pip install pyfvs
```

```python
from pyfvs import Stand

# Initialize a planted stand
stand = Stand.initialize_planted(
    species="LP",
    trees_per_acre=500,
    site_index=70
)

# Simulate 50 years of growth
stand.grow(years=50)

# Get results
metrics = stand.get_metrics()
print(f"Final volume: {metrics['volume']:.0f} ft³/acre")
```

## Growth Models

pyFVS implements individual tree growth models from FVS documentation:

### Height-Diameter (Curtis-Arney)
```
height = 4.5 + p2 × exp(-p3 × DBH^p4)
```

### Large Tree Diameter Growth
```
ln(DDS) = β₁ + β₂×ln(DBH) + β₃×DBH² + β₄×ln(CR) + β₅×RH + β₆×SI + ...
```

### Small Tree Height Growth (Chapman-Richards)
```
height = c1 × SI^c2 × (1 - exp(c3 × age))^(c4 × SI^c5)
```

## Architecture

```
Tree Initial State
       │
       ▼
   DBH >= 3.0? ──No──► Small Tree Model
       │                    │
      Yes                   ▼
       │            Height Growth
       ▼                    │
  Large Tree Model          ▼
       │            Height-Diameter
       ▼                    │
  Predict ln(DDS)           ▼
       │              Update DBH
       ▼                    │
  Calculate DBH Growth      │
       │                    │
       └────────┬───────────┘
                ▼
        Update Crown Ratio
                │
                ▼
        Crown Competition
                │
                ▼
        Updated Tree State
```

## Configuration

Species parameters are stored in YAML configuration files:

```yaml
# cfg/species/lp_loblolly_pine.yaml
species_code: "LP"
common_name: "Loblolly Pine"
height_diameter:
  p2: 243.860648
  p3: 4.28460566
  p4: -0.47130185
bark_ratio:
  b1: -0.48140
  b2: 0.91413
```

## Output

pyFVS generates yield tables with standard forest metrics:

| Age | TPA | QMD | Height | BA | Volume |
|-----|-----|-----|--------|-----|--------|
| 0 | 500 | 0.5 | 1.0 | 0.7 | 0 |
| 10 | 485 | 4.2 | 28.5 | 47.2 | 892 |
| 20 | 420 | 7.8 | 52.1 | 139.8 | 3,241 |
| 30 | 310 | 10.4 | 68.3 | 182.5 | 5,128 |
| ... | ... | ... | ... | ... | ... |

## References

- [FVS Southern Variant Documentation](https://www.fs.usda.gov/fmsc/fvs/)
- Bechtold & Patterson (2005) "The Enhanced Forest Inventory and Analysis Program"

## Citation

```bibtex
@software{pyfvs2025,
  title = {pyFVS: Python Implementation of the Forest Vegetation Simulator},
  author = {Mihiar, Christopher},
  year = {2025},
  url = {https://github.com/mihiarc/pyfvs}
}
```

---

<div align="center">
  <sub>Built by <a href="https://github.com/mihiarc">Chris Mihiar</a></sub>
</div>
