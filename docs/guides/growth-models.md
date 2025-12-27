# Growth Models

PyFVS implements the FVS Southern variant growth models for individual tree simulation.

## Model Overview

Trees use different growth models based on size:

| DBH Range | Model | Driver |
|-----------|-------|--------|
| < 1.0" | Small-tree | Height growth drives DBH |
| 1.0" - 3.0" | Blended | Smooth transition |
| > 3.0" | Large-tree | DBH growth drives height |

## Small-Tree Model

For seedlings and saplings (DBH < 1.0"), height growth follows the Chapman-Richards equation:

```
H(t) = c1 × SI^c2 × (1 - exp(c3 × t))^(c4 × SI^c5)
```

DBH is then derived from the height-diameter relationship.

### Coefficients (Loblolly Pine)

| Parameter | Value | Description |
|-----------|-------|-------------|
| c1 | 1.1421 | Scale factor |
| c2 | 1.0042 | Site index exponent |
| c3 | -0.0374 | Rate parameter |
| c4 | 0.7632 | Asymptote shape |
| c5 | 0.0358 | SI effect on shape |

## Large-Tree Model

For established trees (DBH > 3.0"), diameter growth uses the DDS (diameter squared increment) equation:

```
ln(DDS) = β₁ + β₂×ln(DBH) + β₃×DBH² + β₄×ln(CR) + β₅×RH
        + β₆×SI + β₇×BA + β₈×PBAL + topographic terms
```

### Variables

| Variable | Description |
|----------|-------------|
| DBH | Diameter at breast height (inches) |
| CR | Crown ratio (0-1) |
| RH | Relative height (tree height / dominant height) |
| SI | Site index |
| BA | Stand basal area (ft²/acre) |
| PBAL | Basal area in larger trees (ft²/acre) |

### Height Growth

Large-tree height growth is calculated as:

```
HTG = POTHTG × (0.25 × HGMDCR + 0.75 × HGMDRH)
```

Where:
- POTHTG = Potential height growth from site curves
- HGMDCR = Crown ratio modifier
- HGMDRH = Relative height modifier

## Transition Zone

Trees between 1.0" and 3.0" DBH use a blended model:

```python
# Smoothstep blending
t = (dbh - 1.0) / (3.0 - 1.0)
weight = 3*t² - 2*t³  # Smooth 0→1 transition

growth = (1 - weight) * small_tree_growth + weight * large_tree_growth
```

This prevents discontinuities in growth predictions.

## Time Step Handling

FVS was calibrated for 5-year growth cycles. PyFVS maintains this internally:

- Cycles > 5 years are subdivided into 5-year sub-cycles
- Competition is recalculated at each 5-year interval
- This ensures consistent dynamics regardless of user-specified time step

```python
# These produce equivalent results
stand.grow(years=25, time_step=5)   # 5 cycles of 5 years
stand.grow(years=25, time_step=25)  # Internally uses 5-year sub-cycles
```
