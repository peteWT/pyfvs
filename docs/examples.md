# Examples

Practical examples demonstrating PyFVS capabilities.

## Basic Simulation

### Simple Stand Growth

```python
from pyfvs import Stand

# Create a planted loblolly pine stand
stand = Stand.initialize_planted(
    trees_per_acre=500,
    site_index=70,
    species='LP'
)

# Grow for 30 years
stand.grow(years=30)

# Print final metrics
metrics = stand.get_metrics()
print(f"Age: {stand.age} years")
print(f"Trees/acre: {metrics['tpa']:.0f}")
print(f"Volume: {metrics['volume']:.0f} ft³/acre")
print(f"Mean DBH: {metrics['qmd']:.1f} inches")
```

## Thinning Scenarios

### Single Thin

```python
from pyfvs import Stand

stand = Stand.initialize_planted(
    trees_per_acre=700,
    site_index=70,
    species='LP',
    ecounit='M231'
)

# Grow to thinning age
stand.grow(years=15)
print(f"Pre-thin: {stand.get_metrics()['tpa']:.0f} TPA")

# Thin from below to 300 TPA
stand.thin_from_below(target_tpa=300)
print(f"Post-thin: {stand.get_metrics()['tpa']:.0f} TPA")

# Continue to harvest
stand.grow(years=15)
print(f"Final volume: {stand.get_metrics()['volume']:.0f} ft³/acre")
```

### Multiple Thins

```python
from pyfvs import Stand

stand = Stand.initialize_planted(
    trees_per_acre=800,
    site_index=75,
    species='LP',
    ecounit='M231'
)

# First thin at age 12
stand.grow(years=12)
stand.thin_from_below(target_tpa=400)

# Second thin at age 20
stand.grow(years=8)
stand.thin_from_below(target_tpa=200)

# Final harvest at age 35
stand.grow(years=15)

metrics = stand.get_metrics()
print(f"Final: {metrics['volume']:.0f} ft³/acre, {metrics['qmd']:.1f}\" DBH")
```

## Yield Table Generation

### Single Species

```python
from pyfvs import Stand
import pandas as pd

results = []

for si in [60, 70, 80]:
    for tpa in [400, 600]:
        stand = Stand.initialize_planted(
            trees_per_acre=tpa,
            site_index=si,
            species='LP',
            ecounit='M231'
        )
        stand.grow(years=30)

        m = stand.get_metrics()
        results.append({
            'site_index': si,
            'initial_tpa': tpa,
            'final_volume': m['volume'],
            'final_tpa': m['tpa'],
            'final_qmd': m['qmd']
        })

df = pd.DataFrame(results)
print(df.to_string(index=False))
```

### Using SimulationEngine

```python
from pyfvs.simulation_engine import SimulationEngine
from pathlib import Path

engine = SimulationEngine(output_dir=Path('./output'))

yield_table = engine.simulate_yield_table(
    species=['LP', 'SP'],
    site_indices=[60, 70, 80],
    planting_densities=[400, 600],
    years=40
)

print(yield_table.groupby(['species', 'site_index', 'initial_tpa']).last())
```

## Species Comparison

```python
from pyfvs import Stand

species_list = ['LP', 'SP', 'SA', 'LL']

print("Volume at age 30 (SI=70, 500 TPA):")
print("-" * 40)

for sp in species_list:
    stand = Stand.initialize_planted(
        trees_per_acre=500,
        site_index=70,
        species=sp,
        ecounit='M231'
    )
    stand.grow(years=30)
    m = stand.get_metrics()
    print(f"{sp}: {m['volume']:,.0f} ft³/acre")
```

## Export Results

### To CSV

```python
from pyfvs import Stand

stand = Stand.initialize_planted(trees_per_acre=500, site_index=70, species='LP')
yield_table = stand.grow(years=50)

# Export yield table
yield_table.to_csv('yield_table.csv', index=False)

# Export tree list
tree_list = stand.get_tree_list_dataframe()
tree_list.to_csv('tree_list.csv', index=False)
```

### To Excel

```python
import pandas as pd

with pd.ExcelWriter('simulation_results.xlsx') as writer:
    yield_table.to_excel(writer, sheet_name='Yield Table', index=False)
    tree_list.to_excel(writer, sheet_name='Tree List', index=False)
```
