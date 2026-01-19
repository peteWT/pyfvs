# FVS Variant Coverage Analysis

This document tracks PyFVS implementation status compared to the official USDA Forest Service FVS variants.

## Official FVS Variants (22 Total)

The Forest Vegetation Simulator (FVS) is maintained by the USDA Forest Service Forest Management Service Center (FMSC) in Fort Collins, Colorado. There are 22 geographic variants, each calibrated to specific forest regions.

### Eastern Variants (4)

| Code | Variant Name | States/Region | Species | Base Cycle |
|------|-------------|---------------|---------|------------|
| **CS** | Central States | IL, IN, IA, MO | 91 | 10 years |
| **LS** | Lake States | MI, WI, MN | 67 | 10 years |
| **NE** | Northeast | CT, DE, MA, MD, ME, NH, NJ, NY, OH, PA, RI, VT, WV | 105 | 10 years |
| **SN** | Southern | AL, AR, FL, GA, KY, LA, MS, NC, OK, SC, TN, TX, VA | 90 | 5 years |

### Western Variants (18)

| Code | Variant Name | Region | Species | Base Cycle |
|------|-------------|--------|---------|------------|
| **AK** | Alaska | SE Alaska, Coastal BC | 21 | 10 years |
| **BM** | Blue Mountains | OR/WA Blue Mountains | ~30 | 10 years |
| **CA** | Inland California | CA interior, Southern Cascades | ~40 | 10 years |
| **CI** | Central Idaho | Central Idaho | ~25 | 10 years |
| **CR** | Central Rockies | CO, WY Rocky Mountains | ~25 | 10 years |
| **EC** | East Cascades | OR/WA East Cascades | ~35 | 10 years |
| **EM** | Eastern Montana | Eastern Montana | ~25 | 10 years |
| **IE** | Inland Empire | Northern Idaho/Inland Empire | ~30 | 10 years |
| **KT** | Kootenai-Kaniksu-Tally Lake | MT/ID border region | ~25 | 10 years |
| **NC** | Klamath Mountains | NW California/SW Oregon | ~40 | 10 years |
| **OC** | ORGANON Southwest | Southwest Oregon | ~30 | variable |
| **OP** | ORGANON Pacific NW | Oregon/Washington | ~25 | variable |
| **PN** | Pacific Northwest Coast | WA, OR, N CA coast | 39 | 10 years |
| **SO** | South Central OR/NE CA | South OR, NE California | 31 | 10 years |
| **TT** | Tetons | Teton region (WY/ID/MT) | 16 | 10 years |
| **UT** | Utah | Utah | 22 | 10 years |
| **WC** | Westside Cascades | Western OR/WA Cascades | 37 | 10 years |
| **WS** | Western Sierra Nevada | California Sierra Nevada | 41 | 10 years |

---

## PyFVS Implementation Status

### Implemented Variants ✅

| Code | Name | Species | DDS Equation | Key Features | Status |
|------|------|---------|--------------|--------------|--------|
| **SN** | Southern | 90 | ln(DDS) = f(D, CR, RELHT, SI, BA, ecounit) | Chapman-Richards height, ecounit modifiers, 5-year cycle | Complete |
| **LS** | Lake States | 67 | ln(DDS) = f(D, CR, RELDBH, SI, BA, BAL) | Linear DDS, Curtis-Arney H-D, RELDBH competition | Complete |
| **PN** | Pacific NW Coast | 39 | ln(DDS) = f(D, CR, RELHT, SI, BA, BAL, elev, slope, aspect) | Topographic effects, high SI (100-200) | Complete |
| **WC** | West Cascades | 37 | ln(DDS) = f(D, CR, RELHT, SI, BA, BAL, elev, slope, aspect) | Shares PN equation structure, forest intercepts | Complete |
| **NE** | Northeast | 108 | BA Growth = B1 * SI * (1 - exp(-B2 * DBH)) | NE-TWIGS iterative BA growth, Curtis-Arney H-D, 10-year cycle | Complete |
| **CS** | Central States | 96 | ln(DDS) = f(D, CR, RELDBH, SI, BA, BAL) | Same equation form as LS, Curtis-Arney H-D, 10-year cycle | Complete |
| **CA** | Inland California | 50 | ln(DDS) = f(D, CR, RELHT, SI, BA, BAL, PCCF, elev, slope, aspect) | 13 equation sets, special RW/GS equation, topographic effects | Complete |

### In Progress 🔄

None currently.

### Planned ⏳

| Code | Name | Species | Priority | Notes |
|------|------|---------|----------|-------|
| **WS** | Western Sierra Nevada | 41 | Medium | California Sierra coverage |
| **CR** | Central Rockies | ~25 | Medium | Colorado/Wyoming coverage |

### Not Started ❌

| Code | Name | Priority | Rationale |
|------|------|----------|-----------|
| AK | Alaska | Low | Specialized subarctic species |
| BM | Blue Mountains | Low | Overlap with PN/WC/EC |
| CI | Central Idaho | Low | Overlap with IE |
| EC | East Cascades | Low | Overlap with WC |
| EM | Eastern Montana | Low | Specialized region |
| IE | Inland Empire | Low | Overlap with other western |
| KT | Kootenai-Kaniksu | Low | Specialized region |
| NC | Klamath Mountains | Medium | Unique CA/OR interface |
| OC | ORGANON Southwest | Low | ORGANON-based (different model) |
| OP | ORGANON Pacific NW | Low | ORGANON-based (different model) |
| SO | South Central OR/NE CA | Low | Small region (31 species) |
| TT | Tetons | Low | Small region (16 species) |
| UT | Utah | Low | Small region (22 species) |

---

## Coverage Statistics

| Metric | Official FVS | PyFVS | Coverage |
|--------|-------------|-------|----------|
| Total Variants | 22 | 7 | 32% |
| Eastern Variants | 4 | 4 (SN, LS, NE, CS) | 100% |
| Western Variants | 18 | 3 (PN, WC, CA) | 17% |
| Species (deduplicated est.) | ~400 unique | 487 | ~100% |
| US Forest Area Coverage | 100% | ~75% | By commercial timber volume |

---

## Variant Technical Specifications

### Diameter Growth Equation Forms

**Type 1: ln(DDS) with Ecological Units (SN)**
```
ln(DDS) = CONSPP + b1*ln(DBH) + b2*DBH² + b3*CR + b4*CR² + b5*RELHT
        + b6*SI + b7*BA + ECOUNIT_MODIFIER
```

**Type 2: Linear DDS with RELDBH (LS)**
```
ln(DDS) = CONSPP + INTERC + b1/DBH + b2*DBH + b3*DBH² + b4*RELDBH
        + b5*RELDBH² + b6*CR + b7*CR² + b8*BA + b9*BAL + b10*SI
```

**Type 3: ln(DDS) with Topographic Effects (PN, WC)**
```
ln(DDS) = CONSPP + b1*ln(DBH) + b2*DBH² + b3*CR + b4*CR² + b5*RELHT
        + b6*SI + b7*BA + b8*BAL + b9*ELEV + b10*SLOPE + b11*COS(ASPECT)
        + b12*SIN(ASPECT) + b13*SLOPE*COS(ASPECT)
```

**Type 4: NE-TWIGS Basal Area Growth (NE)**
```
POTBAG = B1 * SI * (1 - exp(-B2 * DBH))
Adjusted Growth = POTBAG * 0.7 * BAGMOD * CR_MODIFIER
Iterates annually (10x for 10-year cycle)
```

**Type 5: CS Central States (same as LS)**
```
ln(DDS) = INTERC + VDBHC/D + DBHC*D + DBH2C*D² + RDBHC*RELDBH
        + RDBHSQC*RELDBH² + CRWNC*CR + CRSQC*CR² + SBAC*BA + BALC*BAL + SITEC*SI
DDS = exp(ln(DDS))
```

**Type 6: CA Inland California (13 equation sets)**
```
ln(DDS) = CONSPP + DGLD*ln(D) + CR*(DGCR + CR*DGCRSQ) + DGDSQ*D²
        + DGDBAL*BAL/ln(D+1) + DGPCCF*PCCF + DGHAH*RELHT + DGLBA*ln(BA)
        + DGBAL*BAL + DGSITE*ln(SI) + DGEL*ELEV + DGELSQ*ELEV²
        + DGSLOP*SLOPE + DGCASP*SLOPE*COS(ASP) + DGSASP*SLOPE*SIN(ASP)
Special: Giant Sequoia/Redwood use ln(DDS) = -3.502444 + 0.415435*ln(SI)
Special: Tanoak uses 5-year model scaled to 10-year
```

### Height-Diameter Relationships

| Variant | Equation Type | Key Parameters |
|---------|---------------|----------------|
| SN | Wykoff (1982) | P2, P3 species-specific |
| LS | Curtis-Arney | P2, P3, P4 with CCF adjustment |
| PN | Curtis-Arney | P2, P3, P4 species-specific |
| WC | Curtis-Arney | P2, P3, P4 species-specific |
| NE | Curtis-Arney | P2, P3, P4 species-specific |
| CS | Curtis-Arney | P2, P3, P4 species-specific |
| CA | Curtis-Arney | P2, P3, P4, Z (breakpoint) species-specific |

### Crown Ratio Models

All variants use Weibull-based crown ratio prediction with species-specific coefficients:
```
CR = 1 - exp(-A * (X^B))
where X = relative position in stand
```

---

## Implementation Roadmap

### Phase 1: Eastern Variants (Complete) ✅
1. ✅ SN - Southern (complete)
2. ✅ LS - Lake States (complete)
3. ✅ NE - Northeast (complete)
4. ✅ CS - Central States (complete)

### Phase 2: California Variants (In Progress)
5. ✅ CA - Inland California (complete)
6. WS - Western Sierra Nevada
7. NC - Klamath Mountains (if needed)

### Phase 3: Rocky Mountain Variants
8. CR - Central Rockies
9. UT - Utah
10. TT - Tetons

### Phase 4: Remaining Western (as needed)
11. BM, CI, EC, EM, IE, KT, SO
12. AK - Alaska (specialized)

---

## Official Documentation Sources

- [FVS Software Downloads](https://www.fs.usda.gov/managing-land/forest-management/fvs/software)
- [FVS Variant Key](https://www.fs.usda.gov/fvs/software/variantkey.shtml)
- [FVS User Guides](https://www.fs.usda.gov/managing-land/forest-management/fvs/documents/guides)
- [FVS GitHub Repository](https://github.com/USDAForestService/ForestVegetationSimulator)
- [Species Crosswalk (Excel)](https://www.fs.usda.gov/fmsc/ftp/fvs/docs/overviews/FVS_SpeciesCrosswalk.pdf)

### Variant Overview PDFs
Each variant has an overview document at:
```
https://www.fs.usda.gov/sites/default/files/forest-management/fvs-{code}-overview.pdf
```
Replace `{code}` with lowercase variant code (e.g., `sn`, `ne`, `ls`).

---

## Adding a New Variant

### Required Files

```
src/pyfvs/cfg/{variant}/
├── {variant}_diameter_growth_coefficients.json
├── {variant}_height_diameter_coefficients.json
├── {variant}_mortality_coefficients.json (optional)
├── {variant}_species_config.yaml
└── species/
    └── {species_code}.yaml (for each species)
```

### Code Changes

1. **config_loader.py**: Add variant to `SUPPORTED_VARIANTS` dict
2. **tree.py**: Add `_grow_large_tree_{variant}()` method if equation differs
3. **height_diameter.py**: Add variant to `VARIANT_COEFFICIENT_FILES`
4. **{variant}_diameter_growth.py**: Create if equation structure differs from existing

### Coefficient Sources

1. Official FVS Fortran source: https://github.com/USDAForestService/ForestVegetationSimulator
2. Variant overview PDFs (tables of coefficients)
3. Original research papers cited in overviews

---

## Notes

- Species counts include "other softwood" and "other hardwood" composite categories
- Some species appear in multiple variants with different coefficients
- ORGANON variants (OC, OP) use a different underlying model architecture
- Canadian variants (BC, ON) are maintained externally and not included

Last updated: 2026-01-19
