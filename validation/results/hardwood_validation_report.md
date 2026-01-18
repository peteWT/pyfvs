# PyFVS Hardwood Species Validation Report

**Generated:** 2025-01-18
**Species Validated:** Yellow-poplar (YP), White Oak (WO), Sweetgum (SU)

## Executive Summary

Validation testing compared PyFVS hardwood growth predictions against published yield tables from Silvics of North America and Schnur (1937). Significant discrepancies were found, primarily due to:

1. **Site index base age differences** - Published hardwood tables use base age 50, while FVS-SN uses different internal equations
2. **Height over-prediction** - Models consistently over-predict height by 25-60%
3. **DBH under-prediction** - Diameter growth is under-predicted by 5-40%

## Species-Specific Results

### Yellow-poplar (YP)

**Source:** Beck & Della-Bianca 1970, Silvics of North America

| Site Index | Mean Height Error | Mean DBH Error | Mean Volume Error |
|------------|------------------|----------------|-------------------|
| SI 82      | 42.2% (over)     | 24.8% (under)  | 323% (over)       |
| SI 98      | 33.6% (over)     | 32.7% (under)  | 221% (over)       |
| SI 115     | 24.1% (over)     | 41.0% (under)  | 199% (over)       |

**Key Observations:**
- Height is systematically over-predicted
- The over-prediction decreases on better sites
- DBH is consistently under-predicted
- Volume errors compound height errors (V ∝ D²H)

**Example (SI 82, Age 50):**
- Simulated: Height 115 ft, DBH 10.0"
- Reference: Height 82 ft, DBH 13.4"

### White Oak / Upland Oak (WO)

**Source:** Schnur 1937, Carmean 1971

| Site Index | Mean Height Error | Mean DBH Error | Mean Volume Error |
|------------|------------------|----------------|-------------------|
| SI 55      | 61.1% (over)     | 5.2% (under)   | 179% (over)       |
| SI 65      | 59.5% (over)     | 10.9% (under)  | 152% (over)       |
| SI 75      | 58.3% (over)     | 16.9% (under)  | 136% (over)       |

**Key Observations:**
- Height over-prediction is severe (~60%)
- DBH predictions are actually quite close (5-17% error)
- Likely due to site index base age mismatch (oak SI uses base 50)

**Example (SI 65, Age 50):**
- Simulated: Height 95.6 ft, DBH 6.9"
- Reference: Height 65 ft, DBH 7.9"

### Sweetgum (SU)

**Source:** Silvics of North America (Mississippi Delta data)

| Site Index | Mean Height Error | Mean DBH Error | Mean Volume Error |
|------------|------------------|----------------|-------------------|
| SI 70      | 52.0% (over)     | 18.5% (under)  | 59% (over)        |
| SI 80      | 51.3% (over)     | 23.6% (under)  | 48% (over)        |
| SI 90      | 46.8% (over)     | 27.8% (under)  | 40% (over)        |

**Key Observations:**
- Height over-predicted by ~50%
- DBH under-predicted by 20-30%
- Volume errors are lower than other species, suggesting better overall calibration

## Root Cause Analysis

### 1. Site Index Base Age Mismatch

The primary issue appears to be a **site index base age mismatch**:
- Published hardwood yield tables use **base age 50**
- FVS-SN may internally convert or assume different base ages
- Height at index age should equal site index by definition

### 2. Height-Diameter Relationship

The FVS-SN Curtis-Arney height-diameter model may not be well-calibrated for hardwoods:
- Oaks are predicted too tall for their diameter
- Yellow-poplar shows similar but less severe pattern

### 3. Growth Rate Differences

Stand dynamics differ between published tables and simulated stands:
- Published tables represent **natural regeneration** with variable density
- Simulated stands start with **uniform planted** conditions
- Different competition environments affect individual tree growth

## Recommendations

### Short-term (Documentation)

1. **Document limitations** - Note that hardwood predictions have higher uncertainty
2. **Use relative comparisons** - Focus on management scenario comparisons rather than absolute predictions
3. **Apply adjustment factors** - Consider species-specific correction factors

### Medium-term (Calibration)

1. **Site index conversion** - Implement base age conversion when loading hardwood configs
2. **Height-diameter recalibration** - Update Curtis-Arney coefficients from recent FIA data
3. **Ecounit effects** - Verify hardwood ecounit modifiers are correctly applied

### Long-term (Validation)

1. **FIA comparison** - Compare against FIA remeasurement data for same species/regions
2. **Official FVS comparison** - Run identical scenarios through official FVS-SN
3. **Regional calibration** - Develop region-specific coefficients

## Acceptance Criteria for Hardwoods

Given the current calibration state, suggested tolerances for hardwood validation:

| Metric | Strict | Moderate | Lenient |
|--------|--------|----------|---------|
| Height | ±15%   | ±25%     | ±40%    |
| DBH    | ±20%   | ±30%     | ±40%    |
| Volume | ±30%   | ±50%     | ±100%   |

**Current status:** Most hardwood species fail even "lenient" height criteria but pass DBH criteria.

## Conclusion

The PyFVS hardwood growth models produce biologically reasonable growth patterns (trees grow, mortality occurs, correct directional relationships) but have systematic biases compared to published yield tables. The models are usable for:

- **Relative comparisons** between management scenarios
- **Trend analysis** over time
- **Species composition** in mixed stands

The models should NOT be used for:
- **Absolute yield predictions** where accuracy matters
- **Timber inventory projections** for hardwood stands
- **Financial analysis** dependent on accurate volume estimates

Further calibration work is needed to bring hardwood predictions in line with published data.

---

## Data Sources

1. Beck, D.E. and Della-Bianca, L. 1970. Yield of unthinned yellow-poplar. USDA FS Research Paper SE-58.
2. Schnur, G.L. 1937. Yield, Stand, and Volume Tables for Even-Aged Upland Oak Forests. USDA Tech. Bull. 560.
3. Burns, R.M. and Honkala, B.H. 1990. Silvics of North America, Vol. 2. USDA FS Agriculture Handbook 654.
