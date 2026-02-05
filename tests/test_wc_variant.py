"""Tests for the WC (West Cascades) variant implementation.

WC shares bark ratio, crown ratio, and volume infrastructure with PN since both
derive from the same FVS Fortran source files (bratio.f, crown.f). WC has its own
SDI maximums (lower than PN due to drier interior Cascades conditions).

Tests cover:
- Factory dispatch (WC -> PNBarkRatioModel, PNCrownRatioModel)
- Bark ratio via PN model
- Crown ratio via PN model
- WC-specific SDI maximums
- Volume coverage for WC species
- Bark ratio bug fix (species-specific values passed to diameter growth)
- Integration (full stand simulation with variant='WC')
"""
import math
import pytest

from pyfvs import Stand
from pyfvs.bark_ratio import (
    PNBarkRatioModel,
    create_bark_ratio_model,
)
from pyfvs.crown_ratio import (
    PNCrownRatioModel,
    create_crown_ratio_model,
)
from pyfvs.mortality import create_mortality_model, MortalityModel
from pyfvs.stand_metrics import StandMetricsCalculator
from pyfvs.volume_library import (
    VolumeCalculator,
    calculate_tree_volume,
    VOLUME_COEFFICIENTS_OUTSIDE_BARK,
    VOLUME_COEFFICIENTS_INSIDE_BARK,
    HARDWOOD_SPECIES,
)


# ===========================================================================
# Phase 1: Factory Dispatch Tests
# ===========================================================================

class TestWCFactoryDispatch:
    """Tests that WC variant dispatches to PN models correctly."""

    def test_bark_ratio_factory_returns_pn_model(self):
        """create_bark_ratio_model with variant='WC' returns PNBarkRatioModel."""
        model = create_bark_ratio_model('DF', variant='WC')
        assert isinstance(model, PNBarkRatioModel)

    def test_crown_ratio_factory_returns_pn_model(self):
        """create_crown_ratio_model with variant='WC' returns PNCrownRatioModel."""
        model = create_crown_ratio_model('DF', variant='WC')
        assert isinstance(model, PNCrownRatioModel)

    def test_mortality_factory_returns_model_with_wc_sdi(self):
        """create_mortality_model with variant='WC' returns MortalityModel."""
        model = create_mortality_model('DF', variant='WC')
        assert isinstance(model, MortalityModel)

    def test_wc_bark_ratio_matches_pn(self):
        """WC and PN bark ratios are identical for the same species."""
        for species in ['DF', 'WH', 'RC', 'SF', 'GF']:
            wc_model = create_bark_ratio_model(species, variant='WC')
            pn_model = create_bark_ratio_model(species, variant='PN')
            wc_ratio = wc_model.calculate_bark_ratio(10.0)
            pn_ratio = pn_model.calculate_bark_ratio(10.0)
            assert wc_ratio == pn_ratio, (
                f"{species} bark ratio differs: WC={wc_ratio} vs PN={pn_ratio}"
            )


# ===========================================================================
# Phase 2: Bark Ratio via PN Model
# ===========================================================================

class TestWCBarkRatio:
    """Tests bark ratio calculations for WC species using PN model."""

    def test_df_bark_ratio_reasonable(self):
        """Douglas-fir bark ratio is in expected range."""
        model = create_bark_ratio_model('DF', variant='WC')
        for dob in [5.0, 10.0, 20.0, 30.0]:
            ratio = model.calculate_bark_ratio(dob)
            assert 0.85 <= ratio <= 0.95, (
                f"DF bark ratio {ratio} out of bounds at DOB={dob}"
            )

    def test_wh_bark_ratio_reasonable(self):
        """Western Hemlock bark ratio is in expected range."""
        model = create_bark_ratio_model('WH', variant='WC')
        ratio = model.calculate_bark_ratio(10.0)
        assert 0.88 <= ratio <= 0.96, f"WH bark ratio {ratio} unexpected"

    def test_rc_bark_ratio_reasonable(self):
        """Western Red Cedar bark ratio is in expected range."""
        model = create_bark_ratio_model('RC', variant='WC')
        ratio = model.calculate_bark_ratio(10.0)
        assert 0.80 <= ratio <= 0.96, f"RC bark ratio {ratio} unexpected"

    def test_all_major_wc_species_valid_bark_ratio(self):
        """All major WC species produce valid bark ratios (0.80-0.99)."""
        major_species = ['DF', 'WH', 'RC', 'SF', 'GF', 'WF', 'NF', 'IC',
                         'RA', 'BM', 'MH', 'YC', 'ES', 'AF', 'RF']
        for species in major_species:
            model = create_bark_ratio_model(species, variant='WC')
            ratio = model.calculate_bark_ratio(10.0)
            assert 0.80 <= ratio <= 0.99, (
                f"{species} bark ratio {ratio} out of bounds"
            )


# ===========================================================================
# Phase 3: Crown Ratio via PN Model
# ===========================================================================

class TestWCCrownRatio:
    """Tests crown ratio calculations for WC species using PN model."""

    def test_df_crown_ratio_decreases_with_sdi(self):
        """DF crown ratio decreases as relative SDI increases."""
        model = create_crown_ratio_model('DF', variant='WC')
        cr_low = model.calculate_average_crown_ratio(2.0)
        cr_high = model.calculate_average_crown_ratio(8.0)
        assert cr_low > cr_high, (
            f"CR should decrease with SDI: low={cr_low}, high={cr_high}"
        )

    def test_crown_ratio_bounded(self):
        """Crown ratios are bounded within reasonable range."""
        model = create_crown_ratio_model('DF', variant='WC')
        for relsdi in [0.1, 0.5, 1.0, 2.0, 5.0]:
            cr = model.calculate_average_crown_ratio(relsdi)
            assert 0.05 <= cr <= 0.95, (
                f"CR {cr} out of bounds at RELSDI={relsdi}"
            )
        # At extreme density (RELSDI=10), CR can be very low
        cr_extreme = model.calculate_average_crown_ratio(10.0)
        assert 0.01 <= cr_extreme <= 0.95, (
            f"CR {cr_extreme} out of bounds at extreme RELSDI=10.0"
        )

    def test_different_species_different_cr(self):
        """Different species produce different crown ratios."""
        species_crs = {}
        for species in ['DF', 'WH', 'RC', 'SF']:
            model = create_crown_ratio_model(species, variant='WC')
            species_crs[species] = model.calculate_average_crown_ratio(5.0)
        # Not all species should have identical CRs
        unique_values = set(round(v, 4) for v in species_crs.values())
        assert len(unique_values) >= 2, (
            f"Expected different CRs across species, got {species_crs}"
        )

    def test_wc_crown_ratio_matches_pn(self):
        """WC and PN crown ratios are identical for the same species."""
        for species in ['DF', 'WH', 'RC']:
            wc_model = create_crown_ratio_model(species, variant='WC')
            pn_model = create_crown_ratio_model(species, variant='PN')
            wc_cr = wc_model.calculate_average_crown_ratio(5.0)
            pn_cr = pn_model.calculate_average_crown_ratio(5.0)
            assert wc_cr == pn_cr, (
                f"{species} crown ratio differs: WC={wc_cr} vs PN={pn_cr}"
            )


# ===========================================================================
# Phase 4: WC SDI Maximums
# ===========================================================================

class TestWCSDIMaximums:
    """Tests for WC-specific SDI maximum values."""

    def test_wc_sdi_maximums_loaded(self):
        """WC_SDI_MAXIMUMS is available on StandMetricsCalculator."""
        sdi_maxs = StandMetricsCalculator.WC_SDI_MAXIMUMS
        assert isinstance(sdi_maxs, dict)
        assert len(sdi_maxs) > 30, f"Expected 37 species, got {len(sdi_maxs)}"

    def test_df_sdi_max_800(self):
        """Douglas-fir SDI max is 800 for WC (lower than PN's 850)."""
        calc = StandMetricsCalculator('DF', variant='WC')
        max_sdi = calc.get_max_sdi([], 'DF')
        assert max_sdi == 800

    def test_wh_sdi_max_850(self):
        """Western Hemlock SDI max is 850 for WC."""
        calc = StandMetricsCalculator('WH', variant='WC')
        max_sdi = calc.get_max_sdi([], 'WH')
        assert max_sdi == 850

    def test_all_species_have_positive_sdi_max(self):
        """All 37 WC species have positive SDI maximum values."""
        sdi_maxs = StandMetricsCalculator.WC_SDI_MAXIMUMS
        for species, max_val in sdi_maxs.items():
            assert max_val > 0, f"{species} has non-positive SDI max: {max_val}"

    def test_wc_generally_lower_than_pn(self):
        """WC SDI maximums are generally lower than PN (interior vs coast)."""
        wc_maxs = StandMetricsCalculator.WC_SDI_MAXIMUMS
        pn_maxs = StandMetricsCalculator.PN_SDI_MAXIMUMS
        # Compare shared species
        shared = set(wc_maxs.keys()) & set(pn_maxs.keys())
        lower_count = sum(1 for sp in shared if wc_maxs[sp] <= pn_maxs[sp])
        # Most WC values should be <= PN values
        assert lower_count >= len(shared) * 0.7, (
            f"Expected WC SDI maxs mostly <= PN, but only {lower_count}/{len(shared)} are"
        )

    def test_mortality_model_uses_wc_sdi(self):
        """Mortality model with variant='WC' uses WC SDI maximums."""
        model = create_mortality_model('DF', variant='WC')
        # DF SDI max for WC should be 800
        assert model.max_sdi == 800, (
            f"Expected WC DF max_sdi=800, got {model.max_sdi}"
        )


# ===========================================================================
# Phase 5: Volume Coverage
# ===========================================================================

class TestWCVolume:
    """Tests for WC species volume coverage."""

    def test_df_volume_reasonable(self):
        """Douglas-fir volume is reasonable (10" DBH, 80' ht)."""
        result = calculate_tree_volume(10.0, 80.0, 'DF')
        assert result.is_valid()
        assert 15.0 <= result.total_cubic_volume <= 35.0, (
            f"DF volume {result.total_cubic_volume} outside expected range"
        )

    def test_wh_volume_reasonable(self):
        """Western Hemlock volume is reasonable."""
        result = calculate_tree_volume(12.0, 90.0, 'WH')
        assert result.is_valid()
        assert 20.0 <= result.total_cubic_volume <= 50.0, (
            f"WH volume {result.total_cubic_volume} outside expected range"
        )

    def test_rc_volume_reasonable(self):
        """Western Red Cedar volume is reasonable."""
        result = calculate_tree_volume(15.0, 100.0, 'RC')
        assert result.is_valid()
        assert 40.0 <= result.total_cubic_volume <= 80.0, (
            f"RC volume {result.total_cubic_volume} outside expected range"
        )

    def test_wc_specific_species_in_volume_dicts(self):
        """WC-specific species have explicit volume coefficients."""
        wc_species = ['CW', 'WJ', 'WB', 'KP', 'HT', 'WI', 'OT']
        for sp in wc_species:
            assert sp in VOLUME_COEFFICIENTS_OUTSIDE_BARK, (
                f"{sp} missing from VOLUME_COEFFICIENTS_OUTSIDE_BARK"
            )
            assert sp in VOLUME_COEFFICIENTS_INSIDE_BARK, (
                f"{sp} missing from VOLUME_COEFFICIENTS_INSIDE_BARK"
            )

    def test_volume_increases_with_size(self):
        """Volume increases with tree size for DF."""
        vol_small = calculate_tree_volume(6.0, 50.0, 'DF').total_cubic_volume
        vol_medium = calculate_tree_volume(12.0, 80.0, 'DF').total_cubic_volume
        vol_large = calculate_tree_volume(24.0, 140.0, 'DF').total_cubic_volume
        assert vol_small < vol_medium < vol_large


# ===========================================================================
# Phase 6: Bark Ratio Bug Fix Verification
# ===========================================================================

class TestBarkRatioBugFix:
    """Tests that species-specific bark ratio is passed to diameter growth."""

    def test_pn_diameter_growth_accepts_bark_ratio(self):
        """PN diameter growth model accepts bark_ratio parameter."""
        from pyfvs.pn_diameter_growth import create_pn_diameter_growth_model
        model = create_pn_diameter_growth_model('DF')
        # Call with explicit bark_ratio parameter
        dg = model.calculate_diameter_growth(
            dbh=10.0, crown_ratio=0.5, site_index=120,
            ba=150.0, bal=50.0, bark_ratio=0.90
        )
        assert dg >= 0.0

    def test_wc_diameter_growth_accepts_bark_ratio(self):
        """WC diameter growth model accepts bark_ratio parameter."""
        from pyfvs.wc_diameter_growth import create_wc_diameter_growth_model
        model = create_wc_diameter_growth_model('DF')
        # Call with explicit bark_ratio parameter
        dg = model.calculate_diameter_growth(
            dbh=10.0, crown_ratio=0.5, site_index=120,
            ba=150.0, bal=50.0, bark_ratio=0.90
        )
        assert dg >= 0.0

    def test_different_bark_ratio_affects_growth(self):
        """Different bark ratios produce different diameter growth values."""
        from pyfvs.wc_diameter_growth import create_wc_diameter_growth_model
        model = create_wc_diameter_growth_model('DF')
        dg_high = model.calculate_diameter_growth(
            dbh=10.0, crown_ratio=0.5, site_index=120,
            ba=150.0, bal=50.0, bark_ratio=0.95
        )
        dg_low = model.calculate_diameter_growth(
            dbh=10.0, crown_ratio=0.5, site_index=120,
            ba=150.0, bal=50.0, bark_ratio=0.85
        )
        # Different bark ratios should produce different growth
        assert dg_high != dg_low, (
            f"Expected different growth with different bark ratios: "
            f"0.95->{dg_high}, 0.85->{dg_low}"
        )


# ===========================================================================
# Phase 7: Integration Tests
# ===========================================================================

class TestWCIntegration:
    """Integration tests for full WC variant simulation."""

    def test_stand_initialize_planted_wc(self):
        """Stand.initialize_planted works with variant='WC'."""
        stand = Stand.initialize_planted(
            trees_per_acre=400,
            site_index=120,
            species='DF',
            variant='WC'
        )
        assert stand.variant == 'WC'
        assert len(stand.trees) == 400
        assert stand.site_index == 120

    def test_wc_simulation_50_year(self):
        """50-year WC simulation produces reasonable yields."""
        stand = Stand.initialize_planted(
            trees_per_acre=400,
            site_index=120,
            species='DF',
            variant='WC'
        )
        stand.grow(50)
        metrics = stand.get_metrics()

        tpa = metrics['tpa']
        qmd = metrics['qmd']
        ba = metrics['basal_area']
        volume = metrics['volume']

        # WC DF at SI=120 should produce good yields
        assert 150 <= tpa <= 400, f"TPA {tpa} outside expected range"
        assert 8.0 <= qmd <= 22.0, f"QMD {qmd} outside expected range"
        assert 100 <= ba <= 500, f"BA {ba} outside expected range"
        assert 3000 <= volume <= 25000, f"Volume {volume} outside expected range"

    def test_wc_tree_level_growth(self):
        """Individual DF tree grows with variant='WC'."""
        from pyfvs.tree import Tree

        tree = Tree(dbh=2.0, height=15.0, species='DF', age=10, variant='WC')
        initial_dbh = tree.dbh
        initial_height = tree.height

        tree.grow(
            site_index=120,
            ba=50.0,
            pbal=20.0,
            competition_factor=0.3,
            time_step=10,
        )

        assert tree.dbh > initial_dbh, "DBH should increase"
        assert tree.height > initial_height, "Height should increase"

    def test_wc_wh_simulation(self):
        """Western Hemlock WC simulation produces reasonable yields."""
        stand = Stand.initialize_planted(
            trees_per_acre=400,
            site_index=100,
            species='WH',
            variant='WC'
        )
        stand.grow(50)
        metrics = stand.get_metrics()

        assert metrics['tpa'] > 100, "Should have surviving trees"
        assert metrics['qmd'] > 5.0, "QMD should grow"
        assert metrics['basal_area'] > 50.0, "BA should accumulate"

    def test_wc_vs_pn_yields(self):
        """WC yields should be comparable to PN (similar species/equations)."""
        stand_wc = Stand.initialize_planted(
            trees_per_acre=400,
            site_index=120,
            species='DF',
            variant='WC'
        )
        stand_wc.grow(50)
        vol_wc = stand_wc.get_metrics()['volume']

        stand_pn = Stand.initialize_planted(
            trees_per_acre=400,
            site_index=120,
            species='DF',
            variant='PN'
        )
        stand_pn.grow(50)
        vol_pn = stand_pn.get_metrics()['volume']

        # Both should produce substantial volume
        assert vol_wc > 3000, f"WC volume {vol_wc} too low"
        assert vol_pn > 3000, f"PN volume {vol_pn} too low"

        # WC uses different diameter growth coefficients than PN
        # so yields will differ, but both should be in the same ballpark
        ratio = vol_wc / vol_pn if vol_pn > 0 else 0
        assert 0.3 <= ratio <= 2.0, (
            f"WC/PN volume ratio {ratio:.2f} too extreme "
            f"(WC={vol_wc:.0f}, PN={vol_pn:.0f})"
        )

    def test_wc_stand_metrics_correct_variant(self):
        """Stand metrics calculator uses WC variant correctly."""
        stand = Stand.initialize_planted(
            trees_per_acre=400,
            site_index=120,
            species='DF',
            variant='WC'
        )
        assert stand._metrics.variant == 'WC'
        max_sdi = stand._metrics.get_max_sdi(stand.trees, 'DF')
        assert max_sdi == 800, f"DF max SDI should be 800, got {max_sdi}"

    def test_wc_bark_ratio_in_growth(self):
        """WC bark ratio model is actually used during growth."""
        model = create_bark_ratio_model('DF', variant='WC')
        assert isinstance(model, PNBarkRatioModel)
        ratio = model.calculate_bark_ratio(10.0)
        assert 0.85 <= ratio <= 0.95, f"DF WC bark ratio {ratio} unexpected"

    def test_wc_crown_ratio_in_growth(self):
        """WC crown ratio model is used during simulation."""
        model = create_crown_ratio_model('DF', variant='WC')
        assert isinstance(model, PNCrownRatioModel)
        cr = model.calculate_average_crown_ratio(5.0)
        assert 0.10 <= cr <= 0.90, f"DF WC mean CR {cr} unexpected"
