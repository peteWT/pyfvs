"""
Shared pytest fixtures for FVS-Python tests.

This module provides commonly used fixtures for testing tree and stand
functionality, reducing code duplication across test files.
"""
import pytest
from pathlib import Path
from pyfvs.tree import Tree
from pyfvs.stand import Stand


# =============================================================================
# Output Directory Setup
# =============================================================================

@pytest.fixture(scope="session")
def test_output_dir():
    """Create and return the test output directory.

    Session-scoped to create directories only once per test run.
    """
    output_dir = Path(__file__).parent.parent / 'test_output'
    for subdir in ['tree_tests', 'stand_tests', 'tree_comprehensive_tests']:
        (output_dir / subdir).mkdir(parents=True, exist_ok=True)
    return output_dir


# =============================================================================
# Tree Fixtures - Individual Trees at Various Development Stages
# =============================================================================

@pytest.fixture
def small_tree():
    """Create a small tree for testing the small-tree growth model.

    Returns a tree with:
    - DBH: 1.0 inch (well below Xmin threshold of 1.0")
    - Height: 6.0 feet
    - Age: 2 years
    - Species: LP (Loblolly Pine, default)

    This tree uses the small-tree model (height-driven growth).
    """
    return Tree(dbh=1.0, height=6.0, age=2)


@pytest.fixture
def large_tree():
    """Create a large tree for testing the large-tree growth model.

    Returns a tree with:
    - DBH: 6.0 inches (above Xmax threshold of 3.0")
    - Height: 40.0 feet
    - Age: 15 years
    - Species: LP (Loblolly Pine, default)

    This tree uses the large-tree model (diameter-driven growth).
    """
    return Tree(dbh=6.0, height=40.0, age=15)


@pytest.fixture
def transition_tree():
    """Create a tree in the transition zone between growth models.

    Returns a tree with:
    - DBH: 2.5 inches (between Xmin=1.0" and Xmax=3.0")
    - Height: 20.0 feet
    - Age: 8 years
    - Species: LP (Loblolly Pine, default)

    This tree uses a weighted blend of small-tree and large-tree models.
    """
    return Tree(dbh=2.5, height=20.0, age=8)


@pytest.fixture
def seedling_tree():
    """Create a very small seedling tree.

    Returns a tree with:
    - DBH: 0.5 inches
    - Height: 2.0 feet
    - Age: 1 year
    - Species: LP (Loblolly Pine, default)
    """
    return Tree(dbh=0.5, height=2.0, age=1)


@pytest.fixture
def pole_tree():
    """Create a pole-sized tree.

    Returns a tree with:
    - DBH: 6.0 inches (pole timber: 5-9" DBH)
    - Height: 40.0 feet
    - Age: 15 years
    - Species: LP (Loblolly Pine, default)
    """
    return Tree(dbh=6.0, height=40.0, age=15)


@pytest.fixture
def sawtimber_tree():
    """Create a sawtimber-sized tree.

    Returns a tree with:
    - DBH: 12.0 inches (sawtimber: >= 9" DBH for softwoods)
    - Height: 70.0 feet
    - Age: 30 years
    - Species: LP (Loblolly Pine, default)
    """
    return Tree(dbh=12.0, height=70.0, age=30)


@pytest.fixture
def mature_tree():
    """Create a mature tree at economic rotation age.

    Returns a tree with:
    - DBH: 14.0 inches
    - Height: 85.0 feet
    - Age: 40 years
    - Species: LP (Loblolly Pine, default)
    """
    return Tree(dbh=14.0, height=85.0, age=40)


# =============================================================================
# Sample Trees Lists - For Stand-Level Testing
# =============================================================================

@pytest.fixture
def sample_trees():
    """Create a sample list of trees for stand-level testing.

    Returns a list of 10 trees with:
    - DBH: 4.0 to 8.5 inches (0.5" increments)
    - Height: 30.0 to 57.0 feet (3.0' increments)
    - Age: 15 years
    - Species: LP (Loblolly Pine)

    Useful for testing metrics calculations, competition, etc.
    """
    trees = []
    for i in range(10):
        dbh = 4.0 + i * 0.5  # DBH from 4.0 to 8.5 inches
        height = 30.0 + i * 3.0  # Heights from 30 to 57 feet
        tree = Tree(dbh=dbh, height=height, species='LP', age=15)
        trees.append(tree)
    return trees


@pytest.fixture
def sample_trees_50():
    """Create a larger sample list of 50 trees for testing.

    Returns a list of 50 trees with:
    - DBH: 4.0 to 13.8 inches (0.2" increments)
    - Height: 30.0 to 69.2 feet (0.8' increments)
    - Age: 15 years
    - Species: LP (Loblolly Pine)
    """
    trees = []
    for i in range(50):
        dbh = 4.0 + i * 0.2  # DBH from 4.0 to 13.8 inches
        height = 30.0 + i * 0.8  # Heights from 30 to 69.2 feet
        trees.append(Tree(dbh=dbh, height=height, species='LP', age=15))
    return trees


@pytest.fixture
def uniform_stand_trees():
    """Create a uniform stand with identical trees.

    Returns a list of 25 identical trees with:
    - DBH: 8.0 inches
    - Height: 50.0 feet
    - Age: 20 years
    - Species: LP (Loblolly Pine)

    Useful for testing calculations that should be uniform.
    """
    return [Tree(dbh=8.0, height=50.0, species='LP', age=20) for _ in range(25)]


@pytest.fixture
def mixed_species_trees():
    """Create trees with mixed species.

    Returns a list of 5 trees:
    - 2 LP (Loblolly Pine)
    - 2 SP (Shortleaf Pine)
    - 1 SA (Slash Pine)
    """
    trees = []
    species_list = ['LP', 'LP', 'SP', 'SP', 'SA']
    for i, species in enumerate(species_list):
        dbh = 5.0 + i * 1.0
        height = 35.0 + i * 5.0
        tree = Tree(dbh=dbh, height=height, species=species, age=20)
        trees.append(tree)
    return trees


# =============================================================================
# Mortality Testing Fixtures - Stands at Various Density Levels
# =============================================================================

@pytest.fixture
def sparse_stand_trees():
    """Create a sparse stand (below 55% SDI threshold).

    Returns a list of 50 trees representing a low-density stand
    where only background mortality applies (below SDI threshold).
    """
    trees = []
    for i in range(50):  # Low TPA
        dbh = 4.0 + i * 0.2
        height = 30.0 + i * 1.0
        trees.append(Tree(dbh=dbh, height=height, species='LP', age=15))
    return trees


@pytest.fixture
def dense_stand_trees():
    """Create a dense stand (above 55% SDI threshold).

    Returns a list of 500 trees representing a high-density stand
    where density-dependent mortality kicks in.
    """
    trees = []
    for i in range(500):  # High TPA
        dbh = 6.0 + (i % 20) * 0.3
        height = 40.0 + (i % 20) * 1.5
        trees.append(Tree(dbh=dbh, height=height, species='LP', age=20))
    return trees


@pytest.fixture
def very_dense_stand_trees():
    """Create a very dense stand (above 85% SDI threshold).

    Returns a list of 800 trees representing a very high-density stand
    where maximum mortality occurs.
    """
    trees = []
    for i in range(800):  # Very high TPA
        dbh = 8.0 + (i % 15) * 0.2
        height = 50.0 + (i % 15) * 1.0
        trees.append(Tree(dbh=dbh, height=height, species='LP', age=25))
    return trees


# =============================================================================
# Stand Fixtures - Complete Stand Objects
# =============================================================================

# Standard TPA values for testing
STANDARD_TPA = 500  # Trees per acre for a typical plantation
LOW_TPA = 300       # Low density plantation
HIGH_TPA = 700      # High density plantation


@pytest.fixture
def young_stand():
    """Create a young planted stand (age 0).

    Returns a Stand with:
    - Trees per acre: 500 (standard planting density)
    - Site index: 70 (default)
    - Age: 0 (newly planted)
    - Species: LP (Loblolly Pine, default)
    """
    return Stand.initialize_planted(trees_per_acre=STANDARD_TPA)


@pytest.fixture
def mature_stand():
    """Create a mature stand by growing for 25 years.

    Returns a Stand with:
    - Initial trees per acre: 500
    - Site index: 70 (default)
    - Final age: 25 years
    - Species: LP (Loblolly Pine, default)

    Note: This fixture grows the stand, so it takes longer to create.
    """
    stand = Stand.initialize_planted(trees_per_acre=STANDARD_TPA)
    stand.grow(years=25)
    return stand


@pytest.fixture
def low_density_stand():
    """Create a low density planted stand.

    Returns a Stand with:
    - Trees per acre: 300
    - Site index: 70 (default)
    - Age: 0 (newly planted)
    """
    return Stand.initialize_planted(trees_per_acre=LOW_TPA)


@pytest.fixture
def high_density_stand():
    """Create a high density planted stand.

    Returns a Stand with:
    - Trees per acre: 700
    - Site index: 70 (default)
    - Age: 0 (newly planted)
    """
    return Stand.initialize_planted(trees_per_acre=HIGH_TPA)


@pytest.fixture
def high_site_stand():
    """Create a stand on a high-quality site.

    Returns a Stand with:
    - Trees per acre: 500
    - Site index: 80 (high)
    - Age: 0 (newly planted)
    """
    return Stand.initialize_planted(trees_per_acre=STANDARD_TPA, site_index=80)


@pytest.fixture
def low_site_stand():
    """Create a stand on a low-quality site.

    Returns a Stand with:
    - Trees per acre: 500
    - Site index: 55 (low)
    - Age: 0 (newly planted)
    """
    return Stand.initialize_planted(trees_per_acre=STANDARD_TPA, site_index=55)


@pytest.fixture
def mountain_ecounit_stand():
    """Create a stand in the M231 (Mountain) ecological unit.

    Returns a Stand with:
    - Trees per acre: 500
    - Site index: 70 (default)
    - Ecounit: M231 (highest growth rates for LP)
    - Age: 0 (newly planted)
    """
    return Stand.initialize_planted(
        trees_per_acre=STANDARD_TPA,
        site_index=70,
        species='LP',
        ecounit='M231'
    )


@pytest.fixture
def empty_stand():
    """Create an empty stand with no trees.

    Returns a Stand with:
    - No trees
    - Site index: 70 (default)

    Useful for testing edge cases.
    """
    return Stand(trees=[], site_index=70)


@pytest.fixture
def small_stand():
    """Create a small stand with fewer than 40 trees.

    Returns a Stand with:
    - Trees per acre: 30
    - Site index: 70 (default)
    - Age: 0 (newly planted)

    Useful for testing top height calculation with n < 40 trees.
    """
    return Stand.initialize_planted(trees_per_acre=30, site_index=70)


# =============================================================================
# Growth Testing Fixtures - Pre-configured growth parameters
# =============================================================================

@pytest.fixture
def standard_growth_params():
    """Return standard growth parameters for testing.

    Returns a dictionary with typical growth parameters.
    """
    return {
        'site_index': 70,
        'competition_factor': 0.2,
        'ba': 100,
        'pbal': 30,
        'slope': 0.05,
        'aspect': 0,
        'time_step': 5
    }


@pytest.fixture
def high_competition_params():
    """Return high competition growth parameters.

    Returns a dictionary with parameters representing high competition.
    """
    return {
        'site_index': 70,
        'competition_factor': 0.9,
        'rank': 0.2,
        'ba': 160,
        'pbal': 90,
        'slope': 0.05,
        'aspect': 0,
        'time_step': 5
    }


@pytest.fixture
def low_competition_params():
    """Return low competition growth parameters.

    Returns a dictionary with parameters representing low competition.
    """
    return {
        'site_index': 70,
        'competition_factor': 0.0,
        'rank': 0.8,
        'ba': 80,
        'pbal': 10,
        'slope': 0.05,
        'aspect': 0,
        'time_step': 5
    }


# =============================================================================
# Pytest Configuration
# =============================================================================

def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
