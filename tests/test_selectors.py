import pytest

from ilmn.pelops import entities, repositories, selectors
from ilmn.pelops.callers import region_callers

# fixture `an_unnamed_genomic_region_set` is definded in conftest


class TestSelectorFactory:
    @pytest.fixture
    def factory(self):
        regions_repository = repositories.BuiltinRegionRepository()
        factory = selectors.SelectorFactory(regions_repository)
        return factory

    def test_build(self, factory):
        observed = factory.build()
        assert len(observed) == 3
        for selector in observed:
            assert isinstance(selector, selectors.RegionSelector)

    def test_build_with_blacklist(self, factory):
        features = frozenset([selectors.SelectorFeature.WITH_BLACKLIST])
        observed = factory.build(features=features)
        assert len(observed) == 4
        for selector in observed:
            assert isinstance(selector, selectors.RegionSelector)


selector_data = [
    pytest.param(
        (entities.RegionsName.IGH, ("chr1", 100, 200), False), id="different_chromosome"
    ),
    pytest.param(
        (entities.RegionsName.IGH, ("chr14", 105587000, 105588000), True),
        id="overlaps_only_target",
    ),
    pytest.param(
        (entities.RegionsName.CoreDUX4, ("chr4", 190023666, 190066934), False),
        id="not_overlapping",
    ),
    pytest.param(
        (entities.RegionsName.CoreDUX4, ("chr4", 190172700, 190172800), True),
        id="overlaps_one_target",
    ),
]


@pytest.fixture(params=selector_data)
def test_case(request):
    region_name, region, expected = request.param
    return region_name, region, expected


@pytest.fixture
def provided(test_case):
    region_name, region, expected = test_case
    result = entities.GenomicRegion(*region)
    return result


@pytest.fixture
def expected(test_case):
    region_name, region, expected = test_case
    return expected


@pytest.fixture
def region_name(test_case):
    region_name, region, expected = test_case
    return region_name


class TestNamedRegionSelector:
    def test_call(self, provided, expected, region_name):
        region_repo = repositories.BuiltinRegionRepository()
        selector = selectors.NamedRegionSelector(region_name, region_repo)
        assert selector(provided) == expected


class TestRegionDeSelector:
    def test_call(self, provided, expected, region_name):
        region_repo = repositories.BuiltinRegionRepository()
        selector = selectors.NamedRegionSelector(region_name, region_repo)
        deselector = selectors.RegionDeSelector(selector)
        assert deselector(provided) == (not expected)


class TestRearrangementSelector:
    test_cases = [
        pytest.param((3, 1, True), id="above_threshold"),
        pytest.param((3, 5, False), id="below_threshold"),
        pytest.param((3, 3.0, True), id="same_as_threshold"),
    ]

    @pytest.fixture(params=test_cases)
    def test_case(self, request):
        srpb, threshold, expected = request.param
        return srpb, threshold, expected

    @pytest.fixture
    def rearrangement(self, an_unnamed_genomic_region_set, test_case):
        srpb, threshold, expected = test_case
        counts = entities.ClassifiedSegmentCount(paired=2, split=5, spanning=7)
        region_pairs = entities.CompoundRegionPair(
            an_unnamed_genomic_region_set, an_unnamed_genomic_region_set
        )
        result = entities.Rearrangement(region_pairs, counts, srpb, segments=set([]))
        return result

    def test_call(self, rearrangement, test_case):
        _, threshold, expected = test_case
        selector = selectors.SRPBRearrangementSelector(threshold)
        assert selector(rearrangement) == expected
