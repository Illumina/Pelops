import pathlib
from unittest import mock

import pytest

from ilmn.pelops import entities, notifications, repositories, selectors
from ilmn.pelops.callers import (
    caller_factories,
    rearrangement_callers,
    region_pair_callers,
    sorters,
)
from ilmn.pelops.infrastructure import blacklist_region_repository, pysam_repositories
from tests import stubs


@pytest.fixture
def region_repository():
    return stubs.SmallRegionRepository()


@pytest.fixture
def notification_service(notification_factory):
    result = notification_factory.build()
    return result


@pytest.fixture
def region_repository_factory(region_repository):
    factory = mock.Mock(spec=repositories.RegionRepositoryFactory, instance=True)
    factory.build = mock.Mock(return_value=region_repository)
    return factory


class TestNamedRegionPairCaller:
    def test_get_compound_region_pairs(self, region_repository):
        igh = entities.RegionsName.IGH
        dux4 = entities.RegionsName.CoreDUX4
        caller = region_pair_callers.NamedRegionPairCaller(region_repository)
        caller.name_pairs = [(dux4, igh)]
        expected = [
            entities.CompoundRegionPair(
                region_repository.get(dux4), region_repository.get(igh)
            )
        ]
        observed = list(caller.get_compound_region_pairs())
        assert observed == expected


class TestNotifyRegionPairCaller:
    @pytest.fixture
    def helper_caller(self, region_repository):
        result = region_pair_callers.NamedRegionPairCaller(region_repository)
        return result

    def test_get_compound_region_pairs(self, notification_service, helper_caller):
        caller = notifications.NotifyRegionPairCaller(
            helper_caller, notification_service
        )
        region_pairs = list(caller.get_compound_region_pairs())
        observed = notification_service.get_notifications()
        expected = [
            "Retrieving region coordinates for CoreDUX4/IGH.",
            "Retrieving region coordinates for ExtendedDUX4/IGH.",
        ]
        assert observed == expected


@pytest.fixture
def segment_repo(region_repository):
    """A stub for a segment repository"""
    # 4 fake segments from "dux4" where the mate read maps to two
    # contiguous regions on chr7
    mate_locations = [
        ("chr7", 1011),
        ("chr7", 1012),
        ("chr7", 2021),
        ("chr7", 2022),
    ]
    fake_paired_reads = [
        entities.PlacedSegment(
            "foo",
            region_repository.get(entities.RegionsName.CoreDUX4).regions,
            entities.ReadOrder.ONE,
            stubs.MockContent(*items),
        )
        for items in mate_locations
    ]
    segment_repo = stubs.StubPlacedSegmentRepository()
    for read in fake_paired_reads:
        segment_repo.add_read(read)
    return segment_repo


@pytest.fixture
def region_pair_caller(
    segment_repo, segment_repo_factory, region_repository_factory, notification_factory
):
    # test the `CompoundRegionPairCaller` built by the factory
    segment_repo_factory.build = mock.Mock(return_value=segment_repo)

    candidate_region_caller_factory = caller_factories.CandidateRegionCallerFactory(
        region_repository_factory,
        segment_repo_factory,
    )

    factory = caller_factories.RegionPairCallerFactory(
        candidate_region_caller_factory, region_repository_factory, notification_factory
    )
    result = factory.build(
        caller_factories.CallerType.BAIT, frozenset(), total_number_of_reads=None
    )
    return result


class TestBaitRegionPairCaller:
    def test_get_region_sets(self, region_pair_caller, region_repository):
        igh = entities.RegionsName.IGH
        dux4 = entities.RegionsName.CoreDUX4
        expected_candidate_region_set = entities.CompoundRegion(
            name=entities.RegionsName.UNNAMED,
            regions=frozenset([entities.GenomicRegion("chr7", 1001, 3000)]),
        )
        expected = [
            entities.CompoundRegionPair(
                region_repository.get(dux4),
                expected_candidate_region_set,
            ),
        ]
        observed = list(region_pair_caller.get_compound_region_pairs())
        for pair in expected:
            for regionset in pair:
                for region in regionset.regions:
                    print(region.chrom, region.start, region.end)
        assert observed == expected


class TestSelectableRearrangementCaller:
    @pytest.fixture
    def provided(self, an_unnamed_genomic_region_set):
        counts = entities.ClassifiedSegmentCount(1, 2, 3)
        region_pair = entities.CompoundRegionPair(
            an_unnamed_genomic_region_set, an_unnamed_genomic_region_set
        )
        result = iter(
            [entities.Rearrangement(region_pair, counts, 3, segments=set([]))]
        )
        return result

    @pytest.fixture
    def input_caller(self, provided):
        result = mock.Mock(
            spec=rearrangement_callers.RearrangementCaller, instance=True
        )
        result.get_rearrangements = mock.Mock(return_value=provided)
        return result

    def test_get_rearrangements(self, input_caller, provided):
        selector = selectors.SRPBRearrangementSelector(4)
        caller = rearrangement_callers.SelectableRearrangementCaller(
            input_caller, selector
        )
        observed = list(caller.get_rearrangements())
        assert observed == []
        selector = selectors.SRPBRearrangementSelector(3)
        caller = rearrangement_callers.SelectableRearrangementCaller(
            input_caller, selector
        )
        observed = list(caller.get_rearrangements())
        assert observed == list(provided)


class TestSortedRearrangementCaller:
    @pytest.fixture
    def provided(self):
        return [4, 5, 3]

    @pytest.fixture
    def expected(self):
        return [5, 4, 3]

    def _create_one_rearrangement(self, srpb, regions):
        counts = entities.ClassifiedSegmentCount(1, 2, 3)
        region_pair = entities.CompoundRegionPair(regions, regions)
        result = entities.Rearrangement(region_pair, counts, srpb, segments=set([]))
        return result

    @pytest.fixture
    def input_caller(self, provided, an_unnamed_genomic_region_set):
        rearrangements = [
            self._create_one_rearrangement(srpb, an_unnamed_genomic_region_set)
            for srpb in provided
        ]
        result = mock.Mock(spec=rearrangement_callers.RearrangementCaller)
        result.get_rearrangements = mock.Mock(return_value=rearrangements)
        return result

    def test_get_rearrangements(
        self, input_caller, expected, an_unnamed_genomic_region_set
    ):
        expected_rearrangements = [
            self._create_one_rearrangement(srpb, an_unnamed_genomic_region_set)
            for srpb in expected
        ]
        sorter = sorters.SrpbRearrangementSorter()
        caller = rearrangement_callers.SortedRearrangementCaller(input_caller, sorter)
        observed = list(caller.get_rearrangements())
        assert observed == expected_rearrangements


class TestRearrangementCallerFromFactories:
    """Test the RearrangementCallers built by the factories.
    This is a sort of "intergration test" for the RearrangementCallers, as it
    won't be easy to test every possible case (tested in various compontent)
    but still it provides some evidence that the instance built by the factory has
    at least some expected behaviour"""

    @pytest.fixture
    def factory(
        self,
        segment_repo_factory,
        segment_repo,
        notification_factory,
        candidate_region_caller_factory,
    ):
        repo_factory = pysam_repositories.PysamSegmentRepositoryFactory(
            pathlib.Path("not/used.bam"), notification_factory
        )
        region_repo_factory = blacklist_region_repository.FileRegionRepositoryFactory()

        region_caller_factory = caller_factories.RegionPairCallerFactory(
            candidate_region_caller_factory, region_repo_factory, notification_factory
        )
        segment_repo_factory.build = mock.Mock(return_value=segment_repo)
        read_caller_factory = caller_factories.ReadCallerFactory(
            segment_repo_factory, notification_factory
        )
        result = caller_factories.RearrangementCallerFactory(
            repo_factory, region_caller_factory, read_caller_factory
        )
        return result

    @pytest.fixture
    def reads_counter(self):
        return repositories.ProvidedSegmentCounter(2000)

    def test_get_rearrangements(self, factory, segment_repo, reads_counter):
        caller = factory.build(
            caller_factories.CallerType.MULTI,
            features=frozenset(
                [caller_factories.CallerFeature.WITH_PROVIDED_READ_COUNT]
            ),
            srpb_threshold=0.2,
            minimum_mapping_quality=1,
            total_number_of_reads=20000000,
        )
        observed = list(caller.get_rearrangements())
        assert len(observed) >= 2
        assert observed[0].region_pair.a.name == entities.RegionsName.CoreDUX4
        assert observed[0].region_pair.b.name == entities.RegionsName.IGH
        assert observed[1].region_pair.a.name == entities.RegionsName.ExtendedDUX4
        assert observed[1].region_pair.b.name == entities.RegionsName.IGH
        for rearranmgement in observed[2:]:
            assert rearranmgement.region_pair.a.name == entities.RegionsName.CoreDUX4
            assert rearranmgement.region_pair.b.name == entities.RegionsName.UNNAMED

    def test_get_rearrangements_named(self, factory, segment_repo, reads_counter):
        caller = factory.build(
            caller_factories.CallerType.MULTI,
            features=frozenset(
                [caller_factories.CallerFeature.WITH_PROVIDED_READ_COUNT]
            ),
            srpb_threshold=0.2,
            minimum_mapping_quality=1,
            total_number_of_reads=20000000,
        )
        observed = list(caller.get_rearrangements())
        assert len(observed) == 2
        assert observed[0].region_pair.a.name == entities.RegionsName.CoreDUX4
        assert observed[0].region_pair.b.name == entities.RegionsName.IGH
        assert observed[1].region_pair.a.name == entities.RegionsName.ExtendedDUX4
        assert observed[1].region_pair.b.name == entities.RegionsName.IGH
