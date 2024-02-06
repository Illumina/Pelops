import pytest

from ilmn.pelops import entities, stores
from ilmn.pelops.callers import caller_factories, read_callers
from tests import stubs

# fixtures placed_segments_repository is located in conftest.py

features = [
    frozenset(),
    frozenset([caller_factories.CallerFeature.WITH_NOTIFICATIONS]),
]


@pytest.fixture(params=features)
def read_caller(read_caller_factory, request):
    # behaviour is the same whether we have notifications or not
    features = request.param
    read_caller_factory.reads_to_exclude = []
    result = read_caller_factory.build(
        minimum_mapping_quality=0, features=features, total_number_of_reads=None
    )
    return result


@pytest.fixture
def core_dux4_regions():
    result = tuple([entities.GenomicRegion("chr4", 190066935, 190093279)])
    return result


@pytest.fixture
def igh_regions():
    result = tuple([entities.GenomicRegion("chr14", 105586937, 106879844)])
    return result


class TestNotifyReadsCaller:
    def test_detect_reads_spanning_regions(
        self, read_caller_factory, core_dux4_regions, igh_regions
    ):
        features = frozenset([caller_factories.CallerFeature.WITH_NOTIFICATIONS])
        read_caller = read_caller_factory.build(
            minimum_mapping_quality=0, features=features, total_number_of_reads=None
        )
        read_caller.detect_reads_spanning_regions(core_dux4_regions, igh_regions)
        notification_service = read_caller.notification_service
        observed = notification_service.get_notifications()
        expected = ["Finding evidence for rearrangement."]
        assert expected == observed


class TestSpanningReadsCaller:
    def test_get_segment_count(self, read_caller, core_dux4_regions, igh_regions):
        read_caller.detect_reads_spanning_regions(core_dux4_regions, igh_regions)
        expected = entities.ClassifiedSegmentCount(paired=1, split=0, spanning=1)
        observed = read_caller.get_segment_count()
        assert expected == observed

    def test_get_segments_of_spanning_reads(
        self, read_caller, core_dux4_regions, igh_regions
    ):
        read_caller.detect_reads_spanning_regions(core_dux4_regions, igh_regions)
        observed = list(read_caller.get_segments_of_spanning_reads())
        assert len(observed) == 2

    @pytest.fixture
    def unnamed_region(self):
        result = tuple([entities.GenomicRegion("chr9", 1, 1000)])
        return result

    @pytest.fixture
    def segment_repo(self, core_dux4_regions, unnamed_region):
        R1 = entities.ReadOrder.ONE
        R2 = entities.ReadOrder.TWO
        presegments = [
            # name, location, order, mapping_quality
            ("foo", core_dux4_regions, R1, 0),
            ("foo", unnamed_region, R2, 0),
            ("bar", core_dux4_regions, R1, 0),
            ("bar", unnamed_region, R2, 1),
            ("baaz", core_dux4_regions, R1, 1),
            ("baaz", unnamed_region, R2, 0),
            ("meah", core_dux4_regions, R1, 1),
            ("meah", unnamed_region, R2, 1),
        ]

        repo = stubs.StubPlacedSegmentRepository()
        for item in presegments:
            segment = entities.PlacedSegment(item[0], item[1], item[2])
            repo.add_read(segment, item[3])
        return repo

    def test_get_segment_count_check_filtering(
        self, segment_repo, core_dux4_regions, unnamed_region
    ):
        read_caller = read_callers.SpanningReadsCaller(
            segment_repo, reads_to_exclude=[], minimum_mapping_quality=0
        )
        read_caller.detect_reads_spanning_regions(core_dux4_regions, unnamed_region)
        observed = read_caller.get_segment_count()
        expected = entities.ClassifiedSegmentCount(paired=4, split=0, spanning=4)
        assert observed == expected
        read_caller = read_callers.SpanningReadsCaller(
            segment_repo, reads_to_exclude=[], minimum_mapping_quality=1
        )
        read_caller.detect_reads_spanning_regions(core_dux4_regions, unnamed_region)
        observed = read_caller.get_segment_count()
        expected = entities.ClassifiedSegmentCount(paired=2, split=0, spanning=2)
        assert observed == expected
