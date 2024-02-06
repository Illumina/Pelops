from unittest import mock

import pytest

from ilmn.pelops import entities, repositories, selectors
from ilmn.pelops.callers import read_callers, region_callers
from tests import stubs


class TestContiguousRegionsCaller:
    find_contiguous_data = [
        pytest.param([["chr1", 1001, 2000]], [], id="one region"),
        pytest.param(
            [("chr1", 1001, 2000), ("chr1", 2000, 3000)],
            [[("chr1", 1001, 2000), ("chr1", 2000, 3000)]],
            id="two contiguous regions, overlapping one position",
        ),
        pytest.param(
            [
                ("chr1", 5001, 6000),
                ("chr1", 7001, 8000),
                ("chr1", 8001, 9000),
                ("chr1", 12001, 13001),
            ],
            [[("chr1", 7001, 8000), ("chr1", 8001, 9000)]],
            id="two contiguous regions, with non contiguous region preceeding and following",
        ),
        pytest.param(
            [("chr1", 1000, 1999), ("chr1", 2000, 3000)],
            [[("chr1", 1000, 1999), ("chr1", 2000, 3000)]],
            id="two contiguous regions, not overlapping",
        ),
        pytest.param(
            [("chr1", 1000, 1999), ("chrX", 2000, 3000)],
            [],
            id="two regions on different chromosomes, but with contiguous positions",
        ),
        pytest.param(
            (("chr1", 1001, 2000), ("chr1", 2001, 3000), ("chr1", 3001, 4000)),
            [[("chr1", 1001, 2000), ("chr1", 2001, 3000), ("chr1", 3001, 4000)]],
            id="three contiguous regions",
        ),
        pytest.param(
            (
                ("chr1", 1001, 2000),
                ("chr1", 2001, 3000),
                ("chr1", 10001, 11000),
                ("chr7", 5001, 6000),
                ("chr7", 6001, 7000),
            ),
            [
                [("chr1", 1001, 2000), ("chr1", 2001, 3000)],
                [("chr7", 5001, 6000), ("chr7", 6001, 7000)],
            ],
            id="two sets of contiguous regions",
        ),
        pytest.param(
            (
                ("chr1", 1001, 2000),
                ("chr1", 2001, 3000),
                ("chr7", 3001, 4000),
                ("chr7", 4001, 5000),
            ),
            [
                [("chr1", 1001, 2000), ("chr1", 2001, 3000)],
                [("chr7", 3001, 4000), ("chr7", 4001, 5000)],
            ],
            id="two sets of contiguous regions on different chromosomes",
        ),
    ]

    @pytest.mark.parametrize("provided, expected", find_contiguous_data)
    def test_find_contiguous(self, provided, expected):
        caller = region_callers.ContiguousRegionsCaller()
        provided_regions = [entities.GenomicRegion(*item) for item in provided]
        expected_regions = []
        for block in expected:
            expected_regions.append([entities.GenomicRegion(*item) for item in block])
        observed = list(caller.find_contiguous(provided_regions))
        assert observed == expected_regions


class TestRegionConsolidationService:
    consolidate_data = [
        pytest.param([["chr1", 1000, 2000]], [], id="one region once"),
        pytest.param(
            [["chr1", 5001, 6000]] * 3
            + [["chr1", 7001, 8000]] * 3
            + [["chr1", 8001, 9000]] * 2
            + [["chr1", 12001, 13001]] * 7,
            [["chr1", 5001, 6000], ["chr1", 7001, 9000], ["chr1", 12001, 13001]],
            id="two contiguous regions, with non contiguous region preceding and following",
        ),
        pytest.param(
            [["chr1", 5001, 6000]] * 3
            + [["chr1", 7001, 8000]] * 3
            + [["chr1", 8001, 9000]] * 1
            + [["chr1", 12001, 13001]] * 7,
            [["chr1", 5001, 6000], ["chr1", 7001, 8000], ["chr1", 12001, 13001]],
            id="one region with insufficient counts",
        ),
    ]

    @pytest.mark.parametrize("provided, expected", consolidate_data)
    def test_consolidate(self, provided, expected):
        provided_regions = [entities.GenomicRegion(*item) for item in provided]
        expected_regions = [entities.GenomicRegion(*item) for item in expected]
        service = region_callers.RegionConsolidationService()
        for candidate in provided_regions:
            service.add(candidate)
        service.consolidate()
        observed = list(service.get_regions())
        assert observed == expected_regions


class TestCandidateRegionCaller:
    locationA = frozenset([entities.GenomicRegion("chr4", 190066935, 190093279)])
    get_candidate_data = [
        pytest.param(
            (
                # provided
                [
                    ["chr1", 12340],
                    ["chr3", 12340],
                    ["chr1", 12350],
                    ["chr1", 44000],
                ],
                # expected
                [["chr1", 12001, 13000]],
                # selectors
                [],
            ),
            id="simple case",
        ),
        pytest.param(
            (
                # provided
                [
                    ["chr1", 12340],
                    ["chr1", 12350],
                    ["chr1", 13010],
                    ["chr1", 13020],
                    ["chr1", 13030],
                ],
                # expected
                [["chr1", 12001, 14000]],
                # selectors
                [],
            ),
            id="reads spanning two contiguous 1kb bins",
        ),
        pytest.param(
            (
                # provided
                [
                    ["chr1", 12340],
                    ["chr1", 12350],
                    ["chr1", 14010],
                    ["chr1", 14020],
                    ["chr1", 14030],
                ],
                # expected
                [["chr1", 12001, 13000], ["chr1", 14001, 15000]],
                # selectors
                [],
            ),
            id="reads spanning two disjoint 1kb bins",
        ),
        pytest.param(
            (
                # provided
                [
                    ["chr1", 12340],
                    ["chr1", 12350],
                    ["chr1", 14010],
                    ["chr1", 14020],
                    ["chr1", 14030],
                ],
                # expected
                [["chr1", 12001, 13000], ["chr1", 14001, 15000]],
                # selectors
                [{"include": [entities.RegionsName.GRCh38]}],
            ),
            id="reads spanning two disjoint 1kb bins, with genome selector",
        ),
        pytest.param(
            (
                # provided
                [
                    ["unplaced_contig", 12340],
                    ["unplaced_contig", 12350],
                    ["unplaced_contig", 14010],
                    ["unplaced_contig", 14020],
                    ["unplaced_contig", 14030],
                ],
                # expected
                [],
                # selectors
                [{"include": [entities.RegionsName.GRCh38]}],
            ),
            id="reads spanning two disjoint unplaced contigs, with genome selector",
        ),
    ]

    def build_mates(self, contig, position):
        result = entities.PlacedSegment(
            "ABCD",  # name does not matter, but required
            self.locationA,
            entities.ReadOrder.ONE,  # read order does not matter, but required
            stubs.MockContent(contig, position),
        )
        return result

    def build_expected(self, contig, start, end):
        result = entities.GenomicRegion(contig, start, end)
        return result

    @pytest.fixture(params=get_candidate_data)
    def one_test_case(self, request):
        """Allow parametrisation of expected and provided"""
        return request.param

    @pytest.fixture
    def reads_with_mate(self, one_test_case):
        provided, _, _ = one_test_case
        reads = [self.build_mates(*item) for item in provided]
        return reads

    @pytest.fixture
    def expected(self, one_test_case):
        _, expected, _ = one_test_case
        result = set([self.build_expected(*item) for item in expected])
        return result

    @pytest.fixture
    def all_selectors(self, one_test_case):
        _, _, all_selectors = one_test_case
        result = []
        for one_selector in all_selectors:
            for name in one_selector["include"]:
                selector = selectors.NamedRegionSelector(
                    name, repositories.BuiltinRegionRepository()
                )
                result.append(selector)
        return result

    @pytest.fixture
    def segment_repo(self, reads_with_mate):
        repo = stubs.StubPlacedSegmentRepository()
        for read in reads_with_mate:
            repo.add_read(read)
        return repo

    def test_get_candidate(self, segment_repo, expected, all_selectors):
        caller = region_callers.CandidateRegionCaller(
            segment_repo, selectors=all_selectors
        )
        provided = entities.CompoundRegion(
            name=entities.RegionsName.CoreDUX4, regions=self.locationA
        )
        observed = set(caller.get_candidates_regions(provided))
        assert observed == expected
