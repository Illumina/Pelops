import dataclasses
from typing import FrozenSet

import pytest

from ilmn.pelops import entities, repositories

A = tuple([entities.GenomicRegion(chrom="chr1", start=123, end=456)])
B = tuple([entities.GenomicRegion(chrom="chr2", start=333, end=444)])
C = tuple([entities.GenomicRegion(chrom="chrX", start=6660, end=7770)])

R1 = entities.ReadOrder.ONE
R2 = entities.ReadOrder.TWO


class TestGenomicRegion:
    overlap_data = [
        pytest.param([("chr1", 100, 200), ("chr1", 100, 200), True], id="identical"),
        pytest.param([("1", 100, 200), ("1", 150, 250), True], id="overlapping"),
        pytest.param([("1", 100, 200), ("2", 100, 200), False], id="different contigs"),
        pytest.param([("1", 100, 200), ("1", 444, 555), False], id="separate, sorted"),
        pytest.param(
            [("1", 444, 555), ("1", 100, 200), False], id="separate, not sorted"
        ),
        pytest.param([("1", 100, 200), ("1", 200, 230), True], id="one base overlap"),
        pytest.param([("1", 100, 200), ("1", 201, 230), False], id="contiguous"),
        pytest.param([("1", 100, 200), ("1", 1, 230), True], id="a inside b"),
        pytest.param([("1", 100, 200), ("1", 120, 180), True], id="b inside a"),
    ]

    @pytest.fixture(params=overlap_data)
    def test_case(self, request):
        return request.param[0], request.param[1], request.param[2]

    @pytest.fixture
    def a_region(self, test_case):
        a_region, another_region, expected = test_case
        result = entities.GenomicRegion(*a_region)
        return result

    @pytest.fixture
    def another_region(self, test_case):
        a_region, another_region, expected = test_case
        result = entities.GenomicRegion(*another_region)
        return result

    @pytest.fixture
    def expected(self, test_case):
        a_region, another_region, expected = test_case
        return expected

    def test_overlap(self, a_region, another_region, expected):
        assert a_region.overlaps(another_region) == expected

    def test_ovelap_special_cases(self):
        a_region = entities.GenomicRegion("chr1", 100, 200)
        another_object = ("chr1", 100, 200)
        assert not a_region.overlaps(another_object)
        assert not a_region.overlaps(None)

    def test_within(self):
        a = entities.GenomicRegion("chr1", 1, 100)
        b = entities.GenomicRegion("chr1", 1, 100)
        c = entities.GenomicRegion("chr1", 1, 200)
        assert a.within(b)
        assert a.within(c)
        d = entities.GenomicRegion("chr1", 1, 99)
        e = entities.GenomicRegion("chr1", 2, 100)
        f = entities.GenomicRegion("chr2", 1, 100)  # different chr
        assert not a.within(d)
        assert not a.within(e)
        assert not a.within(f)


class TestCompoundRegion:
    unnamed = entities.RegionsName.UNNAMED
    dux4 = entities.RegionsName.CoreDUX4
    region_a = entities.GenomicRegion("chr1", 100, 200)
    region_b = entities.GenomicRegion("chr1", 100, 150)
    region_c = entities.GenomicRegion("chr2", 100, 150)

    overlap_data = [
        pytest.param([(unnamed, []), (unnamed, []), False], id="empty regions"),
        pytest.param(
            [(unnamed, [region_a, region_c]), (unnamed, [region_b]), True],
            id="one overlapping region, same name",
        ),
        pytest.param(
            [(dux4, [region_a, region_c]), (unnamed, [region_b]), True],
            id="one overlapping region, different names",
        ),
        pytest.param(
            [(dux4, [region_a]), (dux4, [region_c]), False],
            id="no overlapping region, same name",
        ),
        pytest.param(
            [(dux4, [region_a]), (dux4, []), False], id="one empty region set, a"
        ),
        pytest.param(
            [(dux4, []), (dux4, [region_c]), False], id="one empty region set, b"
        ),
    ]

    @pytest.fixture(params=overlap_data)
    def test_case(self, request):
        return request.param[0], request.param[1], request.param[2]

    @pytest.fixture
    def a_region_set(self, test_case):
        a_region_set, _, expected = test_case
        name, regions = a_region_set
        result = entities.CompoundRegion(name, frozenset(regions))
        return result

    @pytest.fixture
    def another_region_set(self, test_case):
        _, another_region_set, expected = test_case
        name, regions = another_region_set
        result = entities.CompoundRegion(name, frozenset(regions))
        return result

    @pytest.fixture
    def expected(self, test_case):
        a_region_set, another_region_set, expected = test_case
        return expected

    def test_overlap(self, a_region_set, another_region_set, expected):
        assert a_region_set.overlaps(another_region_set) == expected

    def test_overlap_fails(self, an_unnamed_genomic_region_set):
        with pytest.raises(TypeError):
            assert an_unnamed_genomic_region_set.overlaps(self.dux4)

    def test_contain(self):
        genome = entities.CompoundRegion(
            entities.RegionsName.GRCh38,
            frozenset(
                [
                    entities.GenomicRegion("chr1", 1, 10000),
                    entities.GenomicRegion("chr2", 1, 2000),
                ]
            ),
        )
        a = entities.GenomicRegion("chr1", 2, 300)
        b = entities.GenomicRegion("chr2", 2, 2001)
        assert genome.contains(a)
        assert not genome.contains(b)

    def test_contain_fails(self, an_unnamed_genomic_region_set):
        with pytest.raises(TypeError):
            assert an_unnamed_genomic_region_set.contains(self.dux4)

    def test_is_hashable(self, a_region_set, another_region_set):
        A = entities.CompoundRegion(
            self.dux4, frozenset([self.region_a, self.region_c])
        )
        B = entities.CompoundRegion(
            self.unnamed, frozenset([self.region_a, self.region_c])
        )
        hashed = frozenset({A, B})
        assert len(hashed) == 2


@pytest.fixture
def read1():
    return entities.ReadOrder.ONE


@pytest.fixture
def read2():
    return entities.ReadOrder.TWO


@pytest.fixture
def region_a():
    return A


@pytest.fixture
def region_b():
    return B


def new_segment(name, region, read_order):
    result = entities.PlacedSegment(name, region, read_order)
    return result


class TestPlacedSegment:
    def test_equal(self, region_a):
        region = frozenset([entities.GenomicRegion("chr1", 1, 200)])
        s1 = entities.PlacedSegment("A", region, R1)
        s2 = entities.PlacedSegment("A", region, R1)
        s3 = entities.PlacedSegment("A", region, R2)
        assert s1 == s2
        assert s1 != s3

    def test_equal_fails(self, region_a):
        region = frozenset([entities.GenomicRegion("chr1", 1, 200)])
        s1 = entities.PlacedSegment("A", region, R1)
        with pytest.raises(TypeError):
            s1 == region


class TestRead:
    read_data = [
        pytest.param(
            [
                {"name": "ABCD", "region": A, "read_order": R1},
            ],
            False,
            id="only one segment",
        ),
        pytest.param(
            [
                {"name": "ABCD", "region": A, "read_order": R1},
                {"name": "ABCD", "region": B, "read_order": R1},
            ],
            True,
            id="typical split read",
        ),
        pytest.param(
            [
                {"name": "ABCD", "region": A, "read_order": R1},
                {"name": "ABCD", "region": A, "read_order": R1},
            ],
            False,
            id="two segment on the same region",
        ),
        pytest.param(
            [
                {"name": "ABCD", "region": A, "read_order": R1},
                {"name": "ABCD", "region": C, "read_order": R1},
            ],
            False,
            id="split read  across irrelevant locations",
        ),
        pytest.param(
            [
                {"name": "ABCD", "region": A, "read_order": R1},
                {"name": "ABCD", "region": B, "read_order": R1},
                {"name": "ABCD", "region": B, "read_order": R1},
            ],
            True,
            id="Three segments across two locations",
        ),
        pytest.param(
            [
                {"name": "ABCD", "region": A, "read_order": R1},
                {"name": "ABCD", "region": B, "read_order": R1},
                {"name": "ABCD", "region": C, "read_order": R1},
            ],
            True,
            id="Three segments across three locations",
        ),
    ]

    @pytest.mark.parametrize("segments, is_split", read_data)
    def test_is_split(self, segments, is_split, region_a, region_b):
        read = entities.Read()
        for segment in segments:
            read.allocate(new_segment(**segment))
        assert read.is_split_across(region_a, region_b) == is_split

    def test_allocate_raises_different_read_order(
        self, region_a, region_b, read1, read2
    ):
        read = entities.Read()
        segment1 = new_segment(name="ABCD", region=region_a, read_order=read1)
        segment2 = new_segment(name="ABCD", region=region_a, read_order=read2)
        read.allocate(segment1)
        with pytest.raises(entities.InvalidReadOrder):
            read.allocate(segment2)

    def test_allocate_split_raises_different_read_name(self, region_a, read1):
        read = entities.Read()
        segment1 = new_segment(name="ABCD", region=region_a, read_order=read1)
        segment2 = new_segment(name="EFG", region=region_a, read_order=read1)
        read.allocate(segment1)
        with pytest.raises(entities.InvalidReadName):
            read.allocate(segment2)

    def test_len(self, read1, region_a, region_b):
        segment1 = new_segment(name="ABCD", region=region_a, read_order=read1)
        segment2 = new_segment(name="ABCD", region=region_b, read_order=read1)
        read = entities.Read()
        assert len(read) == 0
        read.allocate(segment1)
        assert len(read) == 1
        read.allocate(segment2)
        assert len(read) == 2

    def test_get_locations(self, read1, region_a, region_b):
        segment1 = new_segment(name="ABCD", region=region_a, read_order=read1)
        segment2 = new_segment(name="ABCD", region=region_b, read_order=read1)
        read = entities.Read()
        read.allocate(segment1)
        assert read.get_locations() == set([region_a])
        read.allocate(segment2)
        assert read.get_locations() == set([region_a, region_b])


class TestReadPair:
    def test_has_both_mates(self, region_a, read1, read2):
        read_pair = entities.ReadPair()
        assert not read_pair.has_both_mates()

        placed_segment_r1 = new_segment("A00555", region_a, read1)
        read_pair.allocate(placed_segment_r1)
        assert not read_pair.has_both_mates()

        placed_segment_r2 = new_segment("A00555", region_a, read2)
        read_pair.allocate(placed_segment_r2)
        assert read_pair.has_both_mates()

    def test_allocate_fails_if_name_different(self, region_a, read1, read2):
        placed_segment_r1 = new_segment("A00555", region_a, read1)
        placed_segment_r2 = new_segment("ABCDE", region_a, read2)
        read_pair = entities.ReadPair()
        read_pair.allocate(placed_segment_r1)
        with pytest.raises(entities.InvalidReadName) as exc:
            read_pair.allocate(placed_segment_r2)
        assert str(exc.value) == "Impossible to allocate segment ABCDE to read A00555"

    has_split_read_data = [
        pytest.param(
            [
                {"name": "ABCD", "region": A, "read_order": R1},
                {"name": "ABCD", "region": B, "read_order": R1},
                {"name": "ABCD", "region": B, "read_order": R2},
            ],
            True,
            False,
            id="typical split read case",
        ),
        pytest.param(
            [
                {"name": "ABCD", "region": A, "read_order": R1},
                {"name": "ABCD", "region": C, "read_order": R1},
                {"name": "ABCD", "region": C, "read_order": R2},
            ],
            False,
            False,
            id="split read across irrelevant region",
        ),
        pytest.param(
            [
                {"name": "ABCD", "region": A, "read_order": R1},
                {"name": "ABCD", "region": C, "read_order": R1},
                {"name": "ABCD", "region": B, "read_order": R2},
            ],
            False,
            True,
            id="split read across irrelevant region with pair read",
        ),
        pytest.param(
            [
                {"name": "ABCD", "region": A, "read_order": R1},
                {"name": "ABCD", "region": B, "read_order": R1},
                {"name": "ABCD", "region": C, "read_order": R2},
                {"name": "ABCD", "region": B, "read_order": R2},
            ],
            True,
            False,
            id="split read across relevant and irrelevant region",
        ),
        pytest.param(
            [
                {"name": "ABCD", "region": A, "read_order": R1},
                {"name": "ABCD", "region": B, "read_order": R1},
                {"name": "ABCD", "region": B, "read_order": R2},
                {"name": "ABCD", "region": B, "read_order": R2},
            ],
            True,
            False,
            id="Four segments one split read",
        ),
        pytest.param(
            [
                {"name": "ABCD", "region": A, "read_order": R1},
                {"name": "ABCD", "region": B, "read_order": R1},
                {"name": "ABCD", "region": B, "read_order": R2},
                {"name": "ABCD", "region": A, "read_order": R2},
            ],
            True,
            False,
            id="Four segments two split reads",
        ),
        pytest.param(
            [
                {"name": "ABCD", "region": A, "read_order": R1},
            ],
            False,
            False,
            id="only one read",
        ),
        pytest.param(
            [
                {"name": "ABCD", "region": A, "read_order": R1},
                {"name": "ABCD", "region": A, "read_order": R1},
            ],
            False,
            False,
            id="Two segments from one read on the same region",
        ),
        pytest.param(
            [
                {"name": "ABCD", "region": A, "read_order": R1},
                {"name": "ABCD", "region": B, "read_order": R2},
            ],
            False,
            True,
            id="typical improper pair read",
        ),
        pytest.param(
            [
                {"name": "ABCD", "region": A, "read_order": R1},
                {"name": "ABCD", "region": A, "read_order": R2},
            ],
            False,
            False,
            id="typical proper pair read",
        ),
        pytest.param(
            [
                {"name": "ABCD", "region": A, "read_order": R1},
                {"name": "ABCD", "region": B, "read_order": R2},
                {"name": "ABCD", "region": B, "read_order": R2},
            ],
            False,
            True,
            id="paired read across two regions",
        ),
    ]

    @pytest.mark.parametrize("segments_data, is_split, is_paired", has_split_read_data)
    def test_has_split_read_and_improper_pair(
        self, segments_data, is_split, is_paired, region_a, region_b
    ):
        segments = [new_segment(**item) for item in segments_data]
        read_pair = entities.ReadPair()
        for segment in segments:
            read_pair.allocate(segment)
        assert read_pair.has_split_read_across(region_a, region_b) == is_split
        assert read_pair.has_improper_pair_across(region_a, region_b) == is_paired


class TestCompoundRegionPair:
    def test_get_names(self):
        region_repo = repositories.BuiltinRegionRepository()
        dux4 = region_repo.get(entities.RegionsName.CoreDUX4)
        igh = region_repo.get(entities.RegionsName.IGH)
        region_pair = entities.CompoundRegionPair(dux4, igh)
        observed = region_pair.get_names()
        assert observed == (entities.RegionsName.CoreDUX4, entities.RegionsName.IGH)
