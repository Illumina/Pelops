import dataclasses
from typing import FrozenSet

import pytest

from ilmn.pelops import entities, stores


class TestPairedReadStore:
    @pytest.fixture
    def store(self):
        result = stores.PairedReadStore()
        return result

    @pytest.fixture
    def paired_segments(self):
        region_a = frozenset([entities.GenomicRegion(chrom="chr1", start=123, end=456)])
        region_b = frozenset(
            [entities.GenomicRegion(chrom="chr10", start=345, end=456)]
        )
        segment1 = entities.PlacedSegment("ABCD", region_a, entities.ReadOrder.ONE)
        segment2 = entities.PlacedSegment("ABCD", region_b, entities.ReadOrder.TWO)
        return (segment1, segment2)

    def test_get_update_read_pair(self, store, paired_segments):
        # if there is no known read with a certain name raise an error
        with pytest.raises(stores.UnknowReadPairError):
            store.get_read_pair("foobar")

        # adding an unnamed read_pair is not allowed
        read_pair = entities.ReadPair()
        with pytest.raises(stores.UnnamedReadPairError):
            store.update(read_pair)

        segment1, segment2 = paired_segments
        read_pair.allocate(segment1)
        store.update(read_pair)
        observed = store.get_read_pair("ABCD")
        assert observed == read_pair

    def test_mark_as(self, store, paired_segments):
        # if we didn't add any read_pair, all counts are zeros
        observed = store.get_segment_count()
        assert observed == entities.ClassifiedSegmentCount(
            paired=0, split=0, spanning=0
        )

        segment1, segment2 = paired_segments
        read_pair = entities.ReadPair()
        read_pair.allocate(segment1)
        store.update(read_pair)

        read = store.get_read_pair(segment1.read_name)
        store.mark_as_split(read)
        observed = store.get_segment_count()
        assert observed == entities.ClassifiedSegmentCount(
            paired=0, split=1, spanning=0
        )

        store.mark_as_paired(read)
        observed = store.get_segment_count()
        assert observed == entities.ClassifiedSegmentCount(
            paired=1, split=1, spanning=0
        )

        store.mark_as_spanning(read)
        observed = store.get_segment_count()
        assert observed == entities.ClassifiedSegmentCount(
            paired=1, split=1, spanning=1
        )


class TestRegionStore:
    roundtrip_data = [
        pytest.param(
            [("chr1", 1000, 2000)], [(("chr1", 1000, 2000), 1)], id="just one region"
        ),
        pytest.param(
            [("chr1", 1000, 2000), ("chr2", 1000, 2000)],
            [(("chr1", 1000, 2000), 1), (("chr2", 1000, 2000), 1)],
            id="two sorted regions",
        ),
        pytest.param(
            [("chr2", 1000, 2000), ("chr1", 1000, 2000)],
            [(("chr1", 1000, 2000), 1), (("chr2", 1000, 2000), 1)],
            id="two unsorted regions. Different chromosomes",
        ),
        pytest.param(
            [("chr2", 1000, 2000), ("chr1", 1000, 2000), ("chr2", 0, 1000)],
            [
                (("chr1", 1000, 2000), 1),
                (("chr2", 0, 1000), 1),
                (("chr2", 1000, 2000), 1),
            ],
            id="three unsorted regions. Different chromosomes and starts",
        ),
        pytest.param(
            [("chr2", 1000, 2000), ("chr2", 1000, 2000)],
            [(("chr2", 1000, 2000), 2)],
            id="one repeted region",
        ),
        pytest.param(
            [
                ("chr2", 1000, 2000),
                ("chr2", 0, 1000),
                ("chr2", 0, 1000),
                ("chr1", 1000, 2000),
                ("chr2", 0, 1000),
            ],
            [
                (("chr1", 1000, 2000), 1),
                (("chr2", 0, 1000), 3),
                (("chr2", 1000, 2000), 1),
            ],
            id="unsorted and repeted regions",
        ),
    ]

    @pytest.mark.parametrize("provided, expected", roundtrip_data)
    def test_add_get_roundtrip(self, provided, expected):
        provided_regions = [entities.GenomicRegion(*item) for item in provided]
        expected_counted_regions = [
            (entities.GenomicRegion(*item), counts) for item, counts in expected
        ]
        store = stores.RegionStore()
        for region in provided_regions:
            store.add(region)
        observed = list(store.get_counted_regions())
        assert observed == expected_counted_regions
