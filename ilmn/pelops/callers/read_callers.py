import abc
from typing import FrozenSet, Iterable, List

from ilmn.pelops import entities, repositories, stores


class ReadsCaller(abc.ABC):
    @abc.abstractmethod
    def detect_reads_spanning_regions(
        self,
        a_regions: FrozenSet[entities.GenomicRegion],
        b_regions: FrozenSet[entities.GenomicRegion],
    ) -> None:
        """Detect reads spanning regions A and B"""

    @abc.abstractmethod
    def get_segment_count(self) -> entities.ClassifiedSegmentCount:
        """Provide detected counts classified segments"""

    @abc.abstractmethod
    def get_segments_of_spanning_reads(self) -> Iterable[entities.PlacedSegment]:
        """Retrieve PlacedSegments from reads spanning across the regions"""


class SpanningReadsCaller(ReadsCaller):
    """Find paired or split reads spanning across two regions"""

    def __init__(
        self,
        placed_segments_repository: repositories.PlacedSegmentRepository,
        reads_to_exclude: List[repositories.ReadQuery],
        minimum_mapping_quality: int = 0,
    ):
        self._placed_segment_repo = placed_segments_repository
        self._reads_to_exclude = reads_to_exclude
        self._minimum_mapping_quality = minimum_mapping_quality
        self._store = stores.PairedReadStore()

    def detect_reads_spanning_regions(
        self,
        a_regions: FrozenSet[entities.GenomicRegion],
        b_regions: FrozenSet[entities.GenomicRegion],
    ) -> None:
        """Detect reads spanning regions A and B, storing them into `PairReadsStore`"""
        self._store.clear()
        for placed_segment in self._placed_segment_repo.get(
            a_regions, exclude=self._reads_to_exclude, min_quality=0
        ):
            try:
                read_pair = self._store.get_read_pair(placed_segment.read_name)
            except stores.UnknowReadPairError:
                read_pair = entities.ReadPair()
            read_pair.allocate(placed_segment)
            self._store.update(read_pair)

        for placed_segment in self._placed_segment_repo.get(
            b_regions,
            exclude=self._reads_to_exclude,
            min_quality=self._minimum_mapping_quality,
        ):
            if placed_segment.read_name in self._store.get_read_names():
                read_pair = self._store.get_read_pair(placed_segment.read_name)
                read_pair.allocate(placed_segment)
                self._store.update(read_pair)

        for read_pair in self._store:
            if read_pair.has_split_read_across(a_regions, b_regions):
                self._store.mark_as_split(read_pair)
                self._store.mark_as_spanning(read_pair)
            if read_pair.has_improper_pair_across(a_regions, b_regions):
                self._store.mark_as_paired(read_pair)
                self._store.mark_as_spanning(read_pair)

    def get_segment_count(self) -> entities.ClassifiedSegmentCount:
        return self._store.get_segment_count()

    def get_segments_of_spanning_reads(self) -> Iterable[entities.PlacedSegment]:
        """Get PlacedSegment that consititute a spanning ReadPair from all spanning ReadPairs"""
        for read_pair in self._store.get_spanning_reads():
            for segment in read_pair.get_segments():
                yield segment
                if not read_pair.has_both_mates():
                    # NOTE: this is a rare event. This happens when the read pair
                    # has a split read, and the mate was not retrieved (e.g. is
                    # unmapped, or mapped outside the pair of CompoundRegion, or
                    # did not pass filters...)
                    mate = self._placed_segment_repo.get_mate(segment)
                    if mate is not None:
                        yield mate
