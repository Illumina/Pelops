import itertools
from typing import Iterable, Iterator, List, Tuple

from ilmn.pelops import entities, repositories, selectors, stores

GenomicRegion = entities.GenomicRegion


def pairwise(
    iterable: Iterable[GenomicRegion],
) -> Iterator[Tuple[GenomicRegion, GenomicRegion]]:
    """Return successive overlapping pairs taken from iterable."""
    # part of itertools only from 3.10. Using their documented equivalent behaviour
    # https://docs.python.org/3/library/itertools.html#itertools.pairwise
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


class ContiguousRegionsCaller:
    """Find contiguous GenomicRegions"""

    def __init__(self) -> None:
        self._current_block_of_contiguous_regions: List[GenomicRegion] = []

    def find_contiguous(
        self, regions: Iterable[GenomicRegion]
    ) -> Iterable[List[GenomicRegion]]:
        """Go through an iterable of sorted `GenomicRegion` and returns
        an iterable of blocks of two or more contiguous `GenomicRegion`"""

        for previous, current in pairwise(regions):
            if self._are_contiguous(previous, current):
                self._add_to_block(previous, current)
            elif self._current_block_of_contiguous_regions:
                yield self._current_block_of_contiguous_regions
                self._current_block_of_contiguous_regions = []

        # at the end we want to return contiguous regions
        if len(self._current_block_of_contiguous_regions) > 1:
            yield self._current_block_of_contiguous_regions

    def _add_to_block(self, previous: GenomicRegion, current: GenomicRegion) -> None:
        if not self._current_block_of_contiguous_regions:
            self._current_block_of_contiguous_regions.append(previous)
        self._current_block_of_contiguous_regions.append(current)

    def _are_contiguous(self, a: GenomicRegion, b: GenomicRegion) -> bool:
        return a.chrom == b.chrom and a.end + 1 >= b.start


class RegionConsolidationService:
    """Merge or remove regions"""

    low_pass_filter = 2

    def __init__(self) -> None:
        self._store = stores.RegionStore()
        self._contiguous_region_caller = ContiguousRegionsCaller()

    def add(self, candidate: GenomicRegion) -> None:
        self._store.add(candidate)

    def get_regions(self) -> Iterable[GenomicRegion]:
        for region, count in self._store.get_counted_regions():
            yield region

    def consolidate(self) -> None:
        self._remove_low_count_regions()
        self._merge_consecutive_regions()

    def _merge_consecutive_regions(self) -> None:
        """Replace contiguous regions with a single merged region"""

        regions = self._store.get_regions()
        region_blocks = self._contiguous_region_caller.find_contiguous(regions)
        for contiguous_regions in region_blocks:
            merged_region = self._merge(contiguous_regions)
            self._store.add(merged_region)
            for region in contiguous_regions:
                self._store.delete(region)

    def _merge(self, group: List[GenomicRegion]) -> GenomicRegion:
        previous = group[0]
        for current in group[1:]:
            previous = entities.GenomicRegion(
                previous.chrom, previous.start, current.end
            )
        return previous

    def _remove_low_count_regions(self) -> None:
        for region, count in self._store.get_counted_regions():
            if count < self.low_pass_filter:
                self._store.delete(region)


class CandidateRegionCaller:
    """Find candidate mate GenomicRegion(s) of a given CompoundRegion.

    Given a CompoundRegion, finds candidate GenomicRegion(s) where a
    rearrangement might be occurring, based on evidence from placed segments"""

    region_size = 1000
    reads_to_exclude = [
        repositories.ReadQuery.is_duplicate,
        repositories.ReadQuery.is_not_paired,
        repositories.ReadQuery.is_proper_pair,
        repositories.ReadQuery.is_qcfail,
    ]

    def __init__(
        self,
        segment_repo: repositories.PlacedSegmentRepository,
        selectors: List[selectors.RegionSelector],
    ):
        self._segment_repo = segment_repo
        self._candidates = RegionConsolidationService()
        self._selectors = selectors

    def get_candidates_regions(
        self, origin: entities.CompoundRegion
    ) -> Iterable[GenomicRegion]:
        for segment in self._segment_repo.get(
            origin.regions, exclude=self.reads_to_exclude
        ):
            chrom, pos = self._segment_repo.get_mate_exact_position(segment)
            start, end = self._get_boundaries(pos)
            region = GenomicRegion(chrom, start, end)
            if all([selector(region) for selector in self._selectors]):
                self._candidates.add(region)
        self._candidates.consolidate()
        candidates = set(self._candidates.get_regions())
        return iter(candidates)

    def _get_boundaries(self, pos: int) -> Tuple[int, int]:
        low = int(pos / self.region_size) * self.region_size
        high = low + self.region_size
        return low + 1, high
