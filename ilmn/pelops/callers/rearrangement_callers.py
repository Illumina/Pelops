import abc
from typing import Iterable, Sequence, Set, Tuple

from ilmn.pelops import entities, repositories, selectors
from ilmn.pelops.callers import read_callers, region_pair_callers


class RearrangementCaller(abc.ABC):
    @abc.abstractmethod
    def get_rearrangements(self) -> Iterable[entities.Rearrangement]:
        """Get rearrangement evidence"""


class RegionRearrangementCaller(RearrangementCaller):
    """Finds rearrangements evidence across a region pairs"""

    def __init__(
        self,
        read_caller: read_callers.ReadsCaller,
        region_pair_caller: region_pair_callers.CompoundRegionPairCaller,
        reads_counter: repositories.SegmentCounter,
    ):
        self._read_caller = read_caller
        self._region_pair_caller = region_pair_caller
        self._reads_counter = reads_counter

    def get_rearrangements(self) -> Iterable[entities.Rearrangement]:
        region_pairs = self._region_pair_caller.get_compound_region_pairs()
        for region_pair in region_pairs:
            yield self._get_one_rearrangement(region_pair)

    def _get_one_rearrangement(
        self, region_pair: entities.CompoundRegionPair
    ) -> entities.Rearrangement:
        self._read_caller.detect_reads_spanning_regions(
            region_pair.a.regions, region_pair.b.regions
        )
        segments = set(self._read_caller.get_segments_of_spanning_reads())
        counts = self._read_caller.get_segment_count()
        srpb = (
            counts.spanning
            * 1_000_000_000
            / self._reads_counter.get_number_of_segments()
        )
        result = entities.Rearrangement(region_pair, counts, srpb, segments)
        return result


class MultiRearrangementCaller(RearrangementCaller):
    """A RearrangementCaller made of RearrangementCallers"""

    def __init__(self, callers: Sequence[RearrangementCaller]):
        self._callers = callers

    def get_rearrangements(self) -> Iterable[entities.Rearrangement]:
        for caller in self._callers:
            for rearrangement in caller.get_rearrangements():
                yield rearrangement


class SelectableRearrangementCaller(RearrangementCaller):
    def __init__(
        self, caller: RearrangementCaller, selector: selectors.RearrangementSelector
    ):
        self._caller = caller
        self._selector = selector

    def get_rearrangements(self) -> Iterable[entities.Rearrangement]:
        for item in self._caller.get_rearrangements():
            if self._selector(item):
                yield item


class RearrangementSorter(abc.ABC):
    """Sort Rearrangements"""

    @abc.abstractmethod
    def __call__(
        self, rearrangements: Iterable[entities.Rearrangement]
    ) -> Iterable[entities.Rearrangement]:
        """Sort Rearrangements"""


class SortedRearrangementCaller(RearrangementCaller):
    def __init__(self, caller: RearrangementCaller, sorter: RearrangementSorter):
        self._caller = caller
        self._sorter = sorter

    def get_rearrangements(self) -> Iterable[entities.Rearrangement]:
        for rearrangement in self._sorter(self._caller.get_rearrangements()):
            yield rearrangement
