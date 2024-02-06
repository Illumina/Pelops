"""Stubs, altertive simpler versions of real objects in the codebase"""

import dataclasses
from typing import FrozenSet, Iterable, List, Optional, Set, Tuple

from ilmn.pelops import entities, repositories


@dataclasses.dataclass(frozen=True)
class MockContent:
    next_reference_name: str
    next_reference_start: int


class SmallRegionRepository(repositories.BuiltinRegionRepository):
    locationA = frozenset([entities.GenomicRegion("chr1", 1, 100)])
    locationB = frozenset([entities.GenomicRegion("chr2", 2, 200)])
    locationC = frozenset([entities.GenomicRegion("chr3", 3, 300)])
    locationD = frozenset(
        [
            entities.GenomicRegion("chr1", 1, 300),
            entities.GenomicRegion("chr2", 1, 300),
            entities.GenomicRegion("chr3", 1, 300),
            entities.GenomicRegion("chr7", 1, 3000),
        ]
    )
    lookup = {
        entities.RegionsName.CoreDUX4: locationA,
        entities.RegionsName.IGH: locationB,
        entities.RegionsName.ExtendedDUX4: locationC,
        entities.RegionsName.GRCh38: locationD,
    }

    def get(self, name: entities.RegionsName) -> entities.CompoundRegion:
        regions = self.lookup[name]
        result = entities.CompoundRegion(name=name, regions=regions)
        return result


@dataclasses.dataclass(frozen=True)
class ExtraSegmentProperties:
    """Track extra read properties without changing the segment"""

    mapping_quality: int = 0


class StubPlacedSegmentRepository(repositories.PlacedSegmentRepository):
    """A stub repository of PlacedSegment. This stub repository complies with
    the interface, but it has extra methods to set up the PlacedSegment it
    stores and return"""

    def __init__(self):
        self._store: Set[Tuple[entities.PlacedSegment, ExtraSegmentProperties]] = set()

    def get_number_of_segments(
        self, exclude: Optional[List[repositories.ReadQuery]] = None
    ) -> int:
        result = sum([len(annotated_segment) for annotated_segment in self._store])
        return result

    def get(
        self,
        locations: FrozenSet[entities.GenomicRegion],
        exclude: List[repositories.ReadQuery],
        min_quality: int = 0,
    ) -> Iterable[entities.PlacedSegment]:
        for annotated_segment in self._store:
            for location in annotated_segment[0].location:
                if location in locations:
                    extra = annotated_segment[1]
                    if extra.mapping_quality >= min_quality:
                        yield annotated_segment[0]

    def get_mate(self, read: entities.PlacedSegment) -> entities.PlacedSegment:
        return entities.PlacedSegment("foo", frozenset(), entities.ReadOrder.ONE)

    def get_mate_exact_position(self, read: entities.PlacedSegment) -> Tuple[str, int]:
        """Get contig name and starting position of read mate"""
        if read.content is None:
            raise ValueError("invalid read type")
        mate_name = read.content.next_reference_name
        mate_start = read.content.next_reference_start + 1
        return mate_name, mate_start

    def add_read(self, read: entities.PlacedSegment, mapping_quality: int = 0) -> None:
        """Not required by the interface"""
        extra = ExtraSegmentProperties(mapping_quality=mapping_quality)
        self._store.add((read, extra))
