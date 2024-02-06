import abc
import dataclasses
import enum
from typing import (
    Any,
    Dict,
    FrozenSet,
    Iterable,
    Iterator,
    NamedTuple,
    Optional,
    Set,
    Tuple,
    Union,
)


@dataclasses.dataclass
class ClassifiedSegmentCount:
    """Number of observed segments by classification"""

    paired: int
    split: int
    spanning: int


class RegionsName(enum.Enum):
    """Possible names for a CompoundRegions"""

    CoreDUX4 = enum.auto()
    ExtendedDUX4 = enum.auto()
    GRCh38 = enum.auto()
    IGH = enum.auto()
    UNNAMED = enum.auto()
    BLACKLIST = enum.auto()


class GenomicRegion:
    """An immutable genomic region

    1-based genomic region: i.e. ("chr1", 100, 200) means from base 100 to 200
    of contig "chr1" included
    """

    def __init__(self, chrom: str, start: int, end: int):
        self._chrom = chrom
        self._start = start
        self._end = end

    def __eq__(self, other: object) -> bool:
        result = isinstance(other, type(self)) and (
            self.chrom == other.chrom
            and self.start == other.start
            and self.end == other.end
        )
        return result

    def __hash__(self) -> int:
        return hash((self._chrom, self._start, self._end))

    @property
    def chrom(self) -> str:
        return self._chrom

    @property
    def start(self) -> int:
        return self._start

    @property
    def end(self) -> int:
        return self._end

    def overlaps(self, other: object) -> bool:
        """True if self overlaps other. Contiguous regions do not overlap"""
        if not isinstance(other, type(self)) or self.chrom != other.chrom:
            return False

        first, second = sorted([self, other], key=lambda x: x.start)
        return first.end >= second.start

    def within(self, other: object) -> bool:
        """True if self is within other"""
        return (
            isinstance(other, type(self))
            and self.chrom == other.chrom
            and self.start >= other.start
            and self.end <= other.end
        )


class CompoundRegion:
    """A named, immutable set of GenomicRegions."""

    def __init__(self, name: RegionsName, regions: FrozenSet[GenomicRegion]):
        self._name = name
        self._regions = regions

    def __eq__(self, other: object) -> bool:
        result = isinstance(other, type(self)) and (
            self.name == other.name and self.regions == other.regions
        )
        return result

    def __hash__(self) -> int:
        return hash((self._name, self._regions))

    @property
    def name(self) -> RegionsName:
        return self._name

    @property
    def regions(self) -> FrozenSet[GenomicRegion]:
        return self._regions

    def overlaps(self, other: object) -> bool:
        """True if at least one region of self overlaps one region of other.
        Region name is irrelevant"""

        if not isinstance(other, type(self)):
            raise TypeError(
                f"Incompatible type: expected {type(self)}, got {type(other)}"
            )

        for other_region in other.regions:
            if any([region.overlaps(other_region) for region in self.regions]):
                return True

        return False

    def contains(self, other: object) -> bool:
        """True if `other` is within regions of `self`"""
        if isinstance(other, GenomicRegion):
            is_within_genome = any([other.within(region) for region in self._regions])
            return is_within_genome
        raise TypeError(
            f"Incompatible type: expected {type(GenomicRegion)}, got {type(other)}"
        )


class CompoundRegionPair:
    """A pair of CompoundRegions."""

    def __init__(self, a: CompoundRegion, b: CompoundRegion):
        self.__a = a
        self.__b = b

    def __eq__(self, other: object) -> bool:
        result = (
            isinstance(other, type(self)) and self.a == other.a and self.b == other.b
        )
        return result

    def __iter__(self) -> Iterator[CompoundRegion]:
        yield self.__a
        yield self.__b

    @property
    def a(self) -> CompoundRegion:
        return self.__a

    @property
    def b(self) -> CompoundRegion:
        return self.__b

    def get_names(self) -> Tuple[RegionsName, RegionsName]:
        result = self.a.name, self.b.name
        return result


class ReadOrder(enum.Enum):
    """Either read one or read two"""

    ONE = enum.auto()
    TWO = enum.auto()


class PlacedSegment:
    """A Placed Segment that carries over pysam content"""

    def __init__(
        self,
        read_name: str,
        location: FrozenSet[GenomicRegion],
        read_order: ReadOrder,
        content: Any = None,
    ):
        self.__read_name = read_name
        self.__location = location
        self.__read_order = read_order
        self.__content = content

    @property
    def read_name(self) -> str:
        return self.__read_name

    @property
    def location(self) -> FrozenSet[GenomicRegion]:
        return self.__location

    @property
    def read_order(self) -> ReadOrder:
        return self.__read_order

    @property
    def content(self) -> Any:
        """Read content not used by the business logic but required downstream"""
        return self.__content

    def __eq__(self, other: object) -> bool:
        if isinstance(other, type(self)):
            attributes_to_compare = ("read_name", "location", "read_order", "content")
            for attribute in attributes_to_compare:
                if getattr(self, attribute) != getattr(other, attribute):
                    return False
            return True
        else:
            raise TypeError(
                f"Incompatible type: expected {type(self)}, got {type(other)}"
            )

    def __hash__(self) -> int:
        """Uniquely defined by Pysam Segment content"""
        return hash(self.content)


class InvalidReadName(ValueError):
    def __init__(self, read: Union["Read", "ReadPair"], segment: "PlacedSegment"):
        message = f"Impossible to allocate segment {segment.read_name} to read {read.get_name()}"
        super().__init__(message)


class InvalidReadOrder(ValueError):
    def __init__(self, read: Union["Read", "ReadPair"], segment: "PlacedSegment"):
        message = f"Impossible to allocate segment {segment.read_order} to read {read.get_name()}"
        super().__init__(message)


class Read:
    """Models a single read

    A `Read` is created empty and `PlacedSegment` are then allocated
    """

    def __init__(self) -> None:
        self._segments: Set[PlacedSegment] = set()
        self._read_order: Optional["ReadOrder"] = None
        self._read_name: Optional[str] = None

    def __len__(self) -> int:
        return len(self._segments)

    def get_name(self) -> Optional[str]:
        return self._read_name

    def get_locations(self) -> Set[FrozenSet[GenomicRegion]]:
        """Get the locations of the `PlacedSegment` assigned to this read"""
        result = set([segment.location for segment in self._segments])
        return result

    def is_split_across(
        self, a: FrozenSet[GenomicRegion], b: FrozenSet[GenomicRegion]
    ) -> bool:
        if len(self._segments) > 1:
            locations = set([segment.location for segment in self._segments])
            if set([a, b]).issubset(locations):
                return True
        return False

    def allocate(self, placed_segment: PlacedSegment) -> None:
        """Allocate the `PlacedSegment` to the read"""
        if self._read_order is None:
            self._read_order = placed_segment.read_order
        elif self._read_order != placed_segment.read_order:
            raise InvalidReadOrder(self, placed_segment)

        if self._read_name is None:
            self._read_name = placed_segment.read_name
        elif self._read_name != placed_segment.read_name:
            raise InvalidReadName(self, placed_segment)
        self._segments.add(placed_segment)

    def get_segments(self) -> Iterable[PlacedSegment]:
        for segment in self._segments:
            yield segment


class ReadPair:
    """Models a read pair.

    A `ReadPair` is created empty and `PlacedSegment` are then allocated
    """

    def __init__(self) -> None:
        self._pairs: Dict[ReadOrder, Read] = {
            ReadOrder.ONE: Read(),
            ReadOrder.TWO: Read(),
        }
        self._read_name: Optional[str] = None

    def __iter__(self) -> Iterator[Read]:
        yield self._pairs[ReadOrder.ONE]
        yield self._pairs[ReadOrder.TWO]

    def allocate(self, placed_segment: PlacedSegment) -> None:
        """Allocate the placed segment to the correct read of the read pair"""
        if self.get_name() is None:
            self._read_name = placed_segment.read_name
        elif self.get_name() != placed_segment.read_name:
            raise InvalidReadName(self, placed_segment)
        read = self._pairs[placed_segment.read_order]
        read.allocate(placed_segment)

    def get_name(self) -> Optional[str]:
        return self._read_name

    def has_split_read_across(
        self, a: FrozenSet[GenomicRegion], b: FrozenSet[GenomicRegion]
    ) -> bool:
        result = any([read.is_split_across(a, b) for read in self])
        return result

    def has_improper_pair_across(
        self, a: FrozenSet[GenomicRegion], b: FrozenSet[GenomicRegion]
    ) -> bool:
        has_read_across_two_regions = self._has_reads_across(a, b)
        has_split_reads = self.has_split_read_across(a, b)
        result = has_read_across_two_regions and not has_split_reads
        return result

    def has_both_mates(self) -> bool:
        """True if both reads have at least one aligned segment"""
        result = all([len(read) > 0 for read in self._pairs.values()])
        return result

    def get_segments(self) -> Iterable[PlacedSegment]:
        for read in self._pairs.values():
            for segment in read.get_segments():
                yield segment

    def _has_reads_across(
        self, a: FrozenSet[GenomicRegion], b: FrozenSet[GenomicRegion]
    ) -> bool:
        target = set([a, b])
        r1_location, r2_location = [
            target.intersection(read.get_locations()) for read in self
        ]
        return target == r1_location.union(r2_location)


@dataclasses.dataclass
class Rearrangement:
    region_pair: CompoundRegionPair
    counts: ClassifiedSegmentCount
    srpb: float
    segments: Set[PlacedSegment]
