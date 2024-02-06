import abc
import enum
import functools
from typing import FrozenSet, Iterable, List, Optional, Tuple

from ilmn.pelops import entities


class ReadQuery(enum.Enum):
    """Possible queries to filter placed segments"""

    is_duplicate = enum.auto()
    is_not_paired = enum.auto()
    is_proper_pair = enum.auto()
    is_not_proper_pair = enum.auto()
    is_qcfail = enum.auto()
    is_secondary = enum.auto()
    is_supplementary = enum.auto()
    is_unmapped = enum.auto()


class SegmentCounter(abc.ABC):
    @abc.abstractmethod
    def get_number_of_segments(self, exclude: Optional[List[ReadQuery]] = None) -> int:
        """Count total number of reads, excluding given list of `ReadQuery`"""


class ProvidedSegmentCounter(SegmentCounter):
    def __init__(self, provided_number_of_reads: int):
        self._provided_number_of_reads = provided_number_of_reads

    def get_number_of_segments(self, exclude: Optional[List[ReadQuery]] = None) -> int:
        return self._provided_number_of_reads


class CachedSegmentCounter(SegmentCounter):
    def __init__(self, counter: SegmentCounter):
        self.counter = counter

    def get_number_of_segments(self, exclude: Optional[List[ReadQuery]] = None) -> int:
        if exclude is None:
            hashable_exclude = None
        else:
            hashable_exclude = frozenset(exclude)
        return self._get_number_of_segments(hashable_exclude)

    @functools.lru_cache(maxsize=1024)
    def _get_number_of_segments(
        self, hashable_exclude: Optional[FrozenSet[ReadQuery]]
    ) -> int:
        """Results are cached to prevent recomputing very expensive step"""
        if hashable_exclude is None:
            exclude = None
        else:
            exclude = list(hashable_exclude)
        return self.counter.get_number_of_segments(exclude)


class PlacedSegmentRepository(SegmentCounter):
    @abc.abstractmethod
    def get(
        self,
        locations: FrozenSet[entities.GenomicRegion],
        exclude: List[ReadQuery],
        min_quality: int = 0,
    ) -> Iterable[entities.PlacedSegment]:
        """Retrieve segments placed in `location` and cache them"""

    @abc.abstractmethod
    def get_mate(
        self, read: entities.PlacedSegment
    ) -> Optional[entities.PlacedSegment]:
        """Retrieve paired segment of given read"""

    @abc.abstractmethod
    def get_mate_exact_position(self, read: entities.PlacedSegment) -> Tuple[str, int]:
        """Get contig name and starting position of read mate"""


class RegionRepository(abc.ABC):
    @abc.abstractmethod
    def get(self, name: entities.RegionsName) -> entities.CompoundRegion:
        """Get a CompoundRegion by name"""


class BuiltinRegionRepository(RegionRepository):
    _lookup = {
        entities.RegionsName.CoreDUX4: [
            ("chr4", 190020407, 190023665),
            ("chr4", 190066935, 190093279),
            ("chr4", 190172774, 190176845),
            ("chr10", 133663429, 133685936),
            ("chr10", 133739606, 133762125),
        ],
        entities.RegionsName.GRCh38: [
            ("chr1", 1, 248956422),
            ("chr2", 1, 242193529),
            ("chr3", 1, 198295559),
            ("chr4", 1, 190214555),
            ("chr5", 1, 181538259),
            ("chr6", 1, 170805979),
            ("chr7", 1, 159345973),
            ("chr8", 1, 145138636),
            ("chr9", 1, 138394717),
            ("chr10", 1, 133797422),
            ("chr11", 1, 135086622),
            ("chr12", 1, 133275309),
            ("chr13", 1, 114364328),
            ("chr14", 1, 107043718),
            ("chr15", 1, 101991189),
            ("chr16", 1, 90338345),
            ("chr17", 1, 83257441),
            ("chr18", 1, 80373285),
            ("chr19", 1, 58617616),
            ("chr20", 1, 64444167),
            ("chr21", 1, 46709983),
            ("chr22", 1, 50818468),
            ("chrX", 1, 156040895),
            ("chrY", 1, 57227415),
        ],
        entities.RegionsName.ExtendedDUX4: [
            ("chr4", 189967935, 190204560),
            ("chr10", 133564429, 133787422),
            ("chr3", 75667931, 75671185),
            ("chr5", 31248879, 31251987),
            ("chr9", 63816748, 63819462),
            ("chr12", 34207415, 34210675),
            ("chr12", 61599067, 61601843),
            ("chr16", 34134736, 34137792),
            ("chr16", 34140256, 34143849),
            ("chr20", 29317824, 29320579),
            ("chr20", 29323092, 29326049),
            ("chr20", 29409348, 29412600),
            ("chr20", 29447517, 29450306),
            ("chr20", 29877636, 29880363),
            ("chrY", 10170590, 10173725),
            ("chrY", 11305918, 11309181),
            ("chrY", 11313921, 11317187),
            ("chrY", 11320557, 11323823),
            ("chrY", 11331329, 11334595),
        ],
        entities.RegionsName.IGH: [("chr14", 105586937, 106879844)],
    }

    def get(self, name: entities.RegionsName) -> entities.CompoundRegion:
        regions = frozenset(
            entities.GenomicRegion(*region) for region in self._lookup[name]
        )
        result = entities.CompoundRegion(name=name, regions=regions)
        return result


class RegionRepoType(enum.Enum):
    BUILTIN = enum.auto()
    WITH_BLACKLIST = enum.auto()


class RegionRepositoryFactory(abc.ABC):
    @abc.abstractmethod
    def build(self, repo_type: RegionRepoType) -> RegionRepository:
        """Build a RegionRepository"""


class SegmentRepoFeature(enum.Enum):
    WITH_NOTIFICATION = enum.auto()
    BUILTIN = enum.auto()


class SegmentRepositoryFactory(abc.ABC):
    @abc.abstractmethod
    def build_counter(
        self,
        features: FrozenSet[SegmentRepoFeature],
        total_number_of_reads: Optional[int],
    ) -> SegmentCounter:
        """Build a SegmentCounter"""

    @abc.abstractmethod
    def build(
        self,
        features: FrozenSet[SegmentRepoFeature],
        total_number_of_reads: Optional[int],
    ) -> PlacedSegmentRepository:
        """Build a PlacedSegmentRepository"""
