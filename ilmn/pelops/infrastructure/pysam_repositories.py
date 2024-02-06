# mypy: no_warn_unused_ignores
import abc
import enum
import functools
import pathlib
from typing import FrozenSet, Iterable, Iterator, List, Optional, Tuple

import pysam

from ilmn.pelops import entities, notifications, repositories
from ilmn.pelops.callers import caller_factories


class SamFlag(enum.IntEnum):
    """All known pysam flags and their numeric representation for sam/bam files"""

    is_paired = 0x1
    is_proper_pair = 0x2
    is_unmapped = 0x4
    mate_is_unmapped = 0x8
    is_reverse = 0x10
    mate_is_reverse = 0x20
    is_read1 = 0x40
    is_read2 = 0x80
    is_secondary = 0x100
    is_qcfail = 0x200
    is_duplicate = 0x400
    is_supplementary = 0x800


class SegmentSelector(abc.ABC):
    @abc.abstractmethod
    def __call__(self, segment: pysam.AlignedSegment) -> bool:
        """Select the segment"""


class PysamPropertyError(TypeError):
    def __init__(self, name: str, value: object, expected: str):
        self.name = name
        self.value = value
        self.expected = expected

    def __str__(self) -> str:
        property_name = f"pysam.AlignedSegment.{self.name}"
        observed = type(self.value).__name__
        return f"{property_name} returned unexpected type {observed}. Expected {self.expected}"


class BooleanSegmentSelector(SegmentSelector):
    """Safely select a pysam.AlignedSegment on property supposed to be boolean"""

    def __init__(self, name: str):
        self.__name = name

    def __call__(self, segment: pysam.AlignedSegment) -> bool:
        """Select on attribute name. Raise a PysamPropertyError if the value is not a boolean"""
        value = getattr(segment, self.__name)
        if not isinstance(value, bool):
            raise PysamPropertyError(self.__name, value, str(bool))
        return value


class SegmentDeSelector(SegmentSelector):
    """Invert the behaviour of a SegmentSelector"""

    def __init__(self, selector: SegmentSelector):
        self.selector = selector

    def __call__(self, segment: pysam.AlignedSegment) -> bool:
        return not self.selector(segment)


class AnySegmentSelector(SegmentSelector):
    """Selector on Any ReadQuery property"""

    def __init__(self, queries: Iterable[repositories.ReadQuery]):
        self._queries = queries
        self.__selector_factory = SelectorFactory()

    def __call__(self, segment: pysam.AlignedSegment) -> bool:
        """If any of the queries selects the segment, the segment is selected"""
        selectors = [self.__selector_factory.build(query) for query in self._queries]
        result = any([selector(segment) for selector in selectors])
        return result


class SelectorFactory:
    """Build SegmentSelector depending on ReadQuery"""

    def build(self, query: repositories.ReadQuery) -> SegmentSelector:
        result: SegmentSelector
        selector_lookup = {
            repositories.ReadQuery.is_duplicate: "is_duplicate",
            repositories.ReadQuery.is_qcfail: "is_qcfail",
            repositories.ReadQuery.is_proper_pair: "is_proper_pair",
            repositories.ReadQuery.is_not_paired: "is_paired",
        }
        deselectors = {repositories.ReadQuery.is_not_paired}
        result = BooleanSegmentSelector(selector_lookup[query])
        if query in deselectors:
            result = SegmentDeSelector(result)
        return result


class MapQSelector(SegmentSelector):
    def __init__(self, minimum_quality: int):
        self._minimum_quality = minimum_quality

    def __call__(self, segment: pysam.AlignedSegment) -> bool:
        """True if segment has at least minimum_quality"""
        mapping_quality = segment.mapping_quality
        if not isinstance(mapping_quality, (int, float)):
            raise PysamPropertyError("mapping_quality", mapping_quality, "number")
        result = segment.mapping_quality >= self._minimum_quality
        return result


class BamFileSegmentCounter(repositories.SegmentCounter):
    default_filters = [
        repositories.ReadQuery.is_unmapped,
        repositories.ReadQuery.is_secondary,
        repositories.ReadQuery.is_qcfail,
        repositories.ReadQuery.is_duplicate,
        repositories.ReadQuery.is_supplementary,
    ]

    lookup = {
        repositories.ReadQuery.is_unmapped: SamFlag.is_unmapped,
        repositories.ReadQuery.is_secondary: SamFlag.is_secondary,
        repositories.ReadQuery.is_qcfail: SamFlag.is_qcfail,
        repositories.ReadQuery.is_duplicate: SamFlag.is_duplicate,
        repositories.ReadQuery.is_supplementary: SamFlag.is_supplementary,
        repositories.ReadQuery.is_proper_pair: SamFlag.is_proper_pair,
    }

    def __init__(self, bam_file: pathlib.Path, number_of_threads: int = 1):
        self._bam_file = bam_file
        self._number_of_threads = number_of_threads

    def get_number_of_segments(
        self, exclude: Optional[List[repositories.ReadQuery]] = None
    ) -> int:
        if exclude is None:
            exclude = self.default_filters
        samflags = [self.lookup[read_query] for read_query in exclude]
        exclude_flag = sum(samflags)

        count: str = pysam.view(  # type: ignore
            "-@",
            f"{self._number_of_threads}",
            "--count",
            "--exclude-flag",
            f"{exclude_flag}",
            str(self._bam_file),
        )
        return int(count.strip())


class PysamMateFinder:
    """Search for mate reads"""

    def __init__(self, bam_file: pysam.AlignmentFile):
        self._bam_file = bam_file

    def get_mate(self, read: pysam.AlignedSegment) -> Optional[pysam.AlignedSegment]:
        if read.is_supplementary:
            result = None
        else:
            try:
                result = self._bam_file.mate(read)
            except ValueError as exc:
                result = self._get_unmapped_mate(read)
        return result

    def _get_unmapped_mate(
        self, read: pysam.AlignedSegment
    ) -> Optional[pysam.AlignedSegment]:
        possible_mates = list(self._get_read_from_mate_position(read))
        if len(possible_mates) == 1:
            return possible_mates[0]
        else:
            # this is unexpected, but we silently carry on
            return None

    def _get_read_from_mate_position(
        self, read: pysam.AlignedSegment
    ) -> Iterator[pysam.AlignedSegment]:
        mate_start = read.next_reference_start
        mate_end = read.next_reference_start + 1
        mate_name = read.next_reference_name
        for candidate in self._bam_file.fetch(
            start=mate_start, stop=mate_end, contig=mate_name, until_eof=True
        ):
            if self._is_mate(read, candidate):
                yield candidate

    def _is_mate(
        self, given: pysam.AlignedSegment, candidate: pysam.AlignedSegment
    ) -> bool:
        result = (
            (candidate.query_name == given.query_name)
            and (given.is_read1 != candidate.is_read1)
            and not candidate.is_supplementary
        )
        return result


class FilePlacedSegmentRepository(repositories.PlacedSegmentRepository):
    _read_map = {True: entities.ReadOrder.ONE, False: entities.ReadOrder.TWO}

    def __init__(
        self,
        bam_file: pathlib.Path,
        read_counter: repositories.SegmentCounter,
    ):
        self._bam_file = pysam.AlignmentFile(str(bam_file), "rb")
        self._mate_finder = PysamMateFinder(self._bam_file)
        self._read_counter = read_counter

    def get_number_of_segments(
        self, exclude: Optional[List[repositories.ReadQuery]] = None
    ) -> int:
        """Count total number of reads, excluding given list of `ReadQuery`"""
        return self._read_counter.get_number_of_segments(exclude)

    def get(
        self,
        locations: FrozenSet[entities.GenomicRegion],
        exclude: List[repositories.ReadQuery],
        min_quality: int = 0,
    ) -> Iterable[entities.PlacedSegment]:
        exclude_selector = AnySegmentSelector(exclude)
        quality_selector = MapQSelector(minimum_quality=min_quality)
        result = (
            self._convert_read(read, locations)
            for read in self._get_reads_from(locations)
            if not exclude_selector(read) and quality_selector(read)
        )
        return result

    def get_mate_exact_position(self, read: entities.PlacedSegment) -> Tuple[str, int]:
        """Get contig name and starting position of read mate"""
        if isinstance(read, entities.PlacedSegment):
            mate_name = str(read.content.next_reference_name)
            mate_start = read.content.next_reference_start + 1
            return mate_name, mate_start
        else:
            raise NotImplementedError(
                f"have not implemented methods to handle read of type: {type(read)}"
            )

    def get_mate(
        self, read: entities.PlacedSegment
    ) -> Optional[entities.PlacedSegment]:
        if not isinstance(read.content, pysam.AlignedSegment):
            raise ValueError("invalid read type")
        mate = self._mate_finder.get_mate(read.content)
        if mate is None:
            result = None
        else:
            result = self._convert_read(mate, frozenset())
        return result

    def _get_reads_from(
        self, locations: FrozenSet[entities.GenomicRegion]
    ) -> Iterator[pysam.AlignedSegment]:
        for region in locations:
            for read in self._bam_file.fetch(region.chrom, region.start, region.end):
                yield read

    def _convert_read(
        self, read: pysam.AlignedSegment, locations: FrozenSet[entities.GenomicRegion]
    ) -> entities.PlacedSegment:
        read_order = self._read_map[read.is_read1]
        if read.query_name is None:
            raise ValueError("AlignedSegment has no name")
        result = entities.PlacedSegment(read.query_name, locations, read_order, read)
        return result


class PysamSegmentRepositoryFactory(repositories.SegmentRepositoryFactory):
    def __init__(
        self,
        bam_file: pathlib.Path,
        notification_factory: notifications.NotificationServiceFactory,
        number_of_threads: int = 1,
    ):
        self.notification_factory = notification_factory
        self._bam_file = bam_file
        self._number_of_threads = number_of_threads

    def build_counter(
        self,
        features: FrozenSet[repositories.SegmentRepoFeature],
        total_number_of_reads: Optional[int],
    ) -> repositories.SegmentCounter:
        if repositories.SegmentRepoFeature.BUILTIN in features:
            if total_number_of_reads is None:
                raise caller_factories.MissingArgumentError("total_number_of_reads")
            else:
                return repositories.ProvidedSegmentCounter(total_number_of_reads)
        else:
            if repositories.SegmentRepoFeature.WITH_NOTIFICATION in features:
                with_notification = True
            else:
                with_notification = False

            return self.__build_singleton_counter(self._bam_file, with_notification)

    @functools.lru_cache(maxsize=1024)
    def __build_singleton_counter(
        self, bam_file: pathlib.Path, with_notification: bool
    ) -> repositories.CachedSegmentCounter:
        """Both the BamFileSegmentCounter and the CachedSegmentCounter need to be unique"""
        result: repositories.SegmentCounter
        result = BamFileSegmentCounter(self._bam_file, self._number_of_threads)
        if with_notification:
            notification_service = self.notification_factory.build()
            result = notifications.NotifySegmentCounter(result, notification_service)
        result = repositories.CachedSegmentCounter(result)
        return result

    def build(
        self,
        features: FrozenSet[repositories.SegmentRepoFeature],
        total_number_of_reads: Optional[int],
    ) -> repositories.PlacedSegmentRepository:
        counter = self.build_counter(features, total_number_of_reads)
        result = FilePlacedSegmentRepository(self._bam_file, counter)
        return result
