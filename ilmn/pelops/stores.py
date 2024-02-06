"""Store data in domain specific objects"""

from typing import Container, Dict, Iterable, Iterator, Set, Tuple, cast

from ilmn.pelops import entities


class UnknowReadPairError(KeyError):
    def __init__(self, read_name: str) -> None:
        message = f"Unknown read pair named {read_name}"
        super().__init__(message)


class UnnamedReadPairError(ValueError):
    def __init__(self) -> None:
        message = f"Unable to store unnamed read pair"
        super().__init__(message)


class PairedReadStore:
    """Initialise, Store and Mark ReadPairs"""

    def __init__(self) -> None:
        self._reads: Dict[str, entities.ReadPair] = {}
        self._paired: Set[str] = set()
        self._spanning: Set[str] = set()
        self._split: Set[str] = set()

    def __iter__(self) -> Iterator[entities.ReadPair]:
        return iter(self._reads.values())

    def mark_as_split(self, read: entities.ReadPair) -> None:
        read_name = cast(str, read.get_name())
        self._split.add(read_name)

    def mark_as_paired(self, read: entities.ReadPair) -> None:
        read_name = cast(str, read.get_name())
        self._paired.add(read_name)

    def mark_as_spanning(self, read: entities.ReadPair) -> None:
        read_name = cast(str, read.get_name())
        self._spanning.add(read_name)

    def get_segment_count(self) -> entities.ClassifiedSegmentCount:
        paired = len(self._paired)
        split = len(self._split)
        spanning = len(self._spanning)
        result = entities.ClassifiedSegmentCount(paired, split, spanning)
        return result

    def get_read_names(self) -> Container[str]:
        return self._reads.keys()

    def get_read_pair(self, name: str) -> entities.ReadPair:
        try:
            result = self._reads[name]
        except KeyError:
            raise UnknowReadPairError(name)
        return result

    def update(self, read_pair: entities.ReadPair) -> None:
        read_name = read_pair.get_name()
        if read_name is None:
            raise UnnamedReadPairError()
        self._reads[read_name] = read_pair

    def get_spanning_reads(self) -> Iterable[entities.ReadPair]:
        return (self._reads[name] for name in self._spanning)

    def clear(self) -> None:
        self._reads.clear()
        self._split.clear()
        self._spanning.clear()
        self._paired.clear()


class SingleChromosomeRegionStore:
    """Store GenomicRegions for a particular chromosome and count repeated
    occurrences"""

    def __init__(self, chromosome_name: str):
        self._chrom_name = chromosome_name
        self._store: Dict[entities.GenomicRegion, int] = {}

    def add(self, region: entities.GenomicRegion) -> None:
        if region.chrom != self._chrom_name:
            raise KeyError(f"Unable to add region from chromosome {region.chrom}")
        elif region not in self._store.keys():
            self._store[region] = 0
        self._store[region] += 1

    def get(self) -> Iterable[Tuple[entities.GenomicRegion, int]]:
        regions = sorted(self._store.keys(), key=lambda region: region.start)
        for region in regions:
            yield region, self._store[region]

    def delete(self, region: entities.GenomicRegion) -> None:
        del self._store[region]


class RegionStore:
    """Store GenomicRegions sorted by chromosome and start position and count
    repeated occurrences"""

    def __init__(self) -> None:
        self._store: Dict[str, "SingleChromosomeRegionStore"] = {}

    def add(self, region: entities.GenomicRegion) -> None:
        if region.chrom not in self._store.keys():
            self._store[region.chrom] = SingleChromosomeRegionStore(region.chrom)
        ministore = self._store[region.chrom]
        ministore.add(region)

    def get_counted_regions(self) -> Iterable[Tuple[entities.GenomicRegion, int]]:
        """Get regions and number of observed occurrences sorted by chromosome
        and start"""

        for chrom in sorted(self._store.keys()):
            ministore = self._store[chrom]
            for region, count in ministore.get():
                yield region, count

    def get_regions(self) -> Iterator[entities.GenomicRegion]:
        """Get regions of observed occurrences sorted by chromosome and start"""
        for region, count in self.get_counted_regions():
            yield region

    def delete(self, region: entities.GenomicRegion) -> None:
        ministore = self._store[region.chrom]
        ministore.delete(region)
