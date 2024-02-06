import functools
import pathlib
from typing import Iterable, List, Optional

from ilmn.pelops import entities, repositories


class InvalidRegionError(ValueError):
    def __init__(self, region: entities.GenomicRegion) -> None:
        message = f"invalid region: {region.chrom}\t{region.start}-{region.end}"
        super().__init__(message)


class BlackListRegionRepository(repositories.RegionRepository):
    def __init__(
        self, file: pathlib.Path, performs_check_at_init: bool = False
    ) -> None:
        self.__file = file
        self.__builtin_repo = repositories.BuiltinRegionRepository()
        if performs_check_at_init:
            self.get(entities.RegionsName.BLACKLIST)

    @functools.lru_cache(maxsize=None)
    def get(self, name: entities.RegionsName) -> entities.CompoundRegion:
        if name != entities.RegionsName.BLACKLIST:
            return self.__builtin_repo.get(name)
        else:
            regions = self.__extract_regions()
            return entities.CompoundRegion(name=name, regions=frozenset(regions))

    def __extract_regions(self) -> Iterable[entities.GenomicRegion]:
        with open(self.__file, "r") as fh:
            for line in fh:
                if not (line.startswith("#") or line.startswith("track")):
                    columns = self.__split_on_whitespaces(line)
                    region = entities.GenomicRegion(
                        columns[0], int(columns[1]), int(columns[2])
                    )
                    self.__raise_if_invalid_order(region)
                    self.__raise_if_invalid_region(region)
                    yield region

    def __split_on_whitespaces(self, line: str) -> List[str]:
        """split a line on regex [ \t]+"""
        parts = line.strip().replace("\t", " ").split()
        return parts

    def __raise_if_invalid_order(self, region: entities.GenomicRegion) -> None:
        if region.end <= region.start:
            raise InvalidRegionError(region)

    def __raise_if_invalid_region(self, region: entities.GenomicRegion) -> None:
        genome = self.__builtin_repo.get(entities.RegionsName.GRCh38)
        if not genome.contains(region):
            raise InvalidRegionError(region)


class FileRegionRepositoryFactory(repositories.RegionRepositoryFactory):
    def __init__(self, blacklist_bed_file: Optional[pathlib.Path] = None):
        self.__bed_file = blacklist_bed_file

    def build(
        self, repo_type: repositories.RegionRepoType
    ) -> repositories.RegionRepository:
        if repo_type == repositories.RegionRepoType.WITH_BLACKLIST:
            if self.__bed_file is None:
                raise ValueError(
                    "To build blacklist repository a bed file must be provided"
                )
            else:
                return BlackListRegionRepository(
                    self.__bed_file, performs_check_at_init=True
                )
        else:
            return repositories.BuiltinRegionRepository()
