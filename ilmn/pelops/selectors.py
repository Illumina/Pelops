import abc
import enum
from typing import FrozenSet, List, Optional

from ilmn.pelops import entities, repositories


class RearrangementSelector(abc.ABC):
    """Determine if a RearrangementDTO is to be selected or not"""

    @abc.abstractmethod
    def __call__(self, rearrangement: entities.Rearrangement) -> bool:
        """True if this rearrangement is to be selected/kept"""
        pass


class RegionSelector(abc.ABC):
    """Determine if a GenomicRegion is to be selected or not"""

    @abc.abstractmethod
    def __call__(self, region: entities.GenomicRegion) -> bool:
        """True if this region needs to be selected/kept"""


class NamedRegionSelector(RegionSelector):
    """Select a GenomicRegion if it overlaps a named set of GenomicRegions"""

    def __init__(
        self,
        region_name: entities.RegionsName,
        region_repository: repositories.RegionRepository,
    ):
        self._region_name = region_name
        self._repo = region_repository

        # cache regions, so we use the repository only once
        self._regionset: Optional[entities.CompoundRegion] = None

    def __call__(self, region: entities.GenomicRegion) -> bool:
        if not self._regionset:
            self._regionset = self._repo.get(self._region_name)
        result = any([region.overlaps(target) for target in self._regionset.regions])
        return result


class RegionDeSelector(RegionSelector):
    """Inverts the behaviour of a RegionSelector"""

    def __init__(self, selector: RegionSelector):
        self._selector = selector

    def __call__(self, region: entities.GenomicRegion) -> bool:
        result = not self._selector(region)
        return result


class SRPBRearrangementSelector(RearrangementSelector):
    """Determine if a RearrangementDTO is to be selected or not based on SRPB"""

    def __init__(self, threshold: float):
        self._threshold = threshold

    def __call__(self, rearrangement: entities.Rearrangement) -> bool:
        return rearrangement.srpb >= self._threshold


class SelectorFeature(enum.Enum):
    WITH_BLACKLIST = enum.auto()


class SelectorFactory:
    def __init__(self, regions_repository: repositories.RegionRepository):
        self._regions_repository = regions_repository

    def build(
        self, features: FrozenSet[SelectorFeature] = frozenset()
    ) -> List[RegionSelector]:
        result = [
            NamedRegionSelector(entities.RegionsName.GRCh38, self._regions_repository),
            RegionDeSelector(
                NamedRegionSelector(
                    entities.RegionsName.ExtendedDUX4, self._regions_repository
                )
            ),
            RegionDeSelector(
                NamedRegionSelector(entities.RegionsName.IGH, self._regions_repository)
            ),
        ]
        if SelectorFeature.WITH_BLACKLIST in features:
            result.append(
                RegionDeSelector(
                    NamedRegionSelector(
                        entities.RegionsName.BLACKLIST, self._regions_repository
                    )
                )
            )

        return result
