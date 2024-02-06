"""Module dealing with compound regions handled in pairs."""

import abc
from typing import Iterable

from ilmn.pelops import entities
from ilmn.pelops.callers import region_callers
from ilmn.pelops.repositories import RegionRepository


class CompoundRegionPairCaller(abc.ABC):
    """Interface for CompoundRegionPairCaller.

    Classes in this hierarchy identify pairs of compound regions according to
    various criteria.

    """

    @abc.abstractmethod
    def get_compound_region_pairs(self) -> Iterable[entities.CompoundRegionPair]:
        """Get Pairs of CompoundRegion"""


class NamedRegionPairCaller(CompoundRegionPairCaller):
    """Find CompoundRegion pairs"""

    name_pairs = [
        (entities.RegionsName.CoreDUX4, entities.RegionsName.IGH),
        (entities.RegionsName.ExtendedDUX4, entities.RegionsName.IGH),
    ]

    def __init__(self, region_repository: RegionRepository):
        self._region_repo = region_repository

    def get_compound_region_pairs(self) -> Iterable[entities.CompoundRegionPair]:
        """Convert tuples of `RegionName` to tuples of `CompoundRegion`"""
        for pair in self.name_pairs:
            result = entities.CompoundRegionPair(
                self._region_repo.get(pair[0]),
                self._region_repo.get(pair[1]),
            )
            yield result


class BaitRegionPairCaller(CompoundRegionPairCaller):
    """Find CompoundRegion pairs by looking for candidates to a specific
    "bait" region."""

    def __init__(
        self,
        region_repository: RegionRepository,
        region_caller: region_callers.CandidateRegionCaller,
        bait: entities.RegionsName,
    ):
        self._region_repo = region_repository
        self._region_caller = region_caller
        self._bait = bait

    def get_compound_region_pairs(self) -> Iterable[entities.CompoundRegionPair]:
        bait_compound_region = self._region_repo.get(self._bait)
        for candidate_region in self._region_caller.get_candidates_regions(
            bait_compound_region
        ):
            candidate_compound_region = entities.CompoundRegion(
                entities.RegionsName.UNNAMED, frozenset([candidate_region])
            )

            result = entities.CompoundRegionPair(
                bait_compound_region, candidate_compound_region
            )

            yield result
