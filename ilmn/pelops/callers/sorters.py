from typing import Iterable

from ilmn.pelops import entities
from ilmn.pelops.callers import rearrangement_callers


class SrpbRearrangementSorter(rearrangement_callers.RearrangementSorter):
    def __call__(
        self, rearrangements: Iterable[entities.Rearrangement]
    ) -> Iterable[entities.Rearrangement]:
        """Sort Rearrangements on decreasing SRPB"""
        result = sorted(rearrangements, key=lambda x: x.srpb, reverse=True)
        return iter(result)
