"""Result Model for classify"""

import dataclasses
import enum
from typing import Any, FrozenSet, List, Optional, Set, Tuple

# NOTE: order of class definition in the file matters. Early class definition
# should not contain classes defined below itself, so we can refer to them by
# class name rather than by string


class ReferenceGenome(enum.Enum):
    GRCh38 = enum.auto()


@dataclasses.dataclass(frozen=True)
class GenomicRegionDTO:
    chrom: str
    start: int
    end: int


@dataclasses.dataclass(frozen=True)
class CompoundRegionDTO:
    name: str
    regions: FrozenSet[GenomicRegionDTO]


@dataclasses.dataclass(frozen=True)
class ReadsEvidence:
    paired: int
    split: int
    SRPB: float  # spanning read pairs per billion


@dataclasses.dataclass(frozen=True)
class PlacedSegmentDTO:
    read_name: str
    content: Any = None


@dataclasses.dataclass(frozen=True)
class RearrangementDTO:
    A: CompoundRegionDTO
    B: CompoundRegionDTO
    evidence: ReadsEvidence
    supporting_reads: Set[PlacedSegmentDTO]


@dataclasses.dataclass
class ClassifyResult:
    reference: ReferenceGenome
    unique_mapped_reads: int
    rearrangements: List[RearrangementDTO]
