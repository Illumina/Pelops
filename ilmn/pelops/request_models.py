import dataclasses
import enum
from typing import FrozenSet, Optional


class Feature(enum.Enum):
    DUX4_OTHER = enum.auto()
    PROVIDED_READ_COUNT = enum.auto()
    WITH_BLACKLIST = enum.auto()
    WITH_NOTIFICATIONS = enum.auto()


@dataclasses.dataclass
class ClassifyRequest:
    features: FrozenSet[Feature] = frozenset()
    minimum_mapping_quality: Optional[int] = None
    srpb_threshold: Optional[float] = None
    total_number_of_reads: Optional[int] = None
