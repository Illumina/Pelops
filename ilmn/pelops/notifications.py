import abc
import functools
from typing import FrozenSet, Iterable, List, Optional

from ilmn.pelops import entities, repositories
from ilmn.pelops.callers import read_callers, region_pair_callers


class NotificationService(abc.ABC):
    """A service to dipatch notification"""

    @abc.abstractmethod
    def notify(self, message: str) -> None:
        """Dispatch one notification"""


class NotificationServiceFactory(abc.ABC):
    @abc.abstractmethod
    def build(self) -> NotificationService:
        """Build notification service"""


class SimpleNotificationService(NotificationService):
    """A simple notification service that uses native `print` to notify"""

    def __init__(self, silent: bool = True) -> None:
        self.__store: List[str] = []
        self.__silent = silent

    def notify(self, message: str) -> None:
        if not self.__silent:
            print(message)
        self.__store.append(message)

    def get_notifications(self) -> List[str]:
        """Retrieve all the notification received so far"""
        return self.__store


class SimpleNotificationServiceFactory(NotificationServiceFactory):
    def __init__(self, silent: bool = True):
        self.__silent = silent

    @functools.lru_cache(maxsize=1)
    def build(self) -> NotificationService:
        return SimpleNotificationService(self.__silent)


class NotifySegmentCounter(repositories.SegmentCounter):
    """SegmentCounter composed with NotificationService"""

    def __init__(
        self,
        read_counter: repositories.SegmentCounter,
        notification_service: NotificationService,
    ):
        self.__read_counter = read_counter
        self.__notification_service = notification_service

    def get_number_of_segments(
        self, exclude: Optional[List[repositories.ReadQuery]] = None
    ) -> int:
        message = "Counting number of unique and mapped reads."
        self.__notification_service.notify(message)
        return self.__read_counter.get_number_of_segments(exclude)


class NotifyReadsCaller(read_callers.ReadsCaller):
    def __init__(
        self,
        read_caller: read_callers.ReadsCaller,
        notification_service: NotificationService,
    ):
        self.read_caller = read_caller
        self.notification_service = notification_service

    def detect_reads_spanning_regions(
        self,
        a_regions: FrozenSet[entities.GenomicRegion],
        b_regions: FrozenSet[entities.GenomicRegion],
    ) -> None:
        self.notification_service.notify("Finding evidence for rearrangement.")
        self.read_caller.detect_reads_spanning_regions(a_regions, b_regions)

    def get_segment_count(self) -> entities.ClassifiedSegmentCount:
        return self.read_caller.get_segment_count()

    def get_segments_of_spanning_reads(self) -> Iterable[entities.PlacedSegment]:
        return self.read_caller.get_segments_of_spanning_reads()


class NotifyRegionPairCaller(region_pair_callers.CompoundRegionPairCaller):
    def __init__(
        self,
        region_pair_caller: region_pair_callers.CompoundRegionPairCaller,
        notification_service: NotificationService,
    ):
        self.__region_pair_caller = region_pair_caller
        self.__notification_service = notification_service

    def get_compound_region_pairs(self) -> Iterable[entities.CompoundRegionPair]:
        """Report about compund region pairs retrived"""
        for region_pair in self.__region_pair_caller.get_compound_region_pairs():
            A, B = (item.name for item in region_pair.get_names())
            if B == entities.RegionsName.UNNAMED.name:
                B = "candidate-region"
            message = f"Retrieving region coordinates for {A}/{B}."
            self.__notification_service.notify(message)
            yield region_pair
