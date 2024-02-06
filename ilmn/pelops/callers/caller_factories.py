import abc
import enum
from typing import FrozenSet, Optional

from ilmn.pelops import entities, notifications, repositories, selectors
from ilmn.pelops.callers import (
    read_callers,
    rearrangement_callers,
    region_callers,
    region_pair_callers,
    sorters,
)


class MissingArgumentError(ValueError):
    def __init__(self, argument: str):
        message = f"Argument `{argument}` must be provided"
        super().__init__(message)


class CallerType(enum.Enum):
    NAMED = enum.auto()
    BAIT = enum.auto()
    SELECTABLE = enum.auto()
    MULTI = enum.auto()
    SRPB_SORTED = enum.auto()


class CallerFeature(enum.Enum):
    WITH_PROVIDED_READ_COUNT = enum.auto()
    WITH_BLACKLIST = enum.auto()
    WITH_NOTIFICATIONS = enum.auto()


class CandidateRegionCallerFactory:
    def __init__(
        self,
        region_repo_factory: repositories.RegionRepositoryFactory,
        segment_repo_factory: repositories.SegmentRepositoryFactory,
    ):
        self.region_repo_factory = region_repo_factory
        self.segment_repo_factory = segment_repo_factory

    def build(
        self,
        caller_features: FrozenSet[CallerFeature],
        total_number_of_reads: Optional[int] = None,
    ) -> region_callers.CandidateRegionCaller:
        region_repo_type = self.__get_region_repo_type(caller_features)
        selector_features = self.__get_selector_feature(caller_features)
        segment_repo_features = self.__get_segment_repo_features(caller_features)

        region_repo = self.region_repo_factory.build(region_repo_type)

        # TODO: inject this factory, and initialise with another factory, not the region_repository
        selector_factory = selectors.SelectorFactory(region_repo)
        region_selectors = selector_factory.build(features=selector_features)
        segment_repo = self.segment_repo_factory.build(
            segment_repo_features, total_number_of_reads
        )
        result = region_callers.CandidateRegionCaller(segment_repo, region_selectors)
        return result

    def __get_segment_repo_features(
        self, caller_features: FrozenSet[CallerFeature]
    ) -> FrozenSet[repositories.SegmentRepoFeature]:
        result = []
        if CallerFeature.WITH_NOTIFICATIONS in caller_features:
            result.append(repositories.SegmentRepoFeature.WITH_NOTIFICATION)
        if CallerFeature.WITH_PROVIDED_READ_COUNT in caller_features:
            result.append(repositories.SegmentRepoFeature.BUILTIN)
        return frozenset(result)

    def __get_selector_feature(
        self, caller_features: FrozenSet[CallerFeature]
    ) -> FrozenSet[selectors.SelectorFeature]:
        if CallerFeature.WITH_BLACKLIST in caller_features:
            result = frozenset([selectors.SelectorFeature.WITH_BLACKLIST])
        else:
            result = frozenset()
        return result

    def __get_region_repo_type(
        self, caller_features: FrozenSet[CallerFeature]
    ) -> repositories.RegionRepoType:
        if CallerFeature.WITH_BLACKLIST in caller_features:
            result = repositories.RegionRepoType.WITH_BLACKLIST
        else:
            result = repositories.RegionRepoType.BUILTIN
        return result


class ReadCallerFactory:
    __converter = {
        CallerFeature.WITH_NOTIFICATIONS: repositories.SegmentRepoFeature.WITH_NOTIFICATION,
        CallerFeature.WITH_PROVIDED_READ_COUNT: repositories.SegmentRepoFeature.BUILTIN,
    }
    reads_to_exclude = [
        repositories.ReadQuery.is_duplicate,
        repositories.ReadQuery.is_not_paired,
        repositories.ReadQuery.is_qcfail,
    ]

    def __init__(
        self,
        segment_repo_factory: repositories.SegmentRepositoryFactory,
        notification_service_factory: notifications.NotificationServiceFactory,
    ):
        self.segment_repo_factory = segment_repo_factory
        self.notification_service_factory = notification_service_factory

    def build(
        self,
        minimum_mapping_quality: int,
        features: FrozenSet[CallerFeature],
        total_number_of_reads: Optional[int],
    ) -> read_callers.ReadsCaller:
        result: read_callers.ReadsCaller
        repo_features = [
            self.__converter[item] for item in features if item in self.__converter
        ]
        segment_repo = self.segment_repo_factory.build(
            frozenset(repo_features), total_number_of_reads
        )

        result = read_callers.SpanningReadsCaller(
            segment_repo,
            reads_to_exclude=self.reads_to_exclude,
            minimum_mapping_quality=minimum_mapping_quality,
        )
        if CallerFeature.WITH_NOTIFICATIONS in features:
            notification_service = self.notification_service_factory.build()
            result = notifications.NotifyReadsCaller(result, notification_service)
        return result


class RegionPairCallerFactory:
    def __init__(
        self,
        candidate_region_caller_factory: CandidateRegionCallerFactory,
        region_repo_factory: repositories.RegionRepositoryFactory,
        notification_factory: notifications.NotificationServiceFactory,
    ):
        self.candidate_region_caller_factory = candidate_region_caller_factory
        self.region_repo_factory = region_repo_factory
        self.notification_factory = notification_factory

    def build(
        self,
        caller_type: CallerType,
        caller_features: FrozenSet[CallerFeature],
        total_number_of_reads: Optional[int],
    ) -> region_pair_callers.CompoundRegionPairCaller:
        result: region_pair_callers.CompoundRegionPairCaller
        region_repo = self.__get_region_repo(caller_features)

        if caller_type == CallerType.BAIT:
            result = self.__build_bait_reagion_pair_caller(
                region_repo, caller_features, total_number_of_reads
            )
        elif caller_type == CallerType.NAMED:
            result = region_pair_callers.NamedRegionPairCaller(region_repo)
        else:
            raise ValueError()
        if CallerFeature.WITH_NOTIFICATIONS in caller_features:
            notification_service = self.notification_factory.build()
            result = notifications.NotifyRegionPairCaller(result, notification_service)
        return result

    def __get_region_repo(
        self, caller_features: FrozenSet[CallerFeature]
    ) -> repositories.RegionRepository:
        if CallerFeature.WITH_BLACKLIST in caller_features:
            region_repo_type = repositories.RegionRepoType.WITH_BLACKLIST
        else:
            region_repo_type = repositories.RegionRepoType.BUILTIN
        return self.region_repo_factory.build(region_repo_type)

    def __build_bait_reagion_pair_caller(
        self,
        region_repo: repositories.RegionRepository,
        caller_features: FrozenSet[CallerFeature],
        total_number_of_reads: Optional[int],
    ) -> region_pair_callers.BaitRegionPairCaller:
        bait = entities.RegionsName.CoreDUX4
        region_caller = self.candidate_region_caller_factory.build(
            caller_features, total_number_of_reads
        )
        result = region_pair_callers.BaitRegionPairCaller(
            region_repo, region_caller, bait
        )
        return result


class RearrangementCallerFactory:
    def __init__(
        self,
        repo_factory: repositories.SegmentRepositoryFactory,
        region_caller_factory: RegionPairCallerFactory,
        read_caller_factory: ReadCallerFactory,
    ):
        self._repo_factory = repo_factory
        self._region_pair_caller_factory = region_caller_factory
        self._read_caller_factory = read_caller_factory

    def build(
        self,
        caller_type: CallerType,
        features: FrozenSet[CallerFeature],
        minimum_mapping_quality: Optional[int] = None,
        srpb_threshold: Optional[float] = None,
        total_number_of_reads: Optional[int] = None,
    ) -> rearrangement_callers.RearrangementCaller:
        self._srpb_threshold = srpb_threshold
        self._minimum_mapping_quality = minimum_mapping_quality
        return self._build(caller_type, features, total_number_of_reads)

    def _build(
        self,
        caller_type: CallerType,
        caller_features: FrozenSet[CallerFeature],
        total_number_of_reads: Optional[int],
    ) -> rearrangement_callers.RearrangementCaller:
        if caller_type == CallerType.MULTI:
            named_caller = self._build(
                CallerType.NAMED, caller_features, total_number_of_reads
            )
            selectable_caller = self._build(
                CallerType.SRPB_SORTED, caller_features, total_number_of_reads
            )
            callers = [named_caller, selectable_caller]
            return rearrangement_callers.MultiRearrangementCaller(callers)

        elif caller_type == CallerType.SRPB_SORTED:
            sorter = sorters.SrpbRearrangementSorter()
            caller = self._build(
                CallerType.SELECTABLE, caller_features, total_number_of_reads
            )
            return rearrangement_callers.SortedRearrangementCaller(caller, sorter)

        elif caller_type == CallerType.SELECTABLE:
            caller = self._build(
                CallerType.BAIT, caller_features, total_number_of_reads
            )
            if self._srpb_threshold is None:
                raise MissingArgumentError("srpb_threshold")
            selector = selectors.SRPBRearrangementSelector(self._srpb_threshold)
            return rearrangement_callers.SelectableRearrangementCaller(caller, selector)
        else:
            minimum_mapping_quality = self.__get_minimum_mapping_quality(caller_type)
            reads_caller = self._read_caller_factory.build(
                minimum_mapping_quality, caller_features, total_number_of_reads
            )
            region_pair_caller = self._region_pair_caller_factory.build(
                caller_type, caller_features, total_number_of_reads
            )

            segment_repo_features = self.__get_segment_repo_features(caller_features)
            reads_counter = self._repo_factory.build_counter(
                segment_repo_features, total_number_of_reads
            )
            return rearrangement_callers.RegionRearrangementCaller(
                reads_caller, region_pair_caller, reads_counter
            )

    def __get_segment_repo_features(
        self, features: FrozenSet[CallerFeature]
    ) -> FrozenSet[repositories.SegmentRepoFeature]:
        result = []
        if CallerFeature.WITH_PROVIDED_READ_COUNT in features:
            result.append(repositories.SegmentRepoFeature.BUILTIN)
        if CallerFeature.WITH_NOTIFICATIONS in features:
            result.append(repositories.SegmentRepoFeature.WITH_NOTIFICATION)
        return frozenset(result)

    def __get_minimum_mapping_quality(self, caller_type: CallerType) -> int:
        if caller_type == CallerType.NAMED:
            minimum_mapping_quality = 0
        elif caller_type == CallerType.BAIT:
            if self._minimum_mapping_quality is None:
                raise MissingArgumentError("minimum_mapping_quality")
            else:
                minimum_mapping_quality = self._minimum_mapping_quality
        else:
            raise ValueError(f"Unable to build caller for {caller_type.name}")
        return minimum_mapping_quality
