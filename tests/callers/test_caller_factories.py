import pytest

from ilmn.pelops import notifications, repositories
from ilmn.pelops.callers import (
    caller_factories,
    read_callers,
    rearrangement_callers,
    region_callers,
    region_pair_callers,
)
from ilmn.pelops.infrastructure import blacklist_region_repository, pysam_repositories

# defined in conftest
# - fixture bedfile
# - notification_factory
# - region_repo_factory
# - candidate_region_caller_factory


class TestRegionPairCallerFactory:
    TEST_CASES = [
        pytest.param(
            caller_factories.CallerType.NAMED,
            frozenset(),
            region_pair_callers.NamedRegionPairCaller,
            id="named",
        ),
        pytest.param(
            caller_factories.CallerType.NAMED,
            frozenset([caller_factories.CallerFeature.WITH_NOTIFICATIONS]),
            notifications.NotifyRegionPairCaller,
            id="named, notifications",
        ),
        pytest.param(
            caller_factories.CallerType.BAIT,
            frozenset(),
            region_pair_callers.BaitRegionPairCaller,
            id="bait",
        ),
        pytest.param(
            caller_factories.CallerType.BAIT,
            frozenset([caller_factories.CallerFeature.WITH_NOTIFICATIONS]),
            notifications.NotifyRegionPairCaller,
            id="bait, notifications",
        ),
    ]

    @pytest.mark.parametrize("caller_type, features, expected", TEST_CASES)
    def test_build(self, region_pair_caller_factory, caller_type, features, expected):
        observed = region_pair_caller_factory.build(
            caller_type=caller_type,
            caller_features=features,
            total_number_of_reads=None,
        )
        assert isinstance(observed, expected)


class TestReadCallerFactory:
    testcases = [
        pytest.param(
            frozenset(), read_callers.SpanningReadsCaller, id="no_notifications"
        ),
        pytest.param(
            frozenset([caller_factories.CallerFeature.WITH_NOTIFICATIONS]),
            notifications.NotifyReadsCaller,
            id="with_notifications",
        ),
    ]

    @pytest.mark.parametrize("features, expected", testcases)
    def test_build(self, features, expected, read_caller_factory):
        observed = read_caller_factory.build(
            minimum_mapping_quality=1, features=features, total_number_of_reads=None
        )
        assert isinstance(observed, expected)


class TestRearrangementCallerFactory:
    def test_build(self, rearrangemet_caller_factory):
        caller_type = caller_factories.CallerType.NAMED
        caller = rearrangemet_caller_factory.build(caller_type, features=frozenset())
        assert isinstance(caller, rearrangement_callers.RegionRearrangementCaller)

    def test_build_fails(self, rearrangemet_caller_factory):
        caller_type = caller_factories.CallerType.SELECTABLE
        with pytest.raises(caller_factories.MissingArgumentError):
            rearrangemet_caller_factory.build(caller_type, features=frozenset())


class TestCandidateRegionCallerFactory:
    testcases = [
        frozenset(),
        frozenset([caller_factories.CallerFeature.WITH_BLACKLIST]),
    ]

    @pytest.mark.parametrize("features", testcases)
    def test_build(self, candidate_region_caller_factory, features):
        caller = candidate_region_caller_factory.build(caller_features=features)
        assert isinstance(caller, region_callers.CandidateRegionCaller)
