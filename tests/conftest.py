"""Shared fixtures"""

import pathlib
from unittest import mock

import pytest

from ilmn.pelops import entities, notifications, repositories
from ilmn.pelops.callers import caller_factories
from ilmn.pelops.infrastructure import blacklist_region_repository, pysam_repositories


@pytest.fixture
def exclude_flags():
    result = [
        repositories.ReadQuery.is_unmapped,
        repositories.ReadQuery.is_secondary,
        repositories.ReadQuery.is_qcfail,
        repositories.ReadQuery.is_duplicate,
        repositories.ReadQuery.is_supplementary,
    ]
    return result


@pytest.fixture
def test_folder():
    return pathlib.Path(__file__).parent.absolute()


@pytest.fixture
def bedfile(test_folder):
    return test_folder / "data" / "blacklist.bed"


@pytest.fixture
def bam_file(test_folder):
    return test_folder / "data" / "HCC1187BL_IGH_DUX4.bam"


@pytest.fixture(params=["bam", "cram"])
def alignment_file(test_folder, request):
    extension = request.param
    return test_folder / "data" / ("HCC1187BL_IGH_DUX4." + extension)


@pytest.fixture
def an_unnamed_genomic_region_set():
    region = entities.GenomicRegion("chr1", 1, 200)
    result = entities.CompoundRegion(entities.RegionsName.UNNAMED, frozenset([region]))
    return result


#### FACTORIES ####
# here we create factories that can be used to tests objects as builds by them


@pytest.fixture
def notification_factory():
    result = notifications.SimpleNotificationServiceFactory()
    return result


@pytest.fixture
def segment_repo_factory(bam_file, notification_factory):
    result = pysam_repositories.PysamSegmentRepositoryFactory(
        bam_file, notification_factory
    )
    return result


@pytest.fixture
def region_repo_factory(bedfile):
    result = blacklist_region_repository.FileRegionRepositoryFactory(bedfile)
    return result


@pytest.fixture
def read_caller_factory(segment_repo_factory, notification_factory):
    result = caller_factories.ReadCallerFactory(
        segment_repo_factory, notification_factory
    )
    return result


@pytest.fixture
def candidate_region_caller_factory(segment_repo_factory, region_repo_factory):
    result = caller_factories.CandidateRegionCallerFactory(
        region_repo_factory, segment_repo_factory
    )
    return result


@pytest.fixture
def region_pair_caller_factory(
    candidate_region_caller_factory, region_repo_factory, notification_factory
):
    result = caller_factories.RegionPairCallerFactory(
        candidate_region_caller_factory, region_repo_factory, notification_factory
    )
    return result


@pytest.fixture
def rearrangemet_caller_factory(
    segment_repo_factory, read_caller_factory, region_pair_caller_factory
):
    result = caller_factories.RearrangementCallerFactory(
        segment_repo_factory, region_pair_caller_factory, read_caller_factory
    )
    return result
