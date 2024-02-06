import pathlib
from unittest import mock

import pytest

from ilmn.pelops import entities, repositories, request_models, result_models
from ilmn.pelops.callers import caller_factories
from ilmn.pelops.infrastructure import pysam_repositories
from ilmn.pelops.interactors import classify_interactor
from tests import stubs

locationA = frozenset([entities.GenomicRegion("chr1", 1, 100)])
locationB = frozenset([entities.GenomicRegion("chr2", 2, 200)])


class TestClassifyInteractor:
    @pytest.fixture
    def fake_segments(self):
        result = [
            entities.PlacedSegment("foo", locationA, entities.ReadOrder.ONE),
            entities.PlacedSegment("foo", locationB, entities.ReadOrder.TWO),
        ]
        return result

    @pytest.fixture
    def expected_segments(self, fake_segments):
        result = [
            result_models.PlacedSegmentDTO(item.read_name) for item in fake_segments
        ]
        return result

    @pytest.fixture
    def segment_repository(self, fake_segments):
        result = stubs.StubPlacedSegmentRepository()
        for segment in fake_segments:
            result.add_read(segment)
        return result

    @pytest.fixture
    def caller_factory(
        self,
        segment_repository,
        segment_repo_factory,
        notification_factory,
    ):
        region_repo = stubs.SmallRegionRepository()
        region_repo_factory = mock.Mock(spec=repositories.RegionRepositoryFactory)
        region_repo_factory.build = mock.Mock(return_value=region_repo)
        candidate_region_caller_factory = caller_factories.CandidateRegionCallerFactory(
            region_repo_factory, segment_repo_factory
        )
        region_caller_factory = caller_factories.RegionPairCallerFactory(
            candidate_region_caller_factory, region_repo_factory, notification_factory
        )
        segment_repo_factory.build = mock.Mock(return_value=segment_repository)
        read_caller_factory = caller_factories.ReadCallerFactory(
            segment_repo_factory, notification_factory
        )
        result = caller_factories.RearrangementCallerFactory(
            segment_repo_factory, region_caller_factory, read_caller_factory
        )
        return result

    @pytest.fixture
    def interactor(self, caller_factory, segment_repo_factory):
        presenter = mock.Mock()
        result = classify_interactor.ClassifyInteractor(
            presenter, caller_factory, segment_repo_factory
        )
        return result

    def test_get_rearrangement_evidence(
        self, interactor, fake_segments, expected_segments
    ):
        request = request_models.ClassifyRequest(
            features=frozenset([request_models.Feature.PROVIDED_READ_COUNT]),
            total_number_of_reads=2_000_000_000,
        )
        dux4_regionset = result_models.CompoundRegionDTO(
            name="CoreDUX4",
            regions=frozenset([result_models.GenomicRegionDTO("chr1", 1, 100)]),
        )
        igh_regionset = result_models.CompoundRegionDTO(
            name="IGH",
            regions=frozenset([result_models.GenomicRegionDTO("chr2", 2, 200)]),
        )
        extended_dux4_regionset = result_models.CompoundRegionDTO(
            name="ExtendedDUX4",
            regions=frozenset([result_models.GenomicRegionDTO("chr3", 3, 300)]),
        )

        expected = result_models.ClassifyResult(
            reference=result_models.ReferenceGenome.GRCh38,
            unique_mapped_reads=2_000_000_000,
            rearrangements=[
                result_models.RearrangementDTO(
                    A=dux4_regionset,
                    B=igh_regionset,
                    evidence=result_models.ReadsEvidence(paired=1, split=0, SRPB=0.5),
                    supporting_reads=set(expected_segments),
                ),
                result_models.RearrangementDTO(
                    A=extended_dux4_regionset,
                    B=igh_regionset,
                    evidence=result_models.ReadsEvidence(paired=0, split=0, SRPB=0.0),
                    supporting_reads=set(),
                ),
            ],
        )
        observed = interactor.get_rearrangement_evidence(request)
        assert observed == expected
