"""Interactors covering all usecases for pelops"""

import abc
from typing import FrozenSet, List

from ilmn.pelops import entities, repositories, request_models, result_models
from ilmn.pelops.callers import caller_factories


class ClassifyPresenter(abc.ABC):
    @abc.abstractmethod
    def present_classification(self, result: result_models.ClassifyResult) -> None:
        """Present the result from classify_reads"""


def convert_compound_region(
    domain: entities.CompoundRegion,
) -> result_models.CompoundRegionDTO:
    """Convert the domain object `CompoundRegion` into a data structure"""
    regions = [convert_genomic_region(item) for item in domain.regions]
    result = result_models.CompoundRegionDTO(
        name=domain.name.name, regions=frozenset(regions)
    )
    return result


def convert_genomic_region(
    domain: entities.GenomicRegion,
) -> result_models.GenomicRegionDTO:
    """Convert the domain object `GenomicRegionDTO` into a data structure"""
    result = result_models.GenomicRegionDTO(
        chrom=domain.chrom, start=domain.start, end=domain.end
    )
    return result


def convert_segment(domain: entities.PlacedSegment) -> result_models.PlacedSegmentDTO:
    result = result_models.PlacedSegmentDTO(
        read_name=domain.read_name, content=domain.content
    )
    return result


def convert_rearrangement(
    rearrangement: entities.Rearrangement,
) -> result_models.RearrangementDTO:
    evidence = result_models.ReadsEvidence(
        rearrangement.counts.paired, rearrangement.counts.split, rearrangement.srpb
    )
    segments = [convert_segment(item) for item in rearrangement.segments]
    a, b = rearrangement.region_pair
    result = result_models.RearrangementDTO(
        convert_compound_region(a),
        convert_compound_region(b),
        evidence,
        set(segments),
    )
    return result


class ClassifyInteractor:
    def __init__(
        self,
        presenter: ClassifyPresenter,
        rearrangement_caller_factory: caller_factories.RearrangementCallerFactory,
        repository_factory: repositories.SegmentRepositoryFactory,
    ):
        self._presenter = presenter
        self._caller_factory = rearrangement_caller_factory
        self._repo_factory = repository_factory

    def get_rearrangement_evidence(
        self, request: request_models.ClassifyRequest
    ) -> result_models.ClassifyResult:
        """Classify reads spanning several GenomiRegionSets pairs into
        spanning, paired, and split reads."""
        rearrangement_caller = self._caller_factory.build(
            self.__get_caller_type(request.features),
            self.__get_caller_features(request.features),
            request.minimum_mapping_quality,
            request.srpb_threshold,
            request.total_number_of_reads,
        )

        reads_counter = self._repo_factory.build_counter(
            self.__get_counter_features(request.features),
            request.total_number_of_reads,
        )

        rearrangements = [
            convert_rearrangement(item)
            for item in rearrangement_caller.get_rearrangements()
        ]
        result = result_models.ClassifyResult(
            reference=result_models.ReferenceGenome.GRCh38,
            unique_mapped_reads=reads_counter.get_number_of_segments(),
            rearrangements=list(rearrangements),
        )
        return result

    def present_rearrangement_evidence(
        self, request: request_models.ClassifyRequest
    ) -> None:
        result = self.get_rearrangement_evidence(request)
        self._presenter.present_classification(result)

    def __get_caller_type(
        self, features: FrozenSet[request_models.Feature]
    ) -> caller_factories.CallerType:
        if request_models.Feature.DUX4_OTHER in features:
            return caller_factories.CallerType.MULTI
        else:
            return caller_factories.CallerType.NAMED

    def __get_caller_features(
        self, features: FrozenSet[request_models.Feature]
    ) -> FrozenSet[caller_factories.CallerFeature]:
        converter = {
            request_models.Feature.PROVIDED_READ_COUNT: caller_factories.CallerFeature.WITH_PROVIDED_READ_COUNT,
            request_models.Feature.WITH_BLACKLIST: caller_factories.CallerFeature.WITH_BLACKLIST,
            request_models.Feature.WITH_NOTIFICATIONS: caller_factories.CallerFeature.WITH_NOTIFICATIONS,
        }
        result = [converter[feature] for feature in features if feature in converter]
        return frozenset(result)

    def __get_counter_features(
        self, features: FrozenSet[request_models.Feature]
    ) -> FrozenSet[repositories.SegmentRepoFeature]:
        converter = {
            request_models.Feature.WITH_NOTIFICATIONS: repositories.SegmentRepoFeature.WITH_NOTIFICATION,
            request_models.Feature.PROVIDED_READ_COUNT: repositories.SegmentRepoFeature.BUILTIN,
        }
        result = [converter[feature] for feature in features if feature in converter]
        return frozenset(result)
