import abc
import dataclasses
from typing import Iterable, Tuple

import pysam

from ilmn.pelops import result_models
from ilmn.pelops.interactors import classify_interactor
from ilmn.pelops.ui.cli.presenters import introspection


@dataclasses.dataclass
class ReadViewModel:
    id: str
    region_names: Tuple[str, str]
    segments: Iterable[pysam.AlignedSegment]
    program_name: str
    program_version: str
    cli_command: str


class ReadsView(abc.ABC):
    @abc.abstractmethod
    def update_reads_model(self, model: Iterable[ReadViewModel]) -> None:
        """Update the reads view model"""


class Formatter:
    """Format number coherently"""

    def format_id(self, id_: int) -> str:
        result = format(id_, "02d")
        return result


class ReadViewModelConverter:
    def __init__(
        self,
        cli_introspection: introspection.CliIntrospection,
    ):
        self._formatter = Formatter()
        self._cli_introspection = cli_introspection

    def convert(
        self, i: int, rearrangement: result_models.RearrangementDTO
    ) -> ReadViewModel:
        id_ = self._formatter.format_id(i)
        region_names = (rearrangement.A.name, rearrangement.B.name)
        sorted_segments = sorted(
            rearrangement.supporting_reads, key=lambda segment: segment.read_name
        )
        sorted_reads = [item.content for item in sorted_segments]

        entrypoint_details = self._cli_introspection.get_entrypoint_details()
        result = ReadViewModel(
            id=id_,
            region_names=region_names,
            segments=sorted_reads,
            program_name=entrypoint_details.program_name,
            program_version=entrypoint_details.program_version,
            cli_command=entrypoint_details.cli_command,
        )
        return result


class ReadsPresenter(classify_interactor.ClassifyPresenter):
    def __init__(
        self,
        reads_view: ReadsView,
        converter: ReadViewModelConverter,
    ):
        self._reads_view = reads_view
        self._converter = converter

    def present_classification(self, result: result_models.ClassifyResult) -> None:
        view_model = [
            self._converter.convert(i + 1, rearrangement)
            for i, rearrangement in enumerate(result.rearrangements)
        ]
        self._reads_view.update_reads_model(view_model)


class NullReadPresenter(classify_interactor.ClassifyPresenter):
    """The Reads Presenter implementation that does not present reads"""

    def __init__(self) -> None:
        pass

    def present_classification(self, result: result_models.ClassifyResult) -> None:
        pass


@dataclasses.dataclass
class ViewEvidence:
    paired_reads: int
    split_reads: int
    SRPB: float  # spanning read pairs per billion


@dataclasses.dataclass
class ViewRegion:
    chrom: str
    start: int
    end: int


@dataclasses.dataclass
class ViewCompoundRegion:
    name: str
    regions: Iterable[ViewRegion]


@dataclasses.dataclass
class ViewRearrangement:
    id: str
    A: ViewCompoundRegion
    B: ViewCompoundRegion
    evidence: ViewEvidence


@dataclasses.dataclass
class ClassificationViewModel:
    reference: str
    unique_mapped_reads: int
    rearrangements: Iterable[ViewRearrangement]
    program_name: str
    version: str
    cli_command: str


class ClassificationView(abc.ABC):
    @abc.abstractmethod
    def update_classification_model(self, result: ClassificationViewModel) -> None:
        """Update the classification view model"""


class ResultModelConverter:
    """Convert part of the result model to part of the view model"""

    def __init__(self) -> None:
        self._formatter = Formatter()

    def convert_result_model(
        self,
        result: result_models.ClassifyResult,
        cli_details: introspection.EntrypointDetails,
    ) -> ClassificationViewModel:
        reference = result.reference.name
        unique_mapped_reads = result.unique_mapped_reads
        rearrangements = [
            self.convert_rearrangement(i + 1, item)
            for i, item in enumerate(result.rearrangements)
        ]
        view_model = ClassificationViewModel(
            reference,
            unique_mapped_reads,
            rearrangements,
            program_name=cli_details.program_name,
            version=cli_details.program_version,
            cli_command=cli_details.cli_command,
        )
        return view_model

    def convert_rearrangement(
        self, i: int, rearrangement: result_models.RearrangementDTO
    ) -> "ViewRearrangement":
        id_ = self._formatter.format_id(i)
        A = self.convert_genomic_region_set(rearrangement.A)
        B = self.convert_genomic_region_set(rearrangement.B)
        evidence = self.convert_evidence(rearrangement.evidence)
        result = ViewRearrangement(id_, A, B, evidence)
        return result

    def convert_genomic_region_set(
        self, region_set: result_models.CompoundRegionDTO
    ) -> "ViewCompoundRegion":
        regions = [
            ViewRegion(**dataclasses.asdict(region)) for region in region_set.regions
        ]
        result = ViewCompoundRegion(region_set.name, regions)
        return result

    def convert_evidence(self, evidence: result_models.ReadsEvidence) -> "ViewEvidence":
        SRPB = round(evidence.SRPB, 2)
        result = ViewEvidence(evidence.paired, evidence.split, SRPB)
        return result


class ClassificationPresenter(classify_interactor.ClassifyPresenter):
    def __init__(
        self,
        classification_view: ClassificationView,
        cli_introspection: introspection.CliIntrospection,
    ) -> None:
        self._classification_view = classification_view
        self._cli_introspection = cli_introspection
        self._converter = ResultModelConverter()

    def present_classification(self, result: result_models.ClassifyResult) -> None:
        entrypoint_details = self._cli_introspection.get_entrypoint_details()
        view_model = self._converter.convert_result_model(result, entrypoint_details)
        self._classification_view.update_classification_model(view_model)


class MultiPresenter(classify_interactor.ClassifyPresenter):
    """A concrete Presenter which delegates"""

    def __init__(
        self,
        classification_presenter: classify_interactor.ClassifyPresenter,
        reads_presenter: classify_interactor.ClassifyPresenter,
    ):
        self._classification_presenter = classification_presenter
        self._reads_presenter = reads_presenter

    def present_classification(self, result: result_models.ClassifyResult) -> None:
        """Present the result from classify_reads"""
        self._classification_presenter.present_classification(result)
        self._reads_presenter.present_classification(result)


class CliIntrospectionPresenter:
    def __init__(self, version_caller: introspection.VersionCaller) -> None:
        self.__version_caller = version_caller

    def present_version(self) -> None:
        print(self.__version_caller.get_version())
