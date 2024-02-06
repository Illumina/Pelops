import dataclasses
from typing import FrozenSet
from unittest import mock

import pytest

from ilmn.pelops import entities, result_models
from ilmn.pelops.ui.cli.presenters import cli_presenter, introspection


@pytest.fixture
def segments_view():
    return mock.Mock(cli_presenter.ReadsView, isinstance=True)


@pytest.fixture
def some_segments():
    segment_one = result_models.PlacedSegmentDTO("ABC")
    segment_two = result_models.PlacedSegmentDTO("CDE")
    segment_three = result_models.PlacedSegmentDTO("GHI")
    result = [segment_two, segment_one, segment_three]
    return result


@pytest.fixture
def a_set():
    location = frozenset([result_models.GenomicRegionDTO("chr4", 190020407, 190023665)])
    return result_models.CompoundRegionDTO("CoreDUX4", location)


@pytest.fixture
def b_set():
    location = frozenset(
        [result_models.GenomicRegionDTO("chr14", 105586937, 106879844)]
    )
    return result_models.CompoundRegionDTO("IGH", location)


@pytest.fixture
def sorted_segments(some_segments):
    return [some_segments[i].content for i in (1, 0, 2)]


@pytest.fixture
def result_model(a_set, b_set, some_segments):
    evidence = result_models.ReadsEvidence(1, 0, 100)
    result = result_models.ClassifyResult(
        result_models.ReferenceGenome.GRCh38,
        1000000,
        [result_models.RearrangementDTO(a_set, b_set, evidence, some_segments)],
    )
    return result


class TestReadsPresenter:
    def test_present_classification(
        self, result_model, sorted_segments, segments_view, version_caller
    ):
        cli_args = ["pelops", "dux4r", "/somewhere.bam"]
        cli_introspection = introspection.CliIntrospection(
            version_caller, cli_args=cli_args
        )
        converter = cli_presenter.ReadViewModelConverter(cli_introspection)
        presenter = cli_presenter.ReadsPresenter(segments_view, converter)
        expected = [
            cli_presenter.ReadViewModel(
                id="01",
                region_names=("CoreDUX4", "IGH"),
                segments=sorted_segments,
                program_name="pelops",
                program_version="0.1.2",
                cli_command="pelops dux4r /somewhere.bam",
            )
        ]
        presenter.present_classification(result_model)
        segments_view.update_reads_model.assert_called_with(expected)


@pytest.fixture
def version_caller():
    result = mock.Mock(introspection.VersionCaller, isinstance=True)
    result.get_version = mock.MagicMock(return_value="0.1.2")
    return result


class TestClassificationPresenter:
    @pytest.fixture
    def cli_introspection(self, version_caller):
        cli_args = [
            "/somewhere/pelops",
            "dux4r",
            "data/input.bam",
            "--json",
            "output/results.json",
        ]
        result = introspection.CliIntrospection(version_caller, cli_args=cli_args)
        return result

    @pytest.fixture
    def classification_view(self):
        return mock.Mock(cli_presenter.ClassificationView, isinstance=True)

    @pytest.fixture
    def result_model(self):
        evidence = result_models.ReadsEvidence(1, 0, 100)
        result = result_models.ClassifyResult(
            result_models.ReferenceGenome.GRCh38, 1000000, []
        )
        return result

    def test_present_classification(
        self, result_model, classification_view, cli_introspection
    ):
        expected = cli_presenter.ClassificationViewModel(
            reference="GRCh38",
            unique_mapped_reads=1000000,
            rearrangements=[],
            program_name="pelops",
            version="0.1.2",
            cli_command="/somewhere/pelops dux4r data/input.bam --json output/results.json",
        )
        presenter = cli_presenter.ClassificationPresenter(
            classification_view, cli_introspection
        )
        presenter.present_classification(result_model)
        classification_view.update_classification_model.assert_called_with(expected)


class TestMultiPresenter:
    @pytest.fixture
    def classification_view(self):
        return mock.Mock(cli_presenter.ClassificationView, isinstance=True)

    @pytest.fixture
    def cli_introspection(self, version_caller):
        cli_args = ["pelops", "dux4r", "input.bam"]
        result = introspection.CliIntrospection(version_caller, cli_args=cli_args)
        return result

    @pytest.fixture
    def presenter(self, classification_view, segments_view, cli_introspection):
        classification_presenter = cli_presenter.ClassificationPresenter(
            classification_view, cli_introspection
        )

        converter = cli_presenter.ReadViewModelConverter(cli_introspection)
        reads_presenter = cli_presenter.ReadsPresenter(segments_view, converter)
        return cli_presenter.MultiPresenter(classification_presenter, reads_presenter)

    @pytest.mark.parametrize(
        "evidence, expected_evidence",
        [
            pytest.param(
                result_models.ReadsEvidence(1, 0, 100),
                {"paired_reads": 1, "split_reads": 0, "SRPB": 100},
                id="pair_split_reads",
            ),
            pytest.param(
                result_models.ReadsEvidence(1, 0, 100.0099),
                {"paired_reads": 1, "split_reads": 0, "SRPB": 100.01},
                id="2_decimal_rounding",
            ),
        ],
    )
    def test_present_reads_classification(
        self,
        classification_view,
        presenter,
        evidence,
        expected_evidence,
        a_set,
        b_set,
    ):
        provided_classify_result = result_models.ClassifyResult(
            result_models.ReferenceGenome.GRCh38,
            1000000,
            [result_models.RearrangementDTO(a_set, b_set, evidence, set())],
        )
        presenter.present_classification(provided_classify_result)

        view_a_regions = [cli_presenter.ViewRegion("chr4", 190020407, 190023665)]
        A = cli_presenter.ViewCompoundRegion("CoreDUX4", view_a_regions)
        view_b_regions = [cli_presenter.ViewRegion("chr14", 105586937, 106879844)]
        B = cli_presenter.ViewCompoundRegion("IGH", view_b_regions)
        view_evidence = cli_presenter.ViewEvidence(**expected_evidence)
        view_rearrangement = cli_presenter.ViewRearrangement("01", A, B, view_evidence)

        expected = cli_presenter.ClassificationViewModel(
            "GRCh38",
            1000000,
            [view_rearrangement],
            "pelops",
            "0.1.2",
            "pelops dux4r input.bam",
        )
        classification_view.update_classification_model.assert_called_with(expected)


class TestCliIntrospectionPresenter:
    def test_present_version(self, version_caller, capsys):
        presenter = cli_presenter.CliIntrospectionPresenter(version_caller)
        presenter.present_version()
        captured = capsys.readouterr()
        assert captured.out == "0.1.2\n"
