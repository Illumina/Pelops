import pathlib
from unittest import mock

import pytest

from ilmn.pelops import request_models, result_models
from ilmn.pelops.factories import interactor_factories
from ilmn.pelops.interactors import classify_interactor as classify_interactor_
from ilmn.pelops.ui.cli import controllers
from ilmn.pelops.ui.cli.presenters import cli_presenter


class TestController:
    @pytest.fixture
    def classify_result(self):
        return result_models.ClassifyResult(
            reference=result_models.ReferenceGenome.GRCh38,
            unique_mapped_reads=800000000,
            rearrangements=[],
        )

    @pytest.fixture
    def classify_interactor(self, classify_result):
        result = mock.Mock(spec=classify_interactor_.ClassifyInteractor)
        result.present_rearrangement_evidence = mock.Mock()
        return result

    @pytest.fixture
    def introspection_interactor(self):
        return mock.Mock(spec=cli_presenter.CliIntrospectionPresenter)

    @pytest.fixture
    def presenter(self):
        result = mock.Mock(spec=classify_interactor_.ClassifyPresenter)
        return result

    @pytest.fixture
    def interactor_factory(self, classify_interactor, introspection_interactor):
        result = mock.Mock(spec=interactor_factories.InteractorFactory)
        result.build.return_value = classify_interactor
        result.build_introspection_interactor.return_value = introspection_interactor
        return result

    dispatch_test_cases = [
        pytest.param(
            (
                ["pelops", "dux4r", "bamfile.bam"],
                {
                    "bam_file": pathlib.Path("bamfile.bam"),
                    "output_json": pathlib.Path("pelops_results.json"),
                    "number_of_threads": 1,
                    "silent": False,
                },
                request_models.ClassifyRequest(
                    features=frozenset(
                        [
                            request_models.Feature.DUX4_OTHER,
                            request_models.Feature.WITH_NOTIFICATIONS,
                        ]
                    ),
                    srpb_threshold=20.0,
                    minimum_mapping_quality=10,
                ),
            ),
            id="all_defaults",
        ),
        pytest.param(
            (
                # fmt: off
                [
                    "pelops", "dux4r", "bamfile.bam",
                    "--filter-regions", "/my/regions.bed",
                    "--srpb-threshold", "5.2",
                    "--minimum-mapq", "3",
                    "--total-number-reads", "12345",
                    "--threads", "4",
                    "--export", "/path/to/dir",
                    "--json", "/somefile.json",
                    "--silent",
                ],
                # fmt: on
                {
                    "bam_file": pathlib.Path("bamfile.bam"),
                    "bedfile": pathlib.Path("/my/regions.bed"),
                    "number_of_threads": 1,
                    "output_json": pathlib.Path("pelops_results.json"),
                    "number_of_threads": 4,
                    "output_dir": pathlib.Path("/path/to/dir"),
                    "output_json": pathlib.Path("/somefile.json"),
                    "silent": True,
                },
                request_models.ClassifyRequest(
                    features=frozenset(
                        [
                            request_models.Feature.DUX4_OTHER,
                            request_models.Feature.PROVIDED_READ_COUNT,
                            request_models.Feature.WITH_BLACKLIST,
                        ]
                    ),
                    srpb_threshold=5.2,
                    minimum_mapping_quality=3,
                    total_number_of_reads=12345,
                ),
            ),
            id="specify_most_options",
        ),
        pytest.param(
            (
                ["pelops", "dux4r", "bamfile.bam", "--only-igh-dux4"],
                {
                    "bam_file": pathlib.Path("bamfile.bam"),
                    "output_json": pathlib.Path("pelops_results.json"),
                    "number_of_threads": 1,
                    "silent": False,
                },
                request_models.ClassifyRequest(
                    features=frozenset([request_models.Feature.WITH_NOTIFICATIONS]),
                ),
            ),
            id="no_dux4_other",
        ),
    ]

    @pytest.fixture(params=dispatch_test_cases)
    def test_case(self, request):
        provided, interactor_factory_args, request_model = request.param
        return provided, interactor_factory_args, request_model

    @pytest.fixture
    def provided(self, test_case):
        provided, interactor_factory_args, request_model = test_case
        return provided

    @pytest.fixture
    def interactor_factory_args(self, test_case):
        provided, interactor_factory_args, request_model = test_case
        # the factory always needs to know the cli_args
        interactor_factory_args["cli_args"] = provided
        return interactor_factory_args

    @pytest.fixture
    def request_model(self, test_case):
        provided, interactor_factory_args, request_model = test_case
        return request_model

    def test_dispatch(
        self,
        provided,
        interactor_factory,
        interactor_factory_args,
        classify_interactor,
        request_model,
    ):
        controller = controllers.CliController(interactor_factory)
        controller.dispatch(provided)
        interactor_factory.build.assert_called_with(**interactor_factory_args)
        classify_interactor.present_rearrangement_evidence.assert_called_with(
            request_model
        )

    def test_classify_no_longer_supported(self, interactor_factory):
        provided = ["pelops", "classify", "bamfile.bam"]
        controller = controllers.CliController(interactor_factory)
        with pytest.raises(SystemExit) as exc:
            controller.dispatch(provided)

    def test_version(self, interactor_factory, introspection_interactor):
        controller = controllers.CliController(interactor_factory)
        provided = ["pelops", "version"]
        controller.dispatch(provided)
        interactor_factory.build_introspection_interactor.assert_called_with()
        introspection_interactor.present_version.assert_called_once()
