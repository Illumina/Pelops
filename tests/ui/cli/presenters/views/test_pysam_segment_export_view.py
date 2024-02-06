"""Test for the views.

Views are the only part of the codebase that does not require a unit test, as the
details are very variable and tedious/hard to test. To "test" the view
developers usually "look" at the output and make sure it has the expected
layout. The test here is not really a test - no assert - but it will create a
view content out of a provided view model, so that the manual step of looking
at the display does not require running all the business logic side
"""

import copy

import pytest

from ilmn.pelops.ui.cli.presenters import cli_presenter
from ilmn.pelops.ui.cli.presenters.views import pysam_segment_export_view as pysamview
from ilmn.pelops.ui.cli.presenters.views import view_factories


@pytest.fixture
def output_dir(tmp_path):
    # user can define location of tmp_path when calling pytest to set the
    # directory with the result will be. Use --basetemp dir (Warning: this
    # directory is removed if it exists.)
    result = tmp_path / "out"
    result.mkdir(parents=True, exist_ok=True)
    return result


class TestNotifyPysamReadsView:
    @pytest.fixture
    def reads_view(self, bam_file, output_dir, notification_factory):
        factory = view_factories.ReadsViewFactory(
            bam_file, output_dir, notification_factory
        )
        result = factory.build(features=frozenset())
        return result

    def test_update_reads_model(self, reads_view, notification_factory, output_dir):
        notification_service = notification_factory.build()
        notify_reads_view = pysamview.NotifyPysamReadsView(
            reads_view, notification_service
        )
        notify_reads_view.update_reads_model([])
        observed = notification_service.get_notifications()
        expected = [
            f"SAM file exports complete. Results are available in folder '{output_dir}'."
        ]
        assert observed == expected


class TestPysamReadsView:
    def test_update_reads_model(self, bam_file, output_dir):
        view = pysamview.PysamReadsView(bam_file, output_dir)
        provided = [
            cli_presenter.ReadViewModel(
                id="01",
                region_names=("foo", "bar"),
                segments=[],
                program_name="pelops",
                program_version="1.2.42",
                cli_command="pelops dux4r /path/to/a/file.bam",
            ),
            cli_presenter.ReadViewModel(
                id="02",
                region_names=("CoreDux4", "UNNAMED"),
                segments=[],
                program_name="pelops",
                program_version="1.2.42",
                cli_command="pelops dux4r /path/to/a/file.bam",
            ),
        ]
        view.update_reads_model(provided)


class TestHeaderManipulator:
    test_cases = [
        pytest.param(
            (
                {"PG": [{"ID": "bwa"}]},
                ("pelops"),
                {
                    "ID": "pelops",
                    "PN": "pelops",
                    "VN": "1.2.42",
                    "CL": "pelops dux4r /path/to/a/file.bam",
                    "DS": "region_a:foo  region_b:bar",
                    "PP": "bwa",
                },
            ),
            id="simple_case",
        ),
        pytest.param(
            (
                {"PG": [{"ID": "bwa"}, {"ID": "pelops"}, {"ID": "samtools"}]},
                ("pelops"),
                {
                    "ID": "pelops.1",
                    "PN": "pelops",
                    "VN": "1.2.42",
                    "CL": "pelops dux4r /path/to/a/file.bam",
                    "DS": "region_a:foo  region_b:bar",
                    "PP": "samtools",
                },
            ),
            id="program_previously_run",
        ),
        pytest.param(
            (
                {
                    "PG": [
                        {"ID": "bwa"},
                        {"ID": "pe.lops"},
                        {"ID": "pe.lops.1"},
                        {"ID": "pe.lops.2"},
                    ]
                },
                ("pe.lops"),
                {
                    "ID": "pe.lops.3",
                    "PN": "pe.lops",
                    "VN": "1.2.42",
                    "CL": "pe.lops dux4r /path/to/a/file.bam",
                    "DS": "region_a:foo  region_b:bar",
                    "PP": "pe.lops.2",
                },
            ),
            id="program_previously_run_more_than_once_with_dot_in_name",
        ),
    ]

    @pytest.fixture(params=test_cases)
    def test_case(self, request):
        provided_header, provided_details, expected_extra = request.param
        return provided_header, provided_details, expected_extra

    @pytest.fixture
    def provided_header(self, test_case):
        provided_header, provided_details, expected_extra = test_case
        return provided_header

    @pytest.fixture
    def provided_details(self, test_case):
        provided_header, provided_details, expected_extra = test_case
        program_name = provided_details

        result = cli_presenter.ReadViewModel(
            id="01",
            region_names=("foo", "bar"),
            segments=[],
            program_name=program_name,
            program_version="1.2.42",
            cli_command=program_name + " dux4r /path/to/a/file.bam",
        )
        return result

    @pytest.fixture
    def expected_extra(self, test_case):
        provided_header, provided_details, expected_extra = test_case
        return expected_extra

    def test_add_program_details(
        self, provided_header, provided_details, expected_extra
    ):
        expected = copy.deepcopy(provided_header)
        expected["PG"].append(expected_extra)
        header_manipulator = pysamview.HeaderManipulator(provided_header)
        header_manipulator.add_program_details(provided_details)
        observed = header_manipulator.get_header()
        assert observed == expected
