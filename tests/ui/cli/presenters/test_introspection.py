from unittest import mock

import pytest

from ilmn.pelops.ui.cli.presenters import introspection


class TestCliInstrospection:
    @pytest.fixture
    def version_caller(self):
        result = mock.Mock(spec=introspection.VersionCaller, instance=True)
        result.get_version = mock.Mock(return_value="0.1.3")
        return result

    get_entripoint_test_cases = [
        pytest.param(
            [
                (["pelops", "dux4r", "/somefile.bam"]),
                ("pelops", "pelops dux4r /somefile.bam"),
            ],
            id="relative_path",
        ),
        pytest.param(
            [
                (["/path/to/pelops", "dux4r", "/somefile.bam"]),
                ("pelops", "/path/to/pelops dux4r /somefile.bam"),
            ],
            id="absolute_path",
        ),
    ]

    @pytest.fixture(params=get_entripoint_test_cases)
    def test_case(self, request):
        provided, expected = request.param
        return provided, expected

    @pytest.fixture
    def cli_args(self, test_case):
        provided, expected = test_case
        return provided

    @pytest.fixture
    def expected(self, test_case):
        provided, expected = test_case
        program_name, cli_command = expected
        result = introspection.EntrypointDetails(
            program_name=program_name, program_version="0.1.3", cli_command=cli_command
        )
        return result

    def test_get_entrypoint_details(self, version_caller, cli_args, expected):
        cli_introspection = introspection.CliIntrospection(
            version_caller, cli_args=cli_args
        )
        observed = cli_introspection.get_entrypoint_details()
        assert observed == expected
