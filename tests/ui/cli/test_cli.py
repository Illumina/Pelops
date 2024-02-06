"""End-to-end test"""

import pytest

from ilmn.pelops import cli

# bam_file is define in conftest


class TestMain:
    test_cases = [
        pytest.param({"--export": "output_folder"}, id="with_output_folder"),
        pytest.param(
            {
                "export": "output_folder",
                "total_number_reads": "200000",
                "threads": "2",
                "srpb_threshold": "5.9",
                "minimum_mapq": "5",
                "with_bed_file": True,
            },
            id="with_all_options",
        ),
    ]

    @pytest.fixture
    def json_out(self, tmp_path):
        result = ["--json", str(tmp_path / "output.json")]
        return result

    @pytest.fixture(params=test_cases)
    def extra_args(self, tmp_path, request, bedfile):
        result = []
        if "export" in request.param:
            output_folder = tmp_path / "output_folder"
            output_folder.mkdir()
            result.extend(["--export", str(output_folder)])
        if "total_number_reads" in request.param:
            result.extend(["--total-number-reads", request.param["total_number_reads"]])
        if "threads" in request.param:
            result.extend(["--threads", request.param["threads"]])
        if "srpb_threshold" in request.param:
            result.extend(["--srpb-threshold", request.param["srpb_threshold"]])
        if "minimum_mapq" in request.param:
            result.extend(["--minimum-mapq", request.param["minimum_mapq"]])
        if "with_bed_file" in request.param:
            result.extend(["--filter-regions", str(bedfile)])

        return result

    def test_main_from_args(self, bam_file, json_out, extra_args, tmp_path):
        provided = ["pelops", "dux4r", str(bam_file)] + json_out + extra_args
        result = cli.main_from_args(provided)
        assert result == 0
