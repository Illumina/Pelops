import textwrap

import pytest

from ilmn.pelops import entities, repositories
from ilmn.pelops.infrastructure import blacklist_region_repository

# fixture bedfile is defined in conftest, whereas bed_file is dynamically
# generated here

bedfile_delimiter_testcases = [
    pytest.param("\t", id="tab"),
    pytest.param(" ", id="single space"),
    pytest.param("    ", id="mutiple spaces"),
]


@pytest.fixture
def coordinates():
    result = [
        ("chr1", "213941196", "213942363"),
        ("chr1", "213942363", "213943530"),
        ("chr1", "213943530", "213944697"),
        ("chr2", "158364697", "158365864"),
    ]
    return result


@pytest.fixture(params=bedfile_delimiter_testcases)
def bed_file_content(tmp_path, request, coordinates):
    separator = request.param
    lines = [separator.join(parts) for parts in coordinates]
    result = "\n".join(lines)
    return result


class TestBlackListRegionRepository:
    @pytest.fixture
    def bed_file(self, tmp_path, bed_file_content, header):
        file = tmp_path / "dummy.bed"
        if header:
            file.write_text(f"{header}\n" + bed_file_content)
        else:
            file.write_text(bed_file_content)
        return file

    @pytest.fixture
    def invalid_bed_file(self, tmp_path, content):
        file = tmp_path / "invalid.bed"
        file.write_text(content)
        return file

    @pytest.mark.parametrize(
        "header",
        [
            pytest.param(None, id="standard_case"),
            pytest.param("# This is a Comment", id="comment_case"),
            pytest.param(
                'track name="ItemRGBDemo" description="Item RGB demonstration" itemRgb="On"',
                id="track_line_case",
            ),
        ],
    )
    def test_get(self, bed_file, coordinates):
        repo = blacklist_region_repository.BlackListRegionRepository(bed_file)
        regions = [
            entities.GenomicRegion(chrom, int(start), int(end))
            for chrom, start, end in coordinates
        ]
        expected = entities.CompoundRegion(
            entities.RegionsName.BLACKLIST, frozenset(regions)
        )
        observed = repo.get(entities.RegionsName.BLACKLIST)
        assert expected == observed

    @pytest.mark.parametrize(
        "content",
        [
            pytest.param("chr23\t1\t3", id="invalid_chromosome"),
            pytest.param("chr1\t1\t248956425", id="invalid_range"),
            pytest.param("chr21\t3\t1", id="invalid_size"),
        ],
    )
    def test_error_raised_when_region_invalid(self, invalid_bed_file):
        repo = blacklist_region_repository.BlackListRegionRepository(invalid_bed_file)
        with pytest.raises(blacklist_region_repository.InvalidRegionError):
            repo.get(entities.RegionsName.BLACKLIST)

    def test_get_for_each_named_region_name(self, bedfile):
        repo = blacklist_region_repository.BlackListRegionRepository(bedfile)
        skip_regions = [entities.RegionsName.UNNAMED]
        region_names = [
            item for item in entities.RegionsName if item not in skip_regions
        ]
        for region in region_names:
            observed = repo.get(region)
            assert isinstance(observed, entities.CompoundRegion)


class TestFileRegionRepositoryFactory:
    @pytest.fixture
    def bed_file(self, tmp_path, bed_file_content):
        file = tmp_path / "dummy.bed"
        file.write_text(bed_file_content)
        return file

    def test_build_with_blacklist(self, bed_file, header=None):
        factory = blacklist_region_repository.FileRegionRepositoryFactory(bed_file)
        observed = factory.build(repositories.RegionRepoType.WITH_BLACKLIST)
        assert isinstance(
            observed, blacklist_region_repository.BlackListRegionRepository
        )

    def test_build_no_blacklist(self, bedfile, header=None):
        factory = blacklist_region_repository.FileRegionRepositoryFactory(bedfile)
        observed = factory.build(repositories.RegionRepoType.BUILTIN)
        assert isinstance(observed, repositories.BuiltinRegionRepository)
