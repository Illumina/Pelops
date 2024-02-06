import time
from typing import List

import pysam
import pytest

from ilmn.pelops import entities, repositories
from ilmn.pelops.infrastructure import pysam_repositories

# fixtures exclude_flags, alignment_file, segment_repo_factory,
# test_folder  are located in conftest.py


class TestSamFlag:
    def test_query_flag(self):
        test_flag = [
            pysam_repositories.SamFlag.is_paired,
            pysam_repositories.SamFlag.is_read1,
        ]
        observed = sum(test_flag)
        assert observed == 65


class TestMapQSelector:
    test_cases = [(0, False), (1, True), (2, True)]

    @pytest.mark.parametrize("read_quality, expected", test_cases)
    def test_call(self, read_quality, expected):
        read = pysam.AlignedSegment()
        read.mapping_quality = read_quality

        selector = pysam_repositories.MapQSelector(minimum_quality=1)
        observed = selector(read)
        assert observed == expected


class TestAnySegmentSelector:
    @pytest.fixture
    def duplicated_segment(self):
        result = pysam.AlignedSegment()
        result.is_duplicate = True
        result.is_qcfail = False
        return result

    @pytest.fixture
    def properly_paired_segment(self):
        result = pysam.AlignedSegment()
        result.is_duplicate = False
        result.is_paired = True
        result.is_proper_pair = True
        return result

    def test_call_no_query(self, duplicated_segment, properly_paired_segment):
        # if we query nothing, no segment should be selected
        selector = pysam_repositories.AnySegmentSelector([])
        assert not selector(duplicated_segment)
        assert not selector(properly_paired_segment)

    def test_call(self, duplicated_segment, properly_paired_segment):
        # A segment that has Any of the queried properties should be selected
        queries = [
            repositories.ReadQuery.is_duplicate,
            repositories.ReadQuery.is_qcfail,
            repositories.ReadQuery.is_not_paired,
        ]
        selector = pysam_repositories.AnySegmentSelector(queries)
        assert selector(duplicated_segment)  # it is duplicated
        assert not selector(properly_paired_segment)  # it has no properties queried


class TestPysamPropertyError:
    def test_message(self):
        with pytest.raises(pysam_repositories.PysamPropertyError) as exc:
            raise pysam_repositories.PysamPropertyError("foo", 3, "bool")
        expected = (
            f"pysam.AlignedSegment.foo returned unexpected type int. Expected bool"
        )
        assert str(exc.value) == expected


@pytest.fixture
def reads_counter(alignment_file):
    non_cached_repo = pysam_repositories.BamFileSegmentCounter(alignment_file)
    result = repositories.CachedSegmentCounter(non_cached_repo)
    return result


@pytest.fixture
def segment_repository(reads_counter, alignment_file):
    result = pysam_repositories.FilePlacedSegmentRepository(
        alignment_file, reads_counter
    )
    return result


class TestCachedSegmentCounter:
    def test_get_number_of_segments(self, reads_counter, exclude_flags):
        observed = reads_counter.get_number_of_segments(exclude=exclude_flags)
        assert observed == 9061

    def test_get_number_of_segments_all(self, reads_counter):
        observed = reads_counter.get_number_of_segments(exclude=[])
        assert observed == 9216

    def test_get_number_of_segments_performance(self, reads_counter, exclude_flags):
        """Check the repository can cache results"""
        tic = time.time()
        observed = reads_counter.get_number_of_segments(exclude=[])
        toc = time.time()
        first_time = toc - tic

        tic = time.time()
        observed = reads_counter.get_number_of_segments(exclude=[])
        toc = time.time()
        second_time = toc - tic
        # second time is at least 10 times faster
        assert second_time < first_time / 10


class TestPlacedSegmentRepository:
    @pytest.fixture
    def locations(self):
        result = tuple([entities.GenomicRegion("chr4", 190066935, 190093279)])
        return result

    def test_get(self, locations, segment_repository):
        reads_to_exclude = [
            repositories.ReadQuery.is_duplicate,
            repositories.ReadQuery.is_not_paired,
            repositories.ReadQuery.is_qcfail,
        ]
        observed = list(segment_repository.get(locations, reads_to_exclude))
        assert len(observed) == 5500

    def test_get_proper_pair(self, locations, segment_repository):
        query = [
            repositories.ReadQuery.is_duplicate,
            repositories.ReadQuery.is_not_paired,
            repositories.ReadQuery.is_proper_pair,
            repositories.ReadQuery.is_qcfail,
        ]
        observed = list(segment_repository.get(locations, exclude=query))
        assert len(observed) == 146

    def test_get_proper_pair_min_quality(self, locations, segment_repository):
        observed = list(segment_repository.get(locations, exclude=[], min_quality=20))
        assert len(observed) == 13

    test_cases = [
        pytest.param(
            # this is supplementary, so should return None
            ("unmapped_mate_test.bam", "HKGGYCCXY:4:1205:3412001:0", 0, None)
        ),
        pytest.param(
            (
                # This is unmapped, should return the third one
                ("unmapped_mate_test.bam", "HKGGYCCXY:4:1205:3412001:0", 1, 2)
                # "HKGGYCCXY:4:1205:3412001:0",
            )
        ),
        pytest.param(
            (
                # this requests the unmpaped one (second one with flag 165)
                ("unmapped_mate_test.bam", "HKGGYCCXY:4:1205:3412001:0", 2, 1)
            )
        ),
        pytest.param(
            (
                # this requests the unmpaped one (133)
                ("unmapped_mate_test_2.bam", "HWFMFCCXX:7:2214:1340522:0", 0, 1)
            )
        ),
        pytest.param(
            (
                # This is unmapped, should return the first one
                ("unmapped_mate_test_2.bam", "HWFMFCCXX:7:2214:1340522:0", 1, 0)
            )
        ),
        pytest.param(
            (
                # this is supplementary, so should return None
                ("unmapped_mate_test_2.bam", "HWFMFCCXX:7:2214:1340522:0", 2, None)
            )
        ),
    ]

    @pytest.fixture(params=test_cases)
    def test_case(self, request):
        file, name, provided_index, expected_index = request.param
        return file, name, provided_index, expected_index

    @pytest.fixture
    def bam_file(self, test_case, test_folder):
        bam_file_name, _, _, _ = test_case
        result = test_folder / "data" / bam_file_name
        return result

    @pytest.fixture
    def repository(self, test_case, bam_file):
        counter = pysam_repositories.BamFileSegmentCounter(bam_file)
        repository = pysam_repositories.FilePlacedSegmentRepository(bam_file, counter)
        return repository

    @pytest.fixture
    def potential_unmapped_reads(self, test_case, bam_file):
        _, read_name, _, _ = test_case
        pysam_reader = pysam.AlignmentFile(str(bam_file), "rb")
        potential_unmapped_reads = []
        for segment in pysam_reader.fetch():
            if segment.query_name == read_name:
                potential_unmapped_reads.append(segment)
        return potential_unmapped_reads

    @pytest.fixture
    def expected(self, test_case, potential_unmapped_reads):
        _, _, _, expected_index = test_case
        if expected_index is None:
            return None
        else:
            return potential_unmapped_reads[expected_index]

    @pytest.fixture
    def provided(self, test_case, potential_unmapped_reads):
        _, _, index, _ = test_case
        read_with_unmapped_mate = potential_unmapped_reads[index]
        result = entities.PlacedSegment(
            read_with_unmapped_mate.query_name,
            frozenset(),
            entities.ReadOrder.ONE,
            read_with_unmapped_mate,
        )
        return result

    def test_get_mate_unmapped(self, repository, expected, provided):
        observed = repository.get_mate(provided)
        if observed is None:
            assert expected == observed
        else:
            assert observed.content == expected

    def test_get_mate_exact_position(self, segment_repository, locations):
        # get a specific read we know about
        reads = [
            read
            for read in segment_repository.get(locations, exclude=[])
            if read.read_name == "HSQ1008:146:C0JD1ACXX:2:1206:15914:105605"
        ]
        assert len(reads) == 1
        chrom, pos = segment_repository.get_mate_exact_position(reads[0])
        assert chrom == "chr14"
        assert pos == 105607021


class TestSegmentCounterFactory:
    counter_test_cases = [
        pytest.param(
            frozenset(),
            None,
            repositories.CachedSegmentCounter,
            id="real file",
        ),
        pytest.param(
            frozenset([repositories.SegmentRepoFeature.WITH_NOTIFICATION]),
            None,
            repositories.CachedSegmentCounter,
            id="with logging",
        ),
        pytest.param(
            frozenset([repositories.SegmentRepoFeature.BUILTIN]),
            1234,
            repositories.ProvidedSegmentCounter,
            id="provided counts",
        ),
        pytest.param(
            frozenset(
                [
                    repositories.SegmentRepoFeature.WITH_NOTIFICATION,
                    repositories.SegmentRepoFeature.BUILTIN,
                ]
            ),
            1234,
            repositories.ProvidedSegmentCounter,
            id="no notifications with builtin",
        ),
    ]

    @pytest.mark.parametrize(
        "features, total_number_of_reads, expected", counter_test_cases
    )
    def test_build_counter(
        self, features, total_number_of_reads, expected, segment_repo_factory
    ):
        repo = segment_repo_factory.build_counter(features, total_number_of_reads)
        assert isinstance(repo, expected)

    test_cases = [
        pytest.param(
            frozenset(),
            None,
            pysam_repositories.FilePlacedSegmentRepository,
            id="real file",
        ),
        pytest.param(
            frozenset(
                [
                    repositories.SegmentRepoFeature.WITH_NOTIFICATION,
                    repositories.SegmentRepoFeature.BUILTIN,
                ]
            ),
            1234,
            pysam_repositories.FilePlacedSegmentRepository,
            id="only wrap the counter",
        ),
    ]

    @pytest.mark.parametrize("features, total_number_of_reads, expected", test_cases)
    def test_build(
        self, features, total_number_of_reads, expected, segment_repo_factory
    ):
        repo = segment_repo_factory.build(features, total_number_of_reads)
        assert isinstance(repo, expected)
