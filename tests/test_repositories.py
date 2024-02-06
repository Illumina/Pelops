import pytest

from ilmn.pelops import entities, notifications, repositories

# fixuture exclude_flags are located in conftest.py


@pytest.fixture
def segment_counter():
    result = repositories.ProvidedSegmentCounter(1234)
    return result


class TestProvidedSegmentCounter:
    def test_get_number_of_segments(self, segment_counter, exclude_flags):
        observed = segment_counter.get_number_of_segments(exclude=exclude_flags)
        assert observed == 1234

    def test_get_number_of_segments_all(self, segment_counter):
        observed = segment_counter.get_number_of_segments(exclude=[])
        assert observed == 1234


class TestNotifySegmentCounter:
    @pytest.fixture
    def notification_service(self):
        result = notifications.SimpleNotificationService()
        return result

    def test_get_number_of_segments(
        self, notification_service, segment_counter, exclude_flags
    ):
        counter_with_notification = notifications.NotifySegmentCounter(
            segment_counter, notification_service
        )
        observed = counter_with_notification.get_number_of_segments(
            exclude=exclude_flags
        )
        observed_messages = notification_service.get_notifications()
        assert observed_messages == ["Counting number of unique and mapped reads."]
        assert observed == 1234


class TestBuiltinRegionRepository:
    @pytest.fixture
    def region_repository(self):
        result = repositories.BuiltinRegionRepository()
        return result

    get_test_cases = [
        pytest.param(
            [
                entities.RegionsName.IGH,
                [("chr14", 105586937, 106879844)],
            ],
            id="IGH",
        ),
        pytest.param(
            [
                entities.RegionsName.GRCh38,
                [
                    ("chr1", 1, 248956422),
                    ("chr2", 1, 242193529),
                    ("chr3", 1, 198295559),
                    ("chr4", 1, 190214555),
                    ("chr5", 1, 181538259),
                    ("chr6", 1, 170805979),
                    ("chr7", 1, 159345973),
                    ("chr8", 1, 145138636),
                    ("chr9", 1, 138394717),
                    ("chr10", 1, 133797422),
                    ("chr11", 1, 135086622),
                    ("chr12", 1, 133275309),
                    ("chr13", 1, 114364328),
                    ("chr14", 1, 107043718),
                    ("chr15", 1, 101991189),
                    ("chr16", 1, 90338345),
                    ("chr17", 1, 83257441),
                    ("chr18", 1, 80373285),
                    ("chr19", 1, 58617616),
                    ("chr20", 1, 64444167),
                    ("chr21", 1, 46709983),
                    ("chr22", 1, 50818468),
                    ("chrX", 1, 156040895),
                    ("chrY", 1, 57227415),
                ],
            ],
            id="GRCh38",
        ),
        pytest.param(
            [
                entities.RegionsName.CoreDUX4,
                [
                    ("chr4", 190020407, 190023665),
                    ("chr4", 190066935, 190093279),
                    ("chr4", 190172774, 190176845),
                    ("chr10", 133663429, 133685936),
                    ("chr10", 133739606, 133762125),
                ],
            ],
            id="CoreDUX4",
        ),
        pytest.param(
            [
                entities.RegionsName.ExtendedDUX4,
                [
                    ("chr4", 189967935, 190204560),
                    ("chr10", 133564429, 133787422),
                    ("chr3", 75667931, 75671185),
                    ("chr5", 31248879, 31251987),
                    ("chr9", 63816748, 63819462),
                    ("chr12", 34207415, 34210675),
                    ("chr12", 61599067, 61601843),
                    ("chr16", 34134736, 34137792),
                    ("chr16", 34140256, 34143849),
                    ("chr20", 29317824, 29320579),
                    ("chr20", 29323092, 29326049),
                    ("chr20", 29409348, 29412600),
                    ("chr20", 29447517, 29450306),
                    ("chr20", 29877636, 29880363),
                    ("chrY", 10170590, 10173725),
                    ("chrY", 11305918, 11309181),
                    ("chrY", 11313921, 11317187),
                    ("chrY", 11320557, 11323823),
                    ("chrY", 11331329, 11334595),
                ],
            ],
            id="ExtendedDUX4",
        ),
    ]

    @pytest.fixture(params=get_test_cases)
    def test_case(self, request):
        name, regions = request.param
        return name, regions

    @pytest.fixture
    def name(self, test_case):
        name, regions = test_case
        return name

    @pytest.fixture
    def regionset(self, test_case):
        name, regions = test_case
        regionset = frozenset([entities.GenomicRegion(*region) for region in regions])
        result = entities.CompoundRegion(name=name, regions=regionset)
        return result

    def test_get(self, name, regionset, region_repository):
        observed = region_repository.get(name)
        expected = regionset
        assert observed == expected

    def test_get_for_each_named_region_name(self, region_repository):
        skip_regions = [entities.RegionsName.UNNAMED, entities.RegionsName.BLACKLIST]
        region_names = [
            item for item in entities.RegionsName if item not in skip_regions
        ]
        for region in region_names:
            observed = region_repository.get(region)
            assert isinstance(observed, entities.CompoundRegion)

    def test_get_raises_for_unnamed(self, region_repository):
        unnamed = entities.RegionsName.UNNAMED
        with pytest.raises(KeyError):
            region_repository.get(unnamed)
