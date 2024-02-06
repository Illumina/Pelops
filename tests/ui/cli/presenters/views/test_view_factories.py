import pytest

from ilmn.pelops.ui.cli.presenters.views import (
    json_classification_export_view,
    pysam_segment_export_view,
    view_factories,
)


@pytest.fixture
def output_json(tmp_path):
    return tmp_path / "a_file.json"


class TestClassifyViewFactory:
    test_cases = [
        pytest.param(
            frozenset(),
            json_classification_export_view.JsonClassificationExportView,
            id="no notification",
        ),
        pytest.param(
            frozenset([view_factories.ClassifyViewFeature.WITH_NOTIFICATIONS]),
            json_classification_export_view.NotifyClassificationView,
            id="with notification",
        ),
    ]

    @pytest.mark.parametrize("features, type_", test_cases)
    def test_build(self, output_json, features, type_, notification_factory):
        view_factory = view_factories.ClassifyViewFactory(
            output_json, notification_factory
        )
        observed = view_factory.build(features)
        assert isinstance(observed, type_)

    def test_build_fails_file_exists(self, output_json, notification_factory):
        output_json.touch()
        view_factory = view_factories.ClassifyViewFactory(
            output_json, notification_factory
        )
        with pytest.raises(view_factories.ResultFileExistsError) as exc:
            view_factory.build(frozenset())

    def test_build_fails_directory_missing(
        self, tmp_path, output_json, notification_factory
    ):
        outdir = tmp_path / "do_not_exist_dir"
        output_json = outdir / "a_file.json"
        view_factory = view_factories.ClassifyViewFactory(
            output_json, notification_factory
        )
        with pytest.raises(view_factories.InvalidDirectoryError) as exc:
            view_factory.build(frozenset())


class TestReadsViewFactory:
    test_cases = [
        pytest.param(
            frozenset(),
            pysam_segment_export_view.PysamReadsView,
            id="no notification",
        ),
        pytest.param(
            frozenset([view_factories.ClassifyViewFeature.WITH_NOTIFICATIONS]),
            pysam_segment_export_view.NotifyPysamReadsView,
            id="with notification",
        ),
    ]

    @pytest.mark.parametrize("features, expected", test_cases)
    def test_build(self, bam_file, tmp_path, features, expected, notification_factory):
        view_factory = view_factories.ReadsViewFactory(
            bam_file, tmp_path, notification_factory
        )
        observed = view_factory.build(features=features)
        assert isinstance(observed, expected)

    def test_build_reads_view_fails_directory_missing(
        self, bam_file, tmp_path, notification_factory
    ):
        outdir = tmp_path / "do_not_exist_dir"
        view_factory = view_factories.ReadsViewFactory(
            bam_file, outdir, notification_factory
        )
        with pytest.raises(view_factories.InvalidDirectoryError) as exc:
            view_factory.build(features=frozenset())
