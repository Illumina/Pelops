import pytest

from ilmn.pelops.ui.cli.presenters import cli_presenter, presenter_factories


@pytest.fixture
def output_json(tmp_path):
    return tmp_path / "a_file.json"


@pytest.fixture
def cli_args():
    return ["pelops", "dux4r", "file.bam"]


class TestPresenterFactory:
    @pytest.fixture
    def factory(self, output_json, cli_args, bam_file, tmp_path, notification_factory):
        return presenter_factories.PresenterFactory(
            bam_template_file=bam_file,
            json_output_file=output_json,
            notification_factory=notification_factory,
            output_dir=tmp_path,
            cli_args=cli_args,
        )

    @pytest.mark.parametrize("silent", [True, False])
    def test_build_classify_presenter(self, factory, silent):
        presenter_type = presenter_factories.PresenterType.READ
        observed = factory.build_classify_presenter(presenter_type, silent=silent)
        assert isinstance(observed, cli_presenter.MultiPresenter)

    def test_build_classification_presenter(self, factory):
        observed = factory.build_classification_presenter(silent=True)
        assert isinstance(observed, cli_presenter.ClassificationPresenter)

    def test_build_reads_presenter_no_output(
        self, bam_file, output_json, cli_args, notification_factory
    ):
        factory = presenter_factories.PresenterFactory(
            bam_template_file=bam_file,
            json_output_file=output_json,
            notification_factory=notification_factory,
            cli_args=cli_args,
        )
        observed = factory.build_reads_presenter()
        assert isinstance(observed, cli_presenter.NullReadPresenter)

    def test_build_reads_presenter(self, bam_file, tmp_path, factory):
        observed = factory.build_reads_presenter()
        assert isinstance(observed, cli_presenter.ReadsPresenter)
