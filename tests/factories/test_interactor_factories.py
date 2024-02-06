import pytest

from ilmn.pelops.factories import interactor_factories
from ilmn.pelops.interactors import classify_interactor
from ilmn.pelops.ui.cli.presenters import cli_presenter


@pytest.fixture
def output_json(tmp_path):
    return tmp_path / "a_file.json"


@pytest.fixture
def cli_args():
    return ["pelops", "dux4r", "file.bam"]


class TestInteractorFactory:
    def test_build_classify_interactor(self, bam_file, output_json, cli_args):
        factory = interactor_factories.InteractorFactory()
        observed = factory.build(
            bam_file=bam_file, output_json=output_json, cli_args=cli_args
        )
        assert isinstance(observed, classify_interactor.ClassifyInteractor)

    def test_build_classify_interactor_with_blacklist(
        self, bam_file, output_json, cli_args, bedfile
    ):
        factory = interactor_factories.InteractorFactory()
        observed = factory.build(
            bam_file=bam_file,
            output_json=output_json,
            cli_args=cli_args,
            bedfile=bedfile,
        )
        assert isinstance(observed, classify_interactor.ClassifyInteractor)

    def test_build_introspection_interactor(self):
        factory = interactor_factories.InteractorFactory()
        observed = factory.build_introspection_interactor()
        assert isinstance(observed, cli_presenter.CliIntrospectionPresenter)
