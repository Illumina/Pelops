import enum
import pathlib
from typing import List, Optional

from ilmn.pelops import notifications
from ilmn.pelops.callers import caller_factories
from ilmn.pelops.infrastructure import blacklist_region_repository, pysam_repositories
from ilmn.pelops.interactors import classify_interactor
from ilmn.pelops.ui.cli.presenters import (
    cli_presenter,
    introspection,
    presenter_factories,
)


class InteractorFactory:
    def build_introspection_interactor(
        self,
    ) -> cli_presenter.CliIntrospectionPresenter:
        version_caller = introspection.VersionCaller()
        presenter = cli_presenter.CliIntrospectionPresenter(version_caller)
        return presenter

    def build(
        self,
        bam_file: pathlib.Path,
        output_json: pathlib.Path,
        number_of_threads: int = 1,
        output_dir: Optional[pathlib.Path] = None,
        cli_args: Optional[List[str]] = None,
        bedfile: Optional[pathlib.Path] = None,
        silent: bool = False,
    ) -> classify_interactor.ClassifyInteractor:
        notification_factory = notifications.SimpleNotificationServiceFactory(
            silent=False
        )
        segment_repo_factory = pysam_repositories.PysamSegmentRepositoryFactory(
            bam_file, notification_factory, number_of_threads
        )

        region_repo_factory = blacklist_region_repository.FileRegionRepositoryFactory(
            bedfile
        )
        candidate_region_caller_factory = caller_factories.CandidateRegionCallerFactory(
            region_repo_factory, segment_repo_factory
        )
        region_caller_factory = caller_factories.RegionPairCallerFactory(
            candidate_region_caller_factory, region_repo_factory, notification_factory
        )

        read_caller_factory = caller_factories.ReadCallerFactory(
            segment_repo_factory, notification_factory
        )
        caller_factory = caller_factories.RearrangementCallerFactory(
            segment_repo_factory, region_caller_factory, read_caller_factory
        )
        presenter_factory = presenter_factories.PresenterFactory(
            bam_template_file=bam_file,
            json_output_file=output_json,
            notification_factory=notification_factory,
            cli_args=cli_args,
            output_dir=output_dir,
        )
        presenter = presenter_factory.build_classify_presenter(
            presenter_factories.PresenterType.READ, silent
        )
        result = classify_interactor.ClassifyInteractor(
            presenter, caller_factory, segment_repo_factory
        )
        return result
