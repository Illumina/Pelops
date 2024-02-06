import abc
import enum
import pathlib
from typing import FrozenSet, List, Optional

from ilmn.pelops import notifications
from ilmn.pelops.interactors import classify_interactor
from ilmn.pelops.ui.cli.presenters import cli_presenter, introspection
from ilmn.pelops.ui.cli.presenters.views import view_factories


class PresenterType(enum.Enum):
    CLASSIFICATION = enum.auto()
    NULL = enum.auto()
    MULTI = enum.auto()
    READ = enum.auto()


class PresenterFeature(enum.Enum):
    WITH_NOTIFICATIONS = enum.auto()


class PresenterFactory:
    def __init__(
        self,
        bam_template_file: pathlib.Path,
        json_output_file: pathlib.Path,
        notification_factory: notifications.NotificationServiceFactory,
        cli_args: Optional[List[str]],
        output_dir: Optional[pathlib.Path] = None,
    ) -> None:
        self._cli_args = cli_args
        self._output_dir = output_dir
        self.__classify_view_factory = view_factories.ClassifyViewFactory(
            json_output_file, notification_factory
        )
        self.__reads_view_factory = view_factories.ReadsViewFactory(
            bam_template_file, output_dir, notification_factory
        )

    def build_classify_presenter(
        self, presenter_type: PresenterType, silent: bool
    ) -> classify_interactor.ClassifyPresenter:
        classification_presenter = self.build_classification_presenter(silent)
        reads_presenter = self.build_reads_presenter(silent)
        return cli_presenter.MultiPresenter(classification_presenter, reads_presenter)

    def build_classification_presenter(
        self, silent: bool = True
    ) -> cli_presenter.ClassificationPresenter:
        if self._cli_args is None:
            raise ValueError("cli_args must be provided")
        else:
            view_features = self.__get_classify_view_feature(silent)
            classification_view = self.__classify_view_factory.build(view_features)
            cli_introspection = introspection.CliIntrospection(
                introspection.VersionCaller(), self._cli_args
            )
            return cli_presenter.ClassificationPresenter(
                classification_view, cli_introspection
            )

    def build_reads_presenter(
        self, silent: bool = True
    ) -> classify_interactor.ClassifyPresenter:
        if self._output_dir is None:
            return cli_presenter.NullReadPresenter()
        else:
            if self._cli_args is None:
                raise ValueError("cli_args must be provided")
            else:
                version_caller = introspection.VersionCaller()
                cli_introspection = introspection.CliIntrospection(
                    version_caller, cli_args=self._cli_args
                )
                converter = cli_presenter.ReadViewModelConverter(cli_introspection)
                view_features = self.__get_classify_view_feature(silent)
                reads_view = self.__reads_view_factory.build(view_features)
                return cli_presenter.ReadsPresenter(reads_view, converter)

    def __get_classify_view_feature(
        self, silent: bool
    ) -> FrozenSet[view_factories.ClassifyViewFeature]:
        if silent:
            return frozenset()
        else:
            return frozenset([view_factories.ClassifyViewFeature.WITH_NOTIFICATIONS])
