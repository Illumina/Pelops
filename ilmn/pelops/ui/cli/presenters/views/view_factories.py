import enum
import pathlib
from typing import FrozenSet, Optional

from ilmn.pelops import notifications
from ilmn.pelops.ui.cli.presenters import cli_presenter
from ilmn.pelops.ui.cli.presenters.views import (
    json_classification_export_view,
    pysam_segment_export_view,
)


class ResultFileExistsError(FileExistsError):
    pass


class InvalidDirectoryError(IOError):
    pass


def _raise_missing_directory(directory: pathlib.Path) -> None:
    message = f"Directory {directory} does not exist. Aborting"
    raise InvalidDirectoryError(message)


def _check_output(output_json: pathlib.Path) -> None:
    parent_directory = output_json.parent
    if output_json.exists():
        message = f"Result file {output_json} already exists. Aborting"
        raise ResultFileExistsError(message)
    elif not parent_directory.exists():
        _raise_missing_directory(parent_directory)


class ClassifyViewFeature(enum.Enum):
    WITH_NOTIFICATIONS = enum.auto()


class ClassifyViewFactory:
    def __init__(
        self,
        output_json: pathlib.Path,
        notification_factory: notifications.NotificationServiceFactory,
    ):
        self.__output_json = output_json
        self.__notification_factory = notification_factory

    def build(
        self, features: FrozenSet[ClassifyViewFeature]
    ) -> cli_presenter.ClassificationView:
        _check_output(self.__output_json)
        view = json_classification_export_view.JsonClassificationExportView(
            self.__output_json
        )
        if ClassifyViewFeature.WITH_NOTIFICATIONS in features:
            notification_service = self.__notification_factory.build()
            return json_classification_export_view.NotifyClassificationView(
                view, notification_service
            )
        else:
            return view


class ReadsViewFactory:
    def __init__(
        self,
        bam_file: pathlib.Path,
        output_dir: Optional[pathlib.Path],
        notification_factory: notifications.NotificationServiceFactory,
    ):
        self.__bam_file = bam_file
        self.__output_dir = output_dir
        self.__notification_factory = notification_factory

    def build(
        self, features: FrozenSet[ClassifyViewFeature]
    ) -> cli_presenter.ReadsView:
        if self.__output_dir is None:
            raise ValueError()
        elif not self.__output_dir.exists():
            _raise_missing_directory(self.__output_dir)
        view = pysam_segment_export_view.PysamReadsView(
            self.__bam_file, self.__output_dir
        )
        if ClassifyViewFeature.WITH_NOTIFICATIONS in features:
            notification_service = self.__notification_factory.build()
            result = pysam_segment_export_view.NotifyPysamReadsView(
                view, notification_service
            )
            return result
        else:
            return view
