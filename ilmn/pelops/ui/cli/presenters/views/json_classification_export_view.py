import dataclasses
import json
import pathlib
from typing import Dict

from ilmn.pelops import notifications
from ilmn.pelops.ui.cli.presenters import cli_presenter


class JsonClassificationExportView(cli_presenter.ClassificationView):
    def __init__(self, output: pathlib.Path):
        self._output = output

    def update_classification_model(
        self, result: cli_presenter.ClassificationViewModel
    ) -> None:
        with open(self._output, "w") as fh:
            json.dump(dataclasses.asdict(result), fh, indent=4)

    def get_output_file(self) -> pathlib.Path:
        return self._output


class NotifyClassificationView(cli_presenter.ClassificationView):
    def __init__(
        self,
        view_model: JsonClassificationExportView,
        notification_service: notifications.NotificationService,
    ) -> None:
        self.__view_model = view_model
        self.__notification_service = notification_service

    def update_classification_model(
        self, result: cli_presenter.ClassificationViewModel
    ) -> None:
        self.__view_model.update_classification_model(result)
        output_file = self.__view_model.get_output_file()
        message = f"Analysis complete. Results are available in file '{output_file}'."
        self.__notification_service.notify(message)
