import pathlib
from unittest import mock

import pytest

from ilmn.pelops import notifications
from ilmn.pelops.ui.cli.presenters import cli_presenter
from ilmn.pelops.ui.cli.presenters.views import json_classification_export_view


class TestNotifyClassificationView:
    @pytest.fixture
    def provided(self):
        result = cli_presenter.ClassificationViewModel(
            reference="GRch38",
            unique_mapped_reads=12345,
            rearrangements=[],
            program_name="pelops",
            version="0.3",
            cli_command="pelops dux4r /file.bam",
        )
        return result

    @pytest.fixture
    def output_json(self, tmp_path):
        return tmp_path / "out.json"

    @pytest.fixture
    def inner_view(self, tmp_path, output_json):
        result = mock.Mock(
            spec=json_classification_export_view.JsonClassificationExportView
        )
        result.get_output_file = mock.Mock(return_value=output_json)
        return result

    def test_update_classification_model(self, provided, inner_view, output_json):
        notification_service = notifications.SimpleNotificationService()
        view = json_classification_export_view.NotifyClassificationView(
            inner_view, notification_service
        )
        view.update_classification_model(provided)
        inner_view.update_classification_model.assert_called_with(provided)
        expected_message = (
            f"Analysis complete. Results are available in file '{output_json}'."
        )
        assert notification_service.get_notifications() == [expected_message]
