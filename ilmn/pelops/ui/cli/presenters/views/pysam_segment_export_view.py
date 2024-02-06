import copy
import pathlib
from typing import Any, Dict, Iterable, List, Optional

import pysam

from ilmn.pelops import notifications, result_models
from ilmn.pelops.ui.cli.presenters import cli_presenter


def export_to_sam_file(
    segments: Iterable[pysam.AlignedSegment],
    file_path: pathlib.Path,
    template: Dict[Any, Any],
) -> None:
    outfile = pysam.AlignmentFile(str(file_path), "w", header=template)
    for segment in segments:
        if isinstance(segment, pysam.AlignedSegment):
            outfile.write(segment)
        else:
            raise NotImplementedError(
                f"unable to export segment of type: {type(segment)}"
            )


def create_named_sam_file(
    model_item: cli_presenter.ReadViewModel, output_dir: pathlib.Path
) -> pathlib.Path:
    region_name = "-".join(model_item.region_names)
    file_name = f"{model_item.id}_{region_name}.sam"
    file_path = output_dir / file_name
    return file_path


def find_next_suitable_id(name: str, increment: int, previous_ids: List[str]) -> str:
    while ".".join([name, str(increment)]) in previous_ids:
        increment = increment + 1
    return ".".join([name, str(increment)])


class HeaderManipulator:
    """Manipulate the pysam header to add/edit relevant information"""

    def __init__(self, header: Dict[Any, Any]) -> None:
        self._header = copy.deepcopy(header)  # do not change input header

    def get_header(self) -> Dict[Any, Any]:
        return self._header

    def add_program_details(self, model_item: cli_presenter.ReadViewModel) -> None:
        previous_ids = [program["ID"] for program in self._header["PG"]]
        if model_item.program_name in previous_ids:
            id_ = find_next_suitable_id(model_item.program_name, 1, previous_ids)
        else:
            id_ = model_item.program_name

        program_details = {
            "ID": id_,
            "PN": model_item.program_name,
            "VN": model_item.program_version,
            "CL": model_item.cli_command,
            "DS": f"region_a:{model_item.region_names[0]}  region_b:{model_item.region_names[1]}",
        }
        if len(self._header["PG"]):
            program_details["PP"] = self._header["PG"][-1]["ID"]
        self._header["PG"].append(program_details)


class PysamReadsView(cli_presenter.ReadsView):
    """The Reads View implementation that uses pysam to export reads to file"""

    def __init__(self, template_bam_file: pathlib.Path, output_dir: pathlib.Path):
        self._template = pysam.AlignmentFile(str(template_bam_file), "rb")
        self._output_dir = output_dir

    def update_reads_model(self, model: Iterable[cli_presenter.ReadViewModel]) -> None:
        original_header = self._template.header.to_dict()
        for item in model:
            header_manipulator = HeaderManipulator(original_header)
            header_manipulator.add_program_details(item)
            header = header_manipulator.get_header()
            file_path = create_named_sam_file(item, self._output_dir)
            export_to_sam_file(item.segments, file_path, header)

    def get_ouput_dir(self) -> pathlib.Path:
        return self._output_dir


class NotifyPysamReadsView(cli_presenter.ReadsView):
    def __init__(
        self,
        reads_view: PysamReadsView,
        notification_service: notifications.NotificationService,
    ):
        self.__reads_view = reads_view
        self.__notification_service = notification_service

    def update_reads_model(self, model: Iterable[cli_presenter.ReadViewModel]) -> None:
        self.__reads_view.update_reads_model(model)
        output_dir = self.__reads_view.get_ouput_dir()
        message = f"SAM file exports complete. Results are available in folder '{output_dir}'."
        self.__notification_service.notify(message)
