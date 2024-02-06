import dataclasses
import pathlib
from typing import List

# mypy: no_warn_unused_ignores
try:
    from importlib import metadata as import_metadata  # type: ignore
except ImportError:
    # metadata is added to importlib in CPython 3.8
    # See https://pypi.org/project/importlib-metadata/
    import importlib_metadata as import_metadata  # type: ignore


@dataclasses.dataclass
class EntrypointDetails:
    program_name: str
    program_version: str
    cli_command: str


class VersionCaller:
    """Find the version of this package"""

    def get_version(self) -> str:
        return import_metadata.version("ilmn-pelops")


class CliIntrospection:
    def __init__(self, version_caller: VersionCaller, cli_args: List[str]):
        self._version_caller = version_caller
        self._cli_args = cli_args

    def get_entrypoint_details(self) -> EntrypointDetails:
        version = self._version_caller.get_version()
        command = " ".join(self._cli_args)
        program_name = pathlib.Path(self._cli_args[0]).name
        result = EntrypointDetails(
            program_name=program_name, program_version=version, cli_command=command
        )
        return result
