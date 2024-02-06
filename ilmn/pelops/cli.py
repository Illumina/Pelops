"""Console script for ilmn_pelops."""

import sys
from typing import List

from ilmn.pelops.factories import interactor_factories
from ilmn.pelops.ui.cli import controllers


def main_from_args(args: List[str]) -> int:
    interactor_factory = interactor_factories.InteractorFactory()
    controller = controllers.CliController(interactor_factory)
    controller.dispatch(args)
    return 0


def main() -> int:
    """Console script for ilmn_pelops."""
    return main_from_args(sys.argv)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
