import argparse
import pathlib
import warnings
from typing import Any, Dict, List

from ilmn.pelops import defaults, request_models
from ilmn.pelops.factories import interactor_factories


def add_custom_help(parser: argparse.ArgumentParser) -> None:
    """Add a custom help message, so we have control over details"""
    parser.add_argument(
        "--help",
        "-h",
        help="Show this help message and exit.",
        action="help",
        default=argparse.SUPPRESS,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pelops",
        description=f"Find spanning reads of genomic rearrangements.",
        add_help=False,
    )
    add_custom_help(parser)
    subparsers = parser.add_subparsers(title="Available actions", dest="action")

    version_parser = subparsers.add_parser(
        "version",
        help="""Print out package version""",
        add_help=False,
    )
    dux4r_parser = subparsers.add_parser(
        "dux4r",
        help="""Find evidence of DUX4-rearrangements in BAM/CRAM file. Enter `%(prog)s dux4r --help` for more info.""",
        add_help=False,
    )
    populate_classify_parser(dux4r_parser)
    return parser


def populate_classify_parser(classify_parser: argparse.ArgumentParser) -> None:
    classify_parser.add_argument("infile", help="Path to input BAM/CRAM file.")
    add_custom_help(classify_parser)
    classify_parser.add_argument(
        "--json",
        help="Path to the output json file. [DEFAULT=%(default)s]",
        default="pelops_results.json",
        metavar="FILE",
    )
    classify_parser.add_argument(
        "--export",
        help="""Path to the output folder. If provided, supporting reads for each
        rearrangement will be saved in SAM files in such folder.""",
        metavar="DIR",
    )
    classify_parser.add_argument(
        "--total-number-reads",
        type=int,
        help="""Number of reads to use for normalisation.
            If not provided, the total number of unique and mapped reads in the
            BAM file will be used.""",
        metavar="INT",
    )
    classify_parser.add_argument(
        "--threads",
        type=int,
        help="""Number of threads to use when computing total number of reads
        in the input file. [DEFAULT=%(default)s]""",
        metavar="INT",
        default=1,
    )
    classify_parser.add_argument(
        "--srpb-threshold",
        type=float,
        help="""Minimum number of Spanning Read Pairs per Billion required for
        a non-IGH DUX4-rearrangement to be called. [DEFAULT=%(default)s]""",
        default=defaults.srpb_threshold,
        metavar="FLOAT",
    )
    classify_parser.add_argument(
        "--only-igh-dux4",
        help="""If provided, it will only call rearrangements between IGH and
        DUX4 (Core and Extended).""",
        action="store_true",
    )
    classify_parser.add_argument(
        "--minimum-mapq",
        type=int,
        help="""Minimum mapping quality of reads outside of DUX4 region required
            to be counted as spanning read in non-IGH DUX4 rearrangements.
            [DEFAULT=%(default)s]""",
        metavar="INT",
        default=defaults.minimum_mapq,
    )
    classify_parser.add_argument(
        "--filter-regions",
        help="""BED file of regions to ignore when calling non-IGH DUX4-rearrangements.""",
        metavar="FILE",
    )
    classify_parser.add_argument(
        "--silent", help="Disable logging", default=False, action="store_true"
    )
    classify_parser.add_argument(
        "--with-experimental-features",
        help=argparse.SUPPRESS,
        nargs="*",
        default=[],
        choices=[x.value for x in request_models.Feature],
    )


class CliController:
    def __init__(self, interactor_factory: interactor_factories.InteractorFactory):
        self._interactor_factory = interactor_factory

    def dispatch(self, args: List[str]) -> None:
        def custom_warning_formatting(message, category, filename, lineno, file=None, line=None):  # type: ignore
            return f"\n{category.__name__}:{message}\n\n"

        parser = build_parser()
        parsed_args = parser.parse_args(args[1:])
        if parsed_args.action == "version":
            introspection = self._interactor_factory.build_introspection_interactor()
            introspection.present_version()
        else:
            interactor_factory_args = self._get_factory_args(parsed_args, args)
            request = self._get_request(parsed_args)
            interactor = self._interactor_factory.build(**interactor_factory_args)
            interactor.present_rearrangement_evidence(request)

    def _get_factory_args(
        self, parsed_args: argparse.Namespace, args: List[str]
    ) -> Dict[str, Any]:
        factory_args: Dict[str, Any] = {}
        factory_args["bam_file"] = pathlib.Path(parsed_args.infile)
        factory_args["output_json"] = pathlib.Path(parsed_args.json)
        factory_args["cli_args"] = args
        factory_args["silent"] = parsed_args.silent
        if getattr(parsed_args, "threads") is not None:
            factory_args["number_of_threads"] = int(parsed_args.threads)
        if getattr(parsed_args, "export"):
            factory_args["output_dir"] = pathlib.Path(parsed_args.export)
        if getattr(parsed_args, "filter_regions") is not None:
            factory_args["bedfile"] = pathlib.Path(parsed_args.filter_regions)
        return factory_args

    def _get_request(
        self, parsed_args: argparse.Namespace
    ) -> request_models.ClassifyRequest:
        features = []
        if not parsed_args.silent:
            features.append(request_models.Feature.WITH_NOTIFICATIONS)
        if not parsed_args.only_igh_dux4:
            features.append(request_models.Feature.DUX4_OTHER)
        if request_models.Feature.DUX4_OTHER in features:
            srpb_threshold = parsed_args.srpb_threshold
            minimum_mapping_quality = parsed_args.minimum_mapq
        else:
            srpb_threshold = minimum_mapping_quality = None

        total_number_of_reads = parsed_args.total_number_reads
        if total_number_of_reads is not None:
            features.append(request_models.Feature.PROVIDED_READ_COUNT)
        if getattr(parsed_args, "filter_regions") is not None:
            features.append(request_models.Feature.WITH_BLACKLIST)
        if parsed_args.with_experimental_features:
            pass  # No experimental features are currently supported
        request = request_models.ClassifyRequest(
            features=frozenset(features),
            minimum_mapping_quality=minimum_mapping_quality,
            srpb_threshold=srpb_threshold,
            total_number_of_reads=total_number_of_reads,
        )
        return request
