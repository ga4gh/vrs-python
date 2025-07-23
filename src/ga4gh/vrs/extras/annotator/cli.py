"""Define command-line interface for VRS annotator tool.

$ vrs-annotate vcf input.vcf.gz --vcf_out output.vcf.gz --vrs_pickle_out vrs_objects.pkl

"""

import logging
from collections.abc import Callable
from enum import Enum
from pathlib import Path
from timeit import default_timer as timer

import click

from ga4gh.vrs.dataproxy import create_dataproxy
from ga4gh.vrs.extras.annotator.vcf import VcfAnnotator, VcfAnnotatorArgsError

_logger = logging.getLogger(__name__)


@click.group()
def _cli() -> None:
    """Annotate input files with VRS variation objects."""
    logging.basicConfig(
        filename="vrs-annotate.log",
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )


class _LogLevel(str, Enum):
    """Define legal values for `--log-level` option."""

    DEBUG = "debug"
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    CRITICAL = "critical"


def _log_level_option(func: Callable) -> Callable:
    """Provide reusable log level CLI option decorator.

    Adds a `--log-level` CLI option to any decorated command. Doesn't pass on any
    values, just sets the logging level for this module.

    :param func: incoming click command
    :return: same command, wrapped with log level option
    """

    def _set_log_level(ctx: dict, param: str, value: _LogLevel) -> None:  # noqa: ARG001
        level_map = {
            _LogLevel.DEBUG: logging.DEBUG,
            _LogLevel.INFO: logging.INFO,
            _LogLevel.WARNING: logging.WARNING,
            _LogLevel.ERROR: logging.ERROR,
            _LogLevel.CRITICAL: logging.CRITICAL,
        }
        logging.getLogger(__name__).setLevel(level_map[value])

    return click.option(
        "--log-level",
        type=click.Choice([v.value for v in _LogLevel.__members__.values()]),
        default="info",
        help="Set the logging level.",
        callback=_set_log_level,
        expose_value=False,
        is_eager=True,
    )(func)


class PathOrDash(click.ParamType):
    """click ParamType to support converting to Path and allowing '-' for stdin/out"""

    name = "path-or-dash"

    def __init__(self, **kwargs) -> None:
        """Initialize PathOrDash"""
        self._path_type = kwargs.pop("path_type", Path)
        self._inner = click.Path(**kwargs, path_type=self._path_type)

    def convert(
        self, value: str, param: click.Parameter | None, ctx: click.Context | None
    ) -> str | Path:
        """Convert the value using click.Path if it is not '-'"""
        if value == "-":
            return value
        # Use click.Path for validation
        return self._inner.convert(value, param, ctx)


@_cli.command(name="vcf")
@_log_level_option
@click.argument(
    "vcf-in",
    nargs=1,
    type=PathOrDash(
        allow_dash=True,
        exists=True,
        readable=True,
        dir_okay=False,
        path_type=Path,
    ),
)
@click.option(
    "--vcf-out",
    required=False,
    type=PathOrDash(writable=True, allow_dash=True, path_type=Path),
    help=(
        'Declare save location for output annotated VCF (or "-" to write to stdout). At least one form of output must be declared.'
    ),
)
@click.option(
    "--pkl-out",
    required=False,
    type=click.Path(writable=True, allow_dash=False, path_type=Path),
    help=(
        "Declare save location for output PKL file mapping VRS IDs to alleles. At least one form of output must be declared."
    ),
)
@click.option(
    "--ndjson-out",
    required=False,
    type=click.Path(writable=True, allow_dash=False, path_type=Path),
    help=(
        "Declare save location for output NDJSON file dump of VRS alleles. At least one form of output must be declared."
    ),
)
@click.option(
    "--vrs-attributes",
    is_flag=True,
    default=False,
    help="Include VRS_Start, VRS_End, and VRS_State fields in the VCF output INFO field.",
)
@click.option(
    "--dataproxy-uri",
    required=False,
    help="URI declaring source of sequence data. See subcommand description for more information.",
    show_default=True,
)
@click.option(
    "--assembly",
    required=False,
    default="GRCh38",
    show_default=True,
    help="Specify assembly that was used to create input VCF.",
    type=str,
)
@click.option(
    "--skip-ref",
    is_flag=True,
    default=False,
    help="Skip VRS computation for REF alleles.",
)
@click.option(
    "--require-validation",
    is_flag=True,
    default=False,
    help="Require validation checks to pass to construct a VRS object.",
)
@click.option(
    "--silent",
    "-s",
    is_flag=True,
    default=False,
    help="Suppress messages printed to stdout",
)
def _annotate_vcf_cli(
    vcf_in: Path | str,
    vcf_out: Path | str | None,
    pkl_out: Path | None,
    ndjson_out: Path | None,
    vrs_attributes: bool,
    dataproxy_uri: str,
    assembly: str,
    skip_ref: bool,
    require_validation: bool,
    silent: bool,
) -> None:
    """Extract VRS objects from VCF located at VCF_IN. VCF_IN can be "-" to read from stdin.

        $ vrs-annotate vcf input.vcf.gz --vcf-out output.vcf.gz --pkl-out vrs_objects.pkl

    Note that at least one of --vcf-out, --pkl-out, or --ndjson-out must be selected and
    defined; otherwise, this process will terminate immediately.

    Sequence data from a provider such as SeqRepo is required. Use the `--dataproxy-uri`
    option or the environment variable `GA4GH_VRS_DATAPROXY_URI` to define its location
    (the former will take priority over the latter when both are set). This option
    falls back upon the VRS-Python library default `seqrepo+http://localhost:5000/seqrepo`
    value if it is not otherwise declared.

    Currently accepted URI schemes:

    \b
     * seqrepo+http://localhost:5000/seqrepo
     * seqrepo+https://somewhere:5000/seqrepo
     * seqrepo+file:///path/to/seqrepo/root
     * seqrepo+:../relative/path/to/seqrepo/root
    """  # noqa: D301
    data_proxy = create_dataproxy(dataproxy_uri)
    annotator = VcfAnnotator(data_proxy)
    start = timer()
    msg = f"Annotating {vcf_in} with the VCF Annotator..."
    _logger.info(msg)

    try:
        annotator.annotate(
            vcf_in,
            output_vcf_path=vcf_out,
            vrs_attributes=vrs_attributes,
            assembly=assembly,
            compute_for_ref=(not skip_ref),
            require_validation=require_validation,
            output_pkl_path=pkl_out,
            output_ndjson_path=ndjson_out,
        )
    except VcfAnnotatorArgsError:
        msg = "No VCF, PKL, or NDJSON output path provided -- must set at least one of --vcf-out, --pkl-out, or --ndjson-out"
        if not silent:
            click.echo(msg)
        _logger.exception(msg)
        exit(1)

    end = timer()
    msg = f"VCF Annotator finished in {(end - start):.5f} seconds"
    _logger.info(msg)
    if not silent:
        click.echo(msg)
