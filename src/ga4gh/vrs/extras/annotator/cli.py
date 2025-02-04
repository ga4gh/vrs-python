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
from ga4gh.vrs.extras.annotator.vcf import VCFAnnotator

_logger = logging.getLogger(__name__)


@click.group()
def _cli() -> None:
    """Annotate input files with VRS variation objects."""
    logging.basicConfig(
        filename="vrs-annotator.log",
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )


class _LogLevel(str, Enum):
    """Define legal values for `--log_level` option."""

    DEBUG = "debug"
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    CRITICAL = "critical"


def _log_level_option(func: Callable) -> Callable:
    """Provide reusable log level CLI option decorator.

    Adds a `--log_level` CLI option to any decorated command. Doesn't pass on any
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
        "--log_level",
        type=click.Choice([v.value for v in _LogLevel.__members__.values()]),
        default="info",
        help="Set the logging level.",
        callback=_set_log_level,
        expose_value=False,
        is_eager=True,
    )(func)


@_cli.command(name="vcf")
@_log_level_option
@click.argument(
    "vcf-in",
    nargs=1,
    type=click.Path(exists=True, readable=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--vcf-out",
    required=False,
    type=click.Path(writable=True, allow_dash=False, path_type=Path),
    help=(
        "Declare save location for output annotated VCF. If not provided, must provide --vrs_pickle_out."
    ),
)
@click.option(
    "--pkl-out",
    required=False,
    type=click.Path(writable=True, allow_dash=False, path_type=Path),
    help=(
        "Declare save location for output VCF pickle. If not provided, must provide --vcf_out."
    ),
)
@click.option(
    "--incl-vrs-attrs",
    is_flag=True,
    default=False,
    help="Include VRS_Start, VRS_End, and VRS_State fields in the VCF output INFO field.",
)
@click.option(
    "--dataproxy-uri",
    required=False,
    default="seqrepo+http://localhost:5000/seqrepo",
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
    "--incl-ref-allele",
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
    vcf_in: Path,
    vcf_out: Path | None,
    pkl_out: Path | None,
    dataproxy_uri: str,
    assembly: str,
    incl_vrs_attrs: bool,
    incl_ref_allele: bool,
    require_validation: bool,
    silent: bool,
) -> None:
    """Extract VRS objects from VCF located at VCF_IN.

        $ vrs-annotate vcf input.vcf.gz --vcf_out output.vcf.gz --vrs_pickle_out vrs_objects.pkl

    Note that at least one of --vcf_out or --vrs_pickle_out must be selected and defined.

    Sequence data from a provider such as SeqRepo is required. Use the `--dataproxy_api`
    option or the environment variable `GA4GH_VRS_DATAPROXY_URI` to define its location.
    Currently accepted URI schemes:

    \b
     * seqrepo+file:///path/to/seqrepo/root
     * seqrepo+:../relative/path/to/seqrepo/root
     * seqrepo+http://localhost:5000/seqrepo
     * seqrepo+https://somewhere:5000/seqrepo
    """  # noqa: D301
    data_proxy = create_dataproxy(dataproxy_uri)
    annotator = VCFAnnotator(data_proxy)
    start = timer()
    msg = f"Annotating {vcf_in} with the VCF Annotator..."
    _logger.info(msg)
    if not silent:
        click.echo(msg)
    annotator.annotate(
        vcf_in,
        output_vcf_path=vcf_out,
        output_pkl_path=pkl_out,
        incl_vrs_attrs=incl_vrs_attrs,
        incl_ref_allele=incl_ref_allele,
        assembly=assembly,
        require_validation=require_validation,
    )
    end = timer()
    msg = f"VCF Annotator finished in {(end - start):.5f} seconds"
    _logger.info(msg)
    if not silent:
        click.echo(msg)
