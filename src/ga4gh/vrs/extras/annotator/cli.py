"""Define command-line interface for VRS annotator tool.

$ vrs-annotate vcf input.vcf.gz --vcf_out output.vcf.gz --vrs_pickle_out vrs_objects.pkl

"""

import logging
from collections.abc import Callable
from enum import Enum
from pathlib import Path
from timeit import default_timer as timer

import click

from ga4gh.vrs.extras.annotator.vcf import SeqRepoProxyType, VcfAnnotator

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
    "vcf_in",
    nargs=1,
    type=click.Path(exists=True, readable=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--vcf_out",
    required=False,
    type=click.Path(writable=True, allow_dash=False, path_type=Path),
    help=(
        "Declare save location for output annotated VCF. If not provided, must provide --vrs_pickle_out."
    ),
)
@click.option(
    "--vrs_pickle_out",
    required=False,
    type=click.Path(writable=True, allow_dash=False, path_type=Path),
    help=(
        "Declare save location for output VCF pickle. If not provided, must provide --vcf_out."
    ),
)
@click.option(
    "--vrs_attributes",
    is_flag=True,
    default=False,
    help="Include VRS_Start, VRS_End, and VRS_State fields in the VCF output INFO field.",
)
@click.option(
    "--seqrepo_dp_type",
    required=False,
    default=SeqRepoProxyType.LOCAL,
    type=click.Choice(
        [v.value for v in SeqRepoProxyType.__members__.values()], case_sensitive=True
    ),
    help="Specify type of SeqRepo dataproxy to use.",
    show_default=True,
    show_choices=True,
)
@click.option(
    "--seqrepo_root_dir",
    required=False,
    default=Path("/usr/local/share/seqrepo/latest"),
    type=click.Path(path_type=Path),
    help="Define root directory for local SeqRepo instance, if --seqrepo_dp_type=local.",
    show_default=True,
)
@click.option(
    "--seqrepo_base_url",
    required=False,
    default="http://localhost:5000/seqrepo",
    help="Specify base URL for SeqRepo REST API, if --seqrepo_dp_type=rest.",
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
    "--skip_ref",
    is_flag=True,
    default=False,
    help="Skip VRS computation for REF alleles.",
)
@click.option(
    "--require_validation",
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
    vrs_pickle_out: Path | None,
    vrs_attributes: bool,
    seqrepo_dp_type: SeqRepoProxyType,
    seqrepo_root_dir: Path,
    seqrepo_base_url: str,
    assembly: str,
    skip_ref: bool,
    require_validation: bool,
    silent: bool,
) -> None:
    """Extract VRS objects from VCF located at VCF_IN.

        $ vrs-annotate vcf input.vcf.gz --vcf_out output.vcf.gz --vrs_pickle_out vrs_objects.pkl

    Note that at least one of --vcf_out or --vrs_pickle_out must be selected and defined.
    """
    annotator = VcfAnnotator(
        seqrepo_dp_type, seqrepo_base_url, str(seqrepo_root_dir.absolute())
    )
    start = timer()
    msg = f"Annotating {vcf_in} with the VCF Annotator..."
    _logger.info(msg)
    if not silent:
        click.echo(msg)
    annotator.annotate(
        vcf_in.absolute(),
        output_vcf_path=vcf_out,
        output_pkl_path=vrs_pickle_out,
        vrs_attributes=vrs_attributes,
        assembly=assembly,
        compute_for_ref=(not skip_ref),
        require_validation=require_validation,
    )
    end = timer()
    msg = f"VCF Annotator finished in {(end - start):.5f} seconds"
    _logger.info(msg)
    if not silent:
        click.echo(msg)
