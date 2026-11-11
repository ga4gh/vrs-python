"""Stage authored pages, supported versions, and notebooks for Zensical."""

import json
import re
from pathlib import Path
from shutil import copy2, copytree, rmtree

from nbconvert import MarkdownExporter

from ga4gh.core.version import CORE_VERSION
from ga4gh.vrs.version import VRS_VERSION

ROOT = Path(__file__).resolve().parents[1]


def output_language(text: str, default: str = "text") -> str:
    """Select JSON highlighting for JSON values, otherwise use the fallback.

    :param text: Saved notebook output.
    :param default: Language for non-JSON output.
    :return: Code-fence language.
    """
    try:
        json.loads(text)
    except (ValueError, TypeError):
        return default
    return "json"


def update_supported_versions(path: Path) -> None:
    """Replace the supported-versions marker with the installed versions.

    :param path: Staged compatibility page.
    """
    table = (
        "| Standard | Supported version |\n| --- | --- |\n"
        f"| VRS | {VRS_VERSION} |\n| GKM-Core | {CORE_VERSION} |"
    )
    content = path.read_text(encoding="utf-8")
    path.write_text(
        content.replace("<!-- supported-versions -->", table), encoding="utf-8"
    )


def stage_notebook(source: Path, destination: Path, exporter: MarkdownExporter) -> None:
    """Write a notebook's rendered page, source, and saved output files.

    :param source: Notebook within the repository.
    :param destination: Staged notebook directory.
    :param exporter: Markdown exporter using saved notebook outputs.
    """
    body, resources = exporter.from_filename(str(source))
    body = body.replace("](README.md)", "](index.md)")
    # Notebook links should open the corresponding rendered page.
    body = re.sub(r"\]\(([^:)]+)\.ipynb(#[^)]*)?\)", r"](\1.md\2)", body)

    notice = (
        '!!! info "Notebook source"\n\n'
        f"    Rendered from [{source.relative_to(ROOT)}]({source.name}).\n\n"
    )
    if "Upcoming" in source.name:
        notice += (
            '!!! warning "Exploratory material"\n\n'
            "    This notebook describes upcoming features, not a guarantee of "
            "support in the current release.\n\n"
        )

    (destination / f"{source.stem}.md").write_text(notice + body, encoding="utf-8")
    copy2(source, destination / source.name)
    for name, data in resources.get("outputs", {}).items():
        target = destination / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)


def main() -> None:
    """Rebuild docs-build from source without executing notebooks."""
    output = ROOT / "docs-build"
    if output.exists():
        rmtree(output)
    copytree(ROOT / "docs", output)
    update_supported_versions(output / "compatibility.md")

    exporter = MarkdownExporter(
        extra_template_basedirs=[str(ROOT / "scripts/templates")],
        template_file=str(ROOT / "scripts/templates/notebook.md.j2"),
        filters={"output_language": output_language},
    )
    for source in sorted((ROOT / "notebooks/getting_started").glob("*.ipynb")):
        stage_notebook(source, output / "notebooks", exporter)


if __name__ == "__main__":
    main()
