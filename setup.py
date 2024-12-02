from setuptools import setup
from pathlib import Path
import subprocess


class VrsSubmoduleFetchError(Exception):
    """Raise for errors during submodule metadata extraction"""


def _get_vrs_submodule_info() -> tuple[str, str, int]:
    """Retrieve the commit hash and tag from the vrs submodule.

    :return: tuple containing commit hash, latest version tag, and rev distance between them
    """
    try:
        vrs_path = "submodules/vrs"
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=vrs_path, text=True
        ).strip()
        subprocess.run(["git", "fetch", "--all", "--tags"], cwd=vrs_path)
        tag = subprocess.check_output(
            ["git", "describe", "--tags", "--abbrev=0"], cwd=vrs_path, text=True
        ).strip()
        distance =  subprocess.check_output(
            ["git", "rev-list", f"{tag}..HEAD", "--count"],
            cwd=vrs_path,
            text=True,
        ).strip()
        return commit, tag, int(distance)
    except Exception as e:
        raise VrsSubmoduleFetchError from e


def write_metadata_file():
    """Generate a Python module with submodule metadata."""
    commit, tag, distance = _get_vrs_submodule_info()
    metadata_path = Path("src/ga4gh/vrs/submodule_metadata.py")
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    if distance == 0:
        version = tag
    else:
        version = f"{tag}.dev{distance}+g{commit}"
    metadata_path.write_text(
        f"# Auto-generated metadata from VRS schema submodule\n"
        f"VRS_VERSION = '{version}'\n"
    )


write_metadata_file()

# Use setuptools for the actual build (delegates to pyproject.toml)
setup()
