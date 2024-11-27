from setuptools import setup
from pathlib import Path
import subprocess


class VrsSubmoduleFetchError(Exception):
    """Raise for errors during submodule metadata extraction"""


def get_vrs_submodule_info():
    """Retrieve the commit hash and tag from the vrs submodule."""
    try:
        return subprocess.check_output(
            ["git", "describe", "--tags", "--abbrev=0"], cwd="submodules/vrs", text=True
        ).strip()
    except Exception as e:
        raise VrsSubmoduleFetchError from e


def write_metadata_file():
    """Generate a Python module with submodule metadata."""
    tag = get_vrs_submodule_info()
    metadata_path = Path("src/ga4gh/vrs/submodule_metadata.py")
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(
        f"# Auto-generated metadata from VRS schema submodule\n"
        f"VRS_RELEASE = '{tag}'\n"
    )


write_metadata_file()

# Use setuptools for the actual build (delegates to pyproject.toml)
setup()
