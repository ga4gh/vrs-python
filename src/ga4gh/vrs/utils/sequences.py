"""Provide generic utilities for sequence-related tasks."""

_SEQUENCE_TYPE_PREFIX_MAP = {
    "refseq:NM_": "c",
    "refseq:NC_012920": "m",
    "refseq:NG_": "g",
    "refseq:NC_00": "g",
    "refseq:NW_": "g",
    "refseq:NT_": "g",
    "refseq:NR_": "n",
    "refseq:NP_": "p",
    "refseq:XM_": "c",
    "refseq:XR_": "n",
    "refseq:XP_": "p",
    "GRCh": "g",
}


def extract_sequence_type(alias: str) -> str | None:
    """Provide a convenient way to extract the sequence type from an accession by matching its prefix to a known set of prefixes.

    Args:
        alias (str): The accession string.

    Returns:
        str or None: The sequence type associated with the accession string, or None if no matching prefix is found.

    """
    for prefix, seq_type in _SEQUENCE_TYPE_PREFIX_MAP.items():
        if alias.startswith(prefix):
            return seq_type
    return None
