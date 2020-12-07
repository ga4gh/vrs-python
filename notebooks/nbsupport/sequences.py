from biocommons.seqrepo import SeqRepo
from ga4gh.vrs.dataproxy import SeqRepoDataProxy



_dp = SeqRepoDataProxy(SeqRepo(root_dir="/usr/local/share/seqrepo/latest"))


def translate_sequence_identifier(ir):
    try:
        return _dp.translate_sequence_identifier(ir, "ga4gh")[0]
    except IndexError:
        raise KeyError(f"Unable to translate {ir} to ga4gh sequence identifier")
