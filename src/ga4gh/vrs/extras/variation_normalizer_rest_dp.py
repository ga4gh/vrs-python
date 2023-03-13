"""Module for Variation Normalizer REST Data Proxy. This is needed for the notebooks,
so that users do not have to manually install and setup UTA.
"""
import requests

class VariationNormalizerRESTDataProxy:  # pylint: disable=too-few-public-methods
    """
    Rest data proxy for Variation Normalizer API
    """
    def __init__(self) -> None:
        """
        Initialize class with the API URL
        """
        self.api = "https://normalize.cancervariants.org/variation"

    def to_hgvs(self, vo, namespace="refseq"):
        """
        tranlsate vrs allele object (vo) to hgvs format
        Use this method if you don't have UTA installed locally or are unable
        to reach the public UTA database due to port issues.
        """
        vo = vo.as_dict()
        data = dict(
            variation = vo,
            namespace = namespace
        )
        r = requests.post(
            url = f"{self.api}/vrs_allele_to_hgvs",
            json=data,
            headers={ "Content-Type": "application/json; charset=utf-8" },
            timeout=5
        )
        if r.status_code == 200:  # pylint: disable=no-else-return
            return r.json().get("variations", [])
        else:
            raise requests.HTTPError(f"Variation normalizer returned the status code: {r.status_code}.")
