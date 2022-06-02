import requests

class VariationNormalizerRESTDataProxy:
    """
    Rest data proxy for Variation Normalizer API
    """
    def __init__(self) -> None:
        """
        Initialize class with the API URL
        """
        self.api = "https://normalize.cancervariants.org/variation"

    def to_hgvs(self, vo):
        """
        tranlsate vrs allele object (vo) to hgvs format
        Use this method if you don't have UTA installed locally or are unable
        to reach the public UTA database due to port issues.
        """
        vo = vo.as_dict()
        data = dict(
            variation = vo,
            fmt = "hgvs"
        )
        r = requests.post(
            url = f"{self.api}/translate_to",
            json=data,
            headers={ "Content-Type": "application/json; charset=utf-8" }
        )
        if r.status_code == 200:
            return r.json().get("variations", [])
        else:
            raise requests.HTTPError(f"Variation normalizer returned the status code: {r.status_code}.")

