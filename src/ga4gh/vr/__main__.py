"""simple command-line shell for experimenting with VR

$ ipython -m ga4gh.vr -i

"""

from ga4gh.core import *
from ga4gh.vr import *

from ga4gh.vr.dataproxy import SeqRepoRESTDataProxy
from ga4gh.vr.extras.translator import Translator


seqrepo_rest_service_url = "http://localhost:5000/seqrepo"
hgvs_expr = "NC_000013.11:g.32936732G>C"


dp = SeqRepoRESTDataProxy(base_url=seqrepo_rest_service_url)
tlr = Translator(data_proxy=dp)

allele = tlr.from_hgvs(hgvs_expr)
location = allele.location
interval = location.interval
state = allele.state

location._id = ga4gh_identify(location)
allele._id = ga4gh_identify(allele)

assert allele._id == "ga4gh:VA.n9ax-9x6gOC0OEt73VMYqCBfqfxG1XUH"
assert allele.location._id == "ga4gh:VSL.v9K0mcjQVugxTDIcdi7GBJ_R6fZ1lsYq"
assert allele.location.sequence_id == "ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT"
