import hgvs
import hgvs.parser

from .models import Allele



hp = hgvs.parser.Parser()



def from_hgvs(hgvs_string):
    v = hp.parse_hgvs_variant
