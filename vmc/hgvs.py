from . import models
from .seqrepo import get_vmc_sequence_id

import hgvs
import hgvs.parser

hp = hgvs.parser.Parser()

def from_hgvs(hgvs_string):
    v = hp.parse_hgvs_variant
    
