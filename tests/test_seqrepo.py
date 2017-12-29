from vmc import models, get_vmc_sequence_id
ir = models.Identifier(namespace="NCBI", accession="NC_000019.10")

assert 'VMC:GS_IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl' == get_vmc_sequence_id(ir)
