import re
import hgvs
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.normalizer
import hgvs.variantmapper
import hgvs.sequencevariant

import logging

from ga4gh.vrs.dataproxy import _DataProxy, create_dataproxy
from ga4gh.vrs import models

_logger = logging.getLogger(__name__)

class HgvsTools():
    """
    A convenience class that exposes only the tools needed by vrs-python for working with HGVS (Human Genome Variation Society) notation.

    Attributes:
        parser (hgvs.parser.Parser): The HGVS parser object.
        uta_conn: The UTA (Universal Transcript Archive) connection object.
        normalizer: The HGVS normalizer object.
        variant_mapper: The HGVS variant mapper object.
        data_proxy (_DataProxy): The data proxy object.
    """
    hgvs_re = re.compile(r"[^:]+:[cgmnpr]\.")

    def __init__(self,data_proxy: _DataProxy = None):
        self.parser = hgvs.parser.Parser()
        self.uta_conn = hgvs.dataproviders.uta.connect()
        self.normalizer = hgvs.normalizer.Normalizer(self.uta_conn, validate=True)
        self.variant_mapper = hgvs.variantmapper.VariantMapper(self.uta_conn)
        self.data_proxy = data_proxy

    

    def close(self):
        # TODO These should only be closed if they are owned by this instance
        self.normalizer = None
        self.variant_mapper = None
        self.data_proxy = None
        if self.uta_conn is not None:
            self.uta_conn.close()

    # convenience methods for hgvs parsing, normalization, and some mappings
    def parse(self, hgvs_str):
        """
        Parses the given HGVS string and returns the corresponding variant.

        Args:
            hgvs_str (str): The HGVS string to parse.

        Returns:
            Variant: The parsed variant object, or None if the HGVS string is invalid.
        """
        if not self.hgvs_re.match(hgvs_str):
            return None
        return self.parser.parse_hgvs_variant(hgvs_str)
    
    def is_intronic(self, sv: hgvs.sequencevariant.SequenceVariant):
        """
        Checks if the given SequenceVariant is intronic.

        Args:
            sv (hgvs.sequencevariant.SequenceVariant): The SequenceVariant to check.

        Returns:
            bool: True if the SequenceVariant is intronic, False otherwise.
        """
        if isinstance(sv.posedit.pos, hgvs.location.BaseOffsetInterval):
            return (sv.posedit.pos.start.is_intronic or sv.posedit.pos.end.is_intronic)
        return False
    

    def get_edit_type(self, sv: hgvs.sequencevariant.SequenceVariant):
        if sv is None or sv.posedit is None or sv.posedit.edit is None:
            return None
        return sv.posedit.edit.type
    
    def get_position_and_state(self, sv: hgvs.sequencevariant.SequenceVariant):
        """
        Get the details of a sequence variant.

        Args:
            sv (hgvs.sequencevariant.SequenceVariant): The sequence variant object.

        Returns:
            tuple: A tuple containing the start position, end position, and state of the variant.

        Raises:
            ValueError: If the HGVS variant type is unsupported.
        """

        if sv.posedit.edit.type == "ins":
            start = sv.posedit.pos.start.base
            end = sv.posedit.pos.start.base
            state = sv.posedit.edit.alt

        elif sv.posedit.edit.type in ("sub", "del", "delins", "identity"):
            start = sv.posedit.pos.start.base - 1
            end = sv.posedit.pos.end.base
            if sv.posedit.edit.type == "identity":
                state = self.data_proxy.get_sequence(sv.ac, start=sv.posedit.pos.start.base - 1, end=sv.posedit.pos.end.base)
            else:
                state = sv.posedit.edit.alt or ""

        elif sv.posedit.edit.type == "dup":
            start = sv.posedit.pos.start.base - 1
            end = sv.posedit.pos.end.base
            ref = self.data_proxy.get_sequence(sv.ac, start=sv.posedit.pos.start.base - 1, end=sv.posedit.pos.end.base)
            state = ref + ref

        else:
            raise ValueError(f"HGVS variant type {sv.posedit.edit.type} is unsupported")

        return start, end, state
    def extract_allele_values(self, hgvs_expr: str):
        """parse hgvs into a VRS Allele Object

        kwargs:
            rle_seq_limit Optional(int): If RLE is set as the new state after
                normalization, this sets the limit for the length of the `sequence`.
                To exclude `sequence` from the response, set to 0.
                For no limit, set to `None`.
                Defaults value set in instance variable, `rle_seq_limit`.
            do_normalize (bool): `True` if fully justified normalization should be
                performed. `False` otherwise. Defaults to `True`

        #>>> a = tlr.from_hgvs("NC_000007.14:g.55181320A>T")
        #>>> a.model_dump()
        {
          'location': {
            'end': 55181320,
            'start': 55181319,
            'sequenceReference': {
              'type': 'SequenceReference',
              'refgetAccession': 'SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul'
            },
            'type': 'SequenceLocation'
          },
          'state': {
            'sequence': 'T',
            'type': 'LiteralSequenceExpression'
          },
          'type': 'Allele'
        }

        """
        sv = self.parse(hgvs_expr)
        if not sv:
            return None
        
        if self.is_intronic(sv):
            raise ValueError("Intronic HGVS variants are not supported")

        refget_accession = self.data_proxy.derive_refget_accession(sv.ac)
        if not refget_accession:
            return None
        
        # translate coding coordinates to positional coordinates, if necessary
        if sv.type == "c":
            sv = self.c_to_n(sv)

        (start,end,state) = self.get_position_and_state(sv)

        retval = {"refget_accession": refget_accession, "start": start, "end": end, "literal_sequence": state}
        return retval
    
    def from_allele(self, vo, namespace=None):
        """generates a *list* of HGVS expressions for VRS Allele.

        If `namespace` is not None, returns HGVS strings for the
        specified namespace.

        If `namespace` is None, returns HGVS strings for all alias
        translations.

        If no alias translations are available, an empty list is
        returned.

        If the VRS object cannot be expressed as HGVS, raises ValueError.

        This method assumes that IRIs are dereferenced, providing a `SequenceReference`
        as the `vo.location.sequenceReference`. If a `SequenceReference` is not
        provided, raises TypeError
        """

        if vo is None:
            return []
        if not isinstance(vo, models.Allele): 
            raise ValueError("VRS object must be an Allele")
        if vo.location is None:
            raise ValueError("VRS allele must have a location")
        
        refget_accession = vo.location.get_refget_accession()
        if refget_accession is None:
            raise ValueError("VRS allele location must have a sequence reference")
        
        sequence = f"ga4gh:{refget_accession}"
        aliases = self.data_proxy.translate_sequence_identifier(sequence, namespace)

        hgvs_exprs = []
        for alias in aliases:
            ns, accession = alias.split(":")
            # skip GRCh accessions unless specifically requested
            # because they are ambiguous without their namespace,
            # which can't be included in HGVS expressions
            # TODO: use default_assembly_name here
            if ns.startswith("GRC") and namespace is None:
                continue

            if not (any(accession.startswith(pfx) for pfx in ("NM", "NP", "NC", "NG", "NR", "NW", "NT", "XM", "XR", "XP"))):
                continue

            sequence_type = self.data_proxy.extract_sequence_type(alias)

            # create the hgvs expression object
            var = self._to_sequence_variant(vo, sequence_type, sequence, accession)
            hgvs_exprs += [str(var)]
         
        return list(set(hgvs_exprs))
    
    def _to_sequence_variant(self, vo, sequence_type, sequence, accession):
        """Creates a SequenceVariant object from an Allele object."""
          # build interval and edit depending on sequence type
        if sequence_type == "p":
            raise ValueError("Only nucleic acid variation is currently supported")
            # ival = hgvs.location.Interval(start=start, end=end)
            # edit = hgvs.edit.AARefAlt(ref=None, alt=vo.state.sequence)
        else:                   # pylint: disable=no-else-raise
            start, end = vo.location.start, vo.location.end
            # ib: 0 1 2 3 4 5
            #  h:  1 2 3 4 5
            if start == end:    # insert: hgvs uses *exclusive coords*
                ref = None
                end += 1
            else:               # else: hgvs uses *inclusive coords*
                ref = self.data_proxy.get_sequence(sequence, start, end)
                start += 1

            ival = hgvs.location.Interval(
                start=hgvs.location.SimplePosition(start),
                end=hgvs.location.SimplePosition(end)
            )
            alt = str(vo.state.sequence.root) or None  # "" => None
            edit = hgvs.edit.NARefAlt(ref=ref, alt=alt)

        posedit = hgvs.posedit.PosEdit(pos=ival, edit=edit)
        var = hgvs.sequencevariant.SequenceVariant(
            ac=accession,
            # at this point, use `n.` because the positions are absolute (not CDS),
            # this will subsequently be converted back to `c.` after hgvs normalization
            type='n' if sequence_type == 'c' else sequence_type,
            posedit=posedit)
        
        try:
            # if the namespace is GRC, can't normalize, since hgvs can't deal with it
            parsed = self.parse(str(var))
            var = self.normalize(parsed)
            
            # if sequence_type is coding, convert from "n." to "c." before continuing
            if sequence_type == "c":
                var = self.n_to_c(var)
            
        except hgvs.exceptions.HGVSDataNotAvailableError:
            _logger.warning(f"No data found for accession {accession}")
        
        return var

    def normalize(self, hgvs):
        return self.normalizer.normalize(hgvs)
    
    def n_to_c(self, hgvs):
        return self.variant_mapper.n_to_c(hgvs)

    def c_to_n(self, hgvs):
        return self.variant_mapper.c_to_n(hgvs)

if __name__ == "__main__":

    import os
    import json
    seqrepo_uri = os.environ.get("SEQREPO_URI", "seqrepo+file:///usr/local/share/seqrepo/latest")
    if seqrepo_uri is None:
        raise ValueError("SEQREPO_URI environment variable must be set to a valid seqrepo URI")

    allele_dict = {
        "id": "ga4gh:VA.SmhjExyBS8GxicngJxQLJ8Ww5GrBQk40",
        "type": "Allele",
        "digest": "SmhjExyBS8GxicngJxQLJ8Ww5GrBQk40",
        "location": {
            "id": "ga4gh:SL.HyEEyTvdJrfsElsTXRG-43Y_tID4QYFt",
            "type": "SequenceLocation",
            "digest": "HyEEyTvdJrfsElsTXRG-43Y_tID4QYFt",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.gg9z_DUcgwhfRz29BS3zq2jim59tOrA9"
            },
            "start": 874,
            "end": 875
        },
        "state": {
            "type": "LiteralSequenceExpression",
            "sequence": "G"
        }
    }

    vrs_allele = models.Allele.parse_obj(allele_dict)
    dp = create_dataproxy(seqrepo_uri)
    hgvsTools = HgvsTools(dp)
    hgvs_expr = hgvsTools.from_allele(vrs_allele, namespace="refseq")

    print(hgvs_expr)

