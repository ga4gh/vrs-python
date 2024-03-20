import re
import hgvs
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.normalizer
import hgvs.variantmapper
import hgvs.sequencevariant

from bioutils.accessions import coerce_namespace
from dataproxy import _DataProxy

class HgvsTools():
    """
    A convenience class that exposes only the tools needed by vrs-python for working with HGVS (Human Genome Variation Society) notation.

    Attributes:
        parser (hgvs.parser.Parser): The HGVS parser object.
        uta_conn: The UTA (Unified Transcriptome and Annotation) connection object.
        normalizer: The HGVS normalizer object.
        variantMapper: The HGVS variant mapper object.
    """
    hgvs_re = re.compile(r"[^:]+:[cgmnpr]\.")

    def __init__(self,dp: _DataProxy = None):
        self.parser = hgvs.parser.Parser()
        self.uta_conn = None
        self.normalizer = None
        self.variantMapper = None
        self.dp = dp

    def _fetch_uta_conn(self):
        # just in time instantiation of the UTA connection
        if self.uta_conn is None:
            self.uta_conn = hgvs.dataproviders.uta.connect()
        return self.uta_conn
    
    def _fetch_normalizer(self):
        # just in time instantiation of the normalizer service, with validation set to True
        if self.normalizer is None:
            self.normalizer = hgvs.normalizer.Normalizer(self._fetch_uta_conn, validate=True)
        return self.normalizer
    
    def _fetch_variantMapper(self):
        # just in time instantiation of the variantMapper service
        if self.variantMapper is None:
            self.variantMapper = hgvs.variantmapper.VariantMapper(self._fetch_uta_conn)
        return self.variantMapper
    
    def close(self):
        # TODO These should only be closed if they are owned by this instance
        self.normalizer = None
        self.variantMapper = None
        self.dp = None
        if self._get_uta_conn is not None:
            self._get_uta_conn.close()

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
        if sv == None or sv.posedit == None or sv.posedit.edit == None:
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
                state = self.dp.get_sequence(sv.ac, start_i=sv.posedit.pos.start.base - 1, end_i=sv.posedit.pos.end.base)
            else:
                state = sv.posedit.edit.alt or ""

        elif sv.posedit.edit.type == "dup":
            start = sv.posedit.pos.start.base - 1
            end = sv.posedit.pos.end.base
            ref = self.dp.get_sequence(sv.ac, start_i=sv.posedit.pos.start.base - 1, end_i=sv.posedit.pos.end.base)
            state = ref + ref

        else:
            raise ValueError(f"HGVS variant type {sv.posedit.edit.type} is unsupported")

        return start, end, state
    
    # def get_hgvs_refget_accession(self, accession: str):
    #     """Get refget accession for hgvs sequence variant"""
    #     # prefix accession with namespace
    #     accession_curie = coerce_namespace(accession)
    #     return self.dp.get_refget_accession(accession_curie)

    def normalize(self, hgvs):
        return self._fetch_normalizer.normalize(hgvs)
    
    def n_to_c(self, hgvs):
        return self._fetch_variantMapper.n_to_c(hgvs)

    def c_to_n(self, hgvs):
        return self._fetch_variantMapper.c_to_n(hgvs)