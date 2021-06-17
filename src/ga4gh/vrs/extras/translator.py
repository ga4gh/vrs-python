"""Translates various external formats into VRS models.

Input formats: VRS (serialized), hgvs, spdi, gnomad (vcf), beacon
Output formats: VRS (serialized), hgvs, spdi, gnomad (vcf)

"""

from collections.abc import Mapping
import copy
import logging
import re

from bioutils.accessions import coerce_namespace
import hgvs.parser

import hgvs.location
import hgvs.posedit
import hgvs.edit
import hgvs.sequencevariant
import hgvs.dataproviders.uta

from ga4gh.core import ga4gh_identify
from ga4gh.vrs import models, normalize
from ga4gh.vrs.extras.decorators import lazy_property  # this should be relocated
from ga4gh.vrs.utils.hgvs_tools import HgvsTools

from pysam import VariantFile, VariantHeader, VariantRecord
from difflib import ndiff

_logger = logging.getLogger(__name__)


class Translator:
    """Translates various variation formats to and from GA4GH VRS models

    All `from_` methods follow this pattern:
    * If the argument does not appear to be an appropriate type, None is returned
    * Otherwise, the argument is expected to be of the correct type.  If an error occurs during processing,
      an exception is raised.
    * Otherwise, the VRS object is returned

    """

    beacon_re = re.compile(r"(?P<chr>[^-]+)\s*:\s*(?P<pos>\d+)\s*(?P<ref>\w+)\s*>\s*(?P<alt>\w+)")
    gnomad_re = re.compile(r"(?P<chr>[^-]+)-(?P<pos>\d+)-(?P<ref>\w+)-(?P<alt>\w+)")
    hgvs_re = re.compile(r"[^:]+:[cgnpr]\.")
    spdi_re = re.compile(r"(?P<ac>[^:]+):(?P<pos>\d+):(?P<del_len_or_seq>\w+):(?P<ins_seq>\w+)")


    def __init__(self,
                 data_proxy,
                 default_assembly_name="GRCh38",
                 translate_sequence_identifiers=True,
                 normalize=True,
                 identify=True):
        self.default_assembly_name = default_assembly_name
        self.data_proxy = data_proxy
        self.translate_sequence_identifiers = translate_sequence_identifiers
        self.identify = identify
        self.normalize = normalize
        self.hgvs_tools = None


    def translate_from(self, var, fmt=None):
        """Translate variation `var` to VRS object

        If fmt is None, guess the appropriate format and return the variant.
        If fmt is specified, try only that format.
        See also notes about `from_` and `to_` methods.
        """

        if fmt:
            t = self.from_translators[fmt]
            o = t(self, var)
            if o is None:
                raise ValueError(f"Unable to parse data as {fmt} variation")
            return o

        for fmt, t in self.from_translators.items():
            o = t(self, var)
            if o:
                return o

        formats = list(self.from_translators.keys())
        raise ValueError(f"Unable to parse data as {', '.join(formats)}")

    def translate_to(self, vo, fmt):
        t = self.to_translators[fmt]
        return t(self, vo)


    ############################################################################
    ## INTERNAL

    def _from_beacon(self, beacon_expr, assembly_name=None):
        """Parse beacon expression into VRS Allele

        #>>> a = tlr.from_beacon("13 : 32936732 G > C")
        #>>> a.as_dict()
        {'location': {'interval': {'end': 32936732,
           'start': 32936731,
           'type': 'SimpleInterval'},
          'sequence_id': 'GRCh38:13 ',
          'type': 'SequenceLocation'},
         'state': {'sequence': 'C', 'type': 'SequenceState'},
         'type': 'Allele'}

        """

        if not isinstance(beacon_expr, str):
            return None
        m = self.beacon_re.match(beacon_expr.replace(" ", ""))
        if not m:
            return None

        g = m.groupdict()
        if assembly_name is None:
            assembly_name = self.default_assembly_name
        sequence_id = assembly_name + ":" + g["chr"]
        start = int(g["pos"]) - 1
        ref = g["ref"]
        alt = g["alt"]
        end = start + len(ref)
        ins_seq = alt

        interval = models.SimpleInterval(start=start, end=end)
        location = models.Location(sequence_id=sequence_id, interval=interval)
        sstate = models.SequenceState(sequence=ins_seq)
        allele = models.Allele(location=location, state=sstate)
        allele = self._post_process_imported_allele(allele)
        return allele


    def _from_gnomad(self, gnomad_expr, assembly_name=None):
        """Parse gnomAD-style VCF expression into VRS Allele

        #>>> a = tlr.from_gnomad("1-55516888-G-GA")
        #>>> a.as_dict()
        {'location': {'interval': {'end': 55516888,
           'start': 55516887,
           'type': 'SimpleInterval'},
          'sequence_id': 'GRCh38:1',
          'type': 'SequenceLocation'},
         'state': {'sequence': 'GA', 'type': 'SequenceState'},
         'type': 'Allele'}

        """

        if not isinstance(gnomad_expr, str):
            return None
        m = self.gnomad_re.match(gnomad_expr)
        if not m:
            return None

        g = m.groupdict()
        if assembly_name is None:
            assembly_name = self.default_assembly_name
        sequence_id = assembly_name + ":" + g["chr"]
        start = int(g["pos"]) - 1
        ref = g["ref"]
        alt = g["alt"]
        end = start + len(ref)
        ins_seq = alt

        interval = models.SimpleInterval(start=start, end=end)
        location = models.Location(sequence_id=sequence_id, interval=interval)
        sstate = models.SequenceState(sequence=ins_seq)
        allele = models.Allele(location=location, state=sstate)
        allele = self._post_process_imported_allele(allele)
        return allele


    def _from_hgvs(self, hgvs_expr):
        """parse hgvs into a VRS object (typically an Allele)

        #>>> a = tlr.from_hgvs("NM_012345.6:c.22A>T")
        #>>> a.as_dict()
        {
          'location': {
            'interval': {'end': 22, 'start': 21, 'type': 'SimpleInterval'},
            'sequence_id': 'refseq:NM_012345.6',
            'type': 'SequenceLocation'
          },
          'state': {'sequence': 'T', 'type': 'SequenceState'},
          'type': 'Allele'
        }

        """

        if not isinstance(hgvs_expr, str):
            return None
        if not self.hgvs_re.match(hgvs_expr):
            return None

        sv = self._hgvs_parser.parse_hgvs_variant(hgvs_expr)

        # prefix accession with namespace
        sequence_id = coerce_namespace(sv.ac)

        if isinstance(sv.posedit.pos, hgvs.location.BaseOffsetInterval):
            if sv.posedit.pos.start.is_intronic or sv.posedit.pos.end.is_intronic:
                raise ValueError("Intronic HGVS variants are not supported ({sv.posedit})")

        if sv.posedit.edit.type == 'ins':
            interval = models.SimpleInterval(start=sv.posedit.pos.start.base,
                                        end=sv.posedit.pos.start.base)
            state = sv.posedit.edit.alt
        elif sv.posedit.edit.type in ('sub', 'del', 'delins', 'identity'):
            interval = models.SimpleInterval(start=sv.posedit.pos.start.base - 1,
                                       end=sv.posedit.pos.end.base)
            if sv.posedit.edit.type == 'identity':
                state = self.data_proxy.get_sequence(sv.ac,
                                                     sv.posedit.pos.start.base - 1,
                                                     sv.posedit.pos.end.base)
            else:
                state = sv.posedit.edit.alt or ''
        elif sv.posedit.edit.type == 'dup':

            interval = models.SimpleInterval(start=sv.posedit.pos.start.base - 1,
                                             end=sv.posedit.pos.end.base)

            ref = self.data_proxy.get_sequence(sv.ac,
                                               sv.posedit.pos.start.base - 1,
                                               sv.posedit.pos.end.base)
            state = ref + ref
        else:
            raise ValueError(f"HGVS variant type {sv.posedit.edit.type} is unsupported")

        location = models.Location(sequence_id=sequence_id, interval=interval)
        sstate = models.SequenceState(sequence=state)
        allele = models.Allele(location=location, state=sstate)
        allele = self._post_process_imported_allele(allele)
        return allele


    def _from_spdi(self, spdi_expr):

        """Parse SPDI expression in to a GA4GH Allele

        #>>> a = tlr.from_spdi("NM_012345.6:21:1:T")
        #>>> a.as_dict()
        {
          'location': {
            'interval': {'end': 22, 'start': 21, 'type': 'SimpleInterval'},
            'sequence_id': 'refseq:NM_012345.6',
            'type': 'SequenceLocation'
          },
          'state': {'sequence': 'T', 'type': 'SequenceState'},
          'type': 'Allele'
        }
        """

        if not isinstance(spdi_expr, str):
            return None
        m = self.spdi_re.match(spdi_expr)
        if not m:
            return None

        g = m.groupdict()
        sequence_id = coerce_namespace(g["ac"])
        start = int(g["pos"])
        try:
            del_len = int(g["del_len_or_seq"])
        except ValueError:
            del_len = len(g["del_len_or_seq"])
        end = start + del_len
        ins_seq = g["ins_seq"]

        interval = models.SimpleInterval(start=start, end=end)
        location = models.Location(sequence_id=sequence_id, interval=interval)
        sstate = models.SequenceState(sequence=ins_seq)
        allele = models.Allele(location=location, state=sstate)
        allele = self._post_process_imported_allele(allele)
        return allele


    def _from_vrs(self, var):
        """convert from dict representation of VRS JSON to VRS object"""
        if not isinstance(var, Mapping):
            return None
        if "type" not in var:
            return None
        try:
            model = models[var["type"]]
        except KeyError:
            return None
        return model(**var)

    # patterns for from_VCF
    chrom_re = re.compile(r'(chr)?(?P<chrom>([0-2]?[0-9]|X|Y))', re.IGNORECASE)

    def _from_vcf_record(self, chrom, pos, ref, alts, assembly_name=None):
        """Given provided record attributes, return a VRS Allele object,
        or None if parsing fails. Can optionally pass an assembly version
        name -- if not, will try to use the class `default_assembly_name`
        value.

        Currently working:
         * Basic SNVs
         * Basic insertions
         * Basic deletions

        TODO:
         * Plan some logic branches for SVs (raise not implemented etc)
         * Handle 0th coordinate refs
         * sequence_id: ensembl
         * other types of variations
         * perform inference of reference sequence?
         * microsatellite vs copy number variation?
         * handle multiple alleles -- waiting on future versions of VRS?
        """
        # construct location
        start = int(pos) - 1
        end = start + len(ref)
        interval = models.SimpleInterval(start=start, end=end)

        if assembly_name:
            sequence_id = f'{assembly_name}:{chrom}'
        elif self.default_assembly_name:
            sequence_id = f'{self.default_assembly_name}:{chrom}'
        else:
            # TODO infer assembly by checking ref base against seqrepo?
            # raise exception?
            return None
        if not sequence_id:
            return None

        location = models.Location(sequence_id=sequence_id, interval=interval)

        # construct state
        alternate_bases = list(filter(None, alts))[0]
        seqstate = models.SequenceState(sequence=alternate_bases)

        # construct return object
        allele = models.Allele(location=location, state=seqstate)
        allele = self._post_process_imported_allele(allele)
        return allele

    def _from_vcf(self, vcf_path, assembly_name=None):
        """Given a path to a VCF file (as a str) and optionally an assembly
        name, parse the file and return a List of valid VRS Allele objects.

        TODO:
         * how to handle GT fields? (May need to wait on future VRS versions?)
         * worth trying to parse info fields to get an assembly ID?
         * enforce filter requirements?
        """
        vcf_in = VariantFile(vcf_path)
        vrs_alleles = []
        for record in vcf_in:
            vrs_allele = self._from_vcf_record(record.chrom, record.pos,
                                               record.ref, record.alts,
                                               assembly_name)
            if vrs_allele:
                vrs_alleles.append(vrs_allele)

        return vrs_alleles

    def _get_hgvs_tools(self):
        """ Only create UTA db connection if needed. There will be one connectionn per translator.
        """
        if self.hgvs_tools is None:
            uta_conn = hgvs.dataproviders.uta.connect()
            self.hgvs_tools = HgvsTools(uta_conn)
        return self.hgvs_tools

    def _to_hgvs(self, vo, namespace="refseq"):
        """generates a *list* of HGVS expressions for VRS Allele.

        If `namespace` is not None, returns HGVS strings for the
        specified namespace.

        If `namespace` is None, returns HGVS strings for all alias
        translations.

        If no alias translations are available, an empty list is
        returned.

        If the VRS object cannot be expressed as HGVS, raises ValueError.

        """

        def ir_stype(a):
            if a.startswith("refseq:NM_"):
                return "n"
            if a.startswith("refseq:NP_"):
                return "p"
            if a.startswith("refseq:NG_"):
                return "g"
            if a.startswith("refseq:NC_"):
                return "g"
            if a.startswith("GRCh"):
                return "g"
            return None

        if (type(vo).__name__ != "Allele"
            or type(vo.location).__name__ != "SequenceLocation"
            or type(vo.state).__name__ != "SequenceState"):
            raise ValueError(f"_to_hgvs requires a VRS Allele with SequenceLocation and SequenceState")

        sequence_id = str(vo.location.sequence_id)
        aliases = self.data_proxy.translate_sequence_identifier(sequence_id, namespace)

        # infer type of sequence based on accession
        # TODO: move to bioutils
        stypes = list(set(t for t in (ir_stype(a) for a in aliases) if t))
        if len(stypes) != 1:
            raise ValueError(f"Couldn't infer sequence type for {sequence_id} ({stypes})")
        stype = stypes[0]

        # build interval and edit depending on sequence type
        if stype == "p":
            raise ValueError("Only nucleic acid variation is currently supported")
            # ival = hgvs.location.Interval(start=start, end=end)
            # edit = hgvs.edit.AARefAlt(ref=None, alt=vo.state.sequence)
        else:
            start = vo.location.interval.start
            end = vo.location.interval.end
            # ib: 0 1 2 3 4 5
            #  h:  1 2 3 4 5
            if start == end:    # insert: hgvs uses *exclusive coords*
                ref = None
                end += 1
            else:               # else: hgvs uses *inclusive coords*
                ref = self.data_proxy.get_sequence(sequence_id, start, end)
                start += 1
            ival = hgvs.location.Interval(
                start=hgvs.location.SimplePosition(base=start),
                end=hgvs.location.SimplePosition(base=end))
            alt = str(vo.state.sequence) or None  # "" => None
            edit = hgvs.edit.NARefAlt(ref=ref, alt=alt)

        posedit = hgvs.posedit.PosEdit(pos=ival, edit=edit)
        var = hgvs.sequencevariant.SequenceVariant(
            ac=None,
            type=stype,
            posedit=posedit)

        hgvs_exprs = []
        for alias in aliases:
            ns, a = alias.split(":")
            # skip GRCh accessions unless specifically requested
            # because they are ambiguous without their namespace,
            # which can't be included in HGVS expressions
            # TODO: use default_assembly_name here
            if ns.startswith("GRC") and namespace is None:
                continue
            var.ac = a

            if not namespace.startswith("GRC"):
                # if the namespace is GRC, can't normalize, since hgvs can't deal with it
                hgvs_tools = self._get_hgvs_tools()
                parsed = hgvs_tools.parse(str(var))
                var = hgvs_tools.normalize(parsed)

            hgvs_exprs += [str(var)]

        return list(set(hgvs_exprs))

    def _to_spdi(self, vo, namespace="refseq"):
        """generates a *list* of SPDI expressions for VRS Allele.

        If `namespace` is not None, returns SPDI strings for the
        specified namespace.

        If `namespace` is None, returns SPDI strings for all alias
        translations.

        If no alias translations are available, an empty list is
        returned.

        If the VRS object cannot be expressed as SPDI, raises ValueError.

        SPDI and VRS use identical normalization. The incoming Allele
        is expected to be normalized per VRS spec.

        """

        if (type(vo).__name__ != "Allele"
            or type(vo.location).__name__ != "SequenceLocation"
            or type(vo.state).__name__ != "SequenceState"):
            raise ValueError(f"_to_hgvs requires a VRS Allele with SequenceLocation and SequenceState")

        sequence_id = str(vo.location.sequence_id)
        aliases = self.data_proxy.translate_sequence_identifier(sequence_id, namespace)
        aliases = [a.split(":")[1] for a in aliases]

        start = vo.location.interval.start
        end = vo.location.interval.end
        spdi_tail = f":{start}:{end-start}:{vo.state.sequence}"
        spdis = [a + spdi_tail for a in aliases]
        return spdis

    def _allele_to_vcf(self, vcfh, vo, namespace="refseq"):
        """Given a Pysam VariantHeader object `vcfh`, and a VRS Allele object
        `vo`, return a Pysam VariantRecord

        TODO
         * more elegant way of getting chrom
        """
        if (type(vo).__name__ != "Allele"
            or type(vo.location).__name__ != "SequenceLocation"
            or type(vo.state).__name__ != "SequenceState"):
            raise ValueError(f"_to_vcf requires a VRS Allele with SequenceLocation and SequenceState")

        start = int(vo.location.interval.start)
        end = int(vo.location.interval.end)
        alt = str(vo.state.sequence)
        alt_sequence_length = len(alt)
        interval_length = end - start
        sequence_id = str(vo.location.sequence_id)
        ref_sequence = self.data_proxy.get_sequence(sequence_id, start, end)
        aliases = self.data_proxy.translate_sequence_identifier(sequence_id,
                                                                namespace)

        if interval_length == alt_sequence_length:
            # SNVs/MNVs
            pos = int(start)
            ref = ref_sequence
        elif interval_length > alt_sequence_length:
            # del
            length_diff = interval_length - alt_sequence_length
            ref = self.data_proxy.get_sequence(sequence_id, start - 1, start + length_diff)
            alt = self.data_proxy.get_sequence(sequence_id, start - 1, start)
            pos = start - 1
        else:
            # ins
            diff_bases = [i[2] for i in ndiff(ref_sequence[::-1], alt[::-1])
                          if i.startswith('+ ')][::-1]
            diff = ''.join(diff_bases)
            pos = start - 1
            ref = self.data_proxy.get_sequence(sequence_id, start - 1, start)
            alt = ref + diff

        if namespace == 'refseq':
            chrom = str(int(aliases[0][14:16]))  # drop leading 0
        else:
            raise NotImplementedError  # TODO
        if chrom not in vcfh.contigs.keys():
            vcfh.add_meta('contig', items=[('ID', chrom)])

        record = vcfh.new_record(contig=chrom, start=pos,
                                 alleles=(ref, alt))
        return record

    def _to_vcf(self, vrs_objects, file_path, namespace="refseq"):
        """Given an iterable collection of VRS Alleles, write to `file_path`
        a basic VCF file.

        Raises ValueError if provided list contains non-Allele objects, or
        if they don't have SequenceLocation and SequenceState attributes.

        WORKING
         * basic SNVs, insertions, deletions

        TODO
         * be more intelligent about getting namespace IDs
         * other variation types
        """
        vcfh = VariantHeader()
        records = []

        for vo in vrs_objects:
            records.append(self._allele_to_vcf(vcfh, vo, namespace))
        records.sort(key=lambda r: ((r.chrom).zfill(2), r.pos))

        vcf = VariantFile(file_path, 'w', header=vcfh)
        for record in records:
            vcf.write(record)
        vcf.close()


    @lazy_property
    def _hgvs_parser(self):
        """instantiates and returns an hgvs parser instance"""
        _logger.info("Creating  parser")
        return hgvs.parser.Parser()


    def _post_process_imported_allele(self, allele):
        """Provide common post-processing for imported Alleles IN-PLACE.

        """

        if self.translate_sequence_identifiers:
            seq_id = self.data_proxy.translate_sequence_identifier(allele.location.sequence_id._value, "ga4gh")[0]
            allele.location.sequence_id = seq_id

        if self.normalize:
            allele = normalize(allele, self.data_proxy)

        if self.identify:
            allele._id = ga4gh_identify(allele)

        return allele


    def _seq_id_mapper(self, ir):
        if self.translate_sequence_identifiers:
            return self.data_proxy.translate_sequence_identifier(ir, "ga4gh")[0]
        return ir


    from_translators = {
        "beacon": _from_beacon,
        "gnomad": _from_gnomad,
        "hgvs": _from_hgvs,
        "spdi": _from_spdi,
        "vrs": _from_vrs,
        "vcf": _from_vcf,
    }

    to_translators = {
        "hgvs": _to_hgvs,
        "spdi": _to_spdi,
        #"gnomad": to_gnomad,
        "vcf": _to_vcf,
    }


if __name__ == "__main__":
    import coloredlogs
    coloredlogs.install(level="INFO")

    from ga4gh.vrs.dataproxy import create_dataproxy
    dp = create_dataproxy("seqrepo+file:///usr/local/share/seqrepo/latest")
    tlr = Translator(data_proxy=dp)

    expressions = [
        "bogus",
        "1-55516888-G-GA",
        "13 : 32936732 G > C",
        "NC_000013.11:g.32936732G>C",
        "NM_000551.3:21:1:T", {
            'location': {
                'interval': {
                    'end': 22,
                    'start': 21,
                    'type': 'SimpleInterval'
                },
                'sequence_id': 'ga4gh:SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_',
                'type': 'SequenceLocation'
            },
            'state': {
                'sequence': 'T',
                'type': 'SequenceState'
            },
            'type': 'Allele'
        }, {
            'end': 22,
            'start': 21,
            'type': 'SimpleInterval'
        }
    ]
    formats = ["hgvs", "gnomad", "beacon", "spdi", "vrs", None]

    for e in expressions:
        print(f"* {e}")
        for f in formats:
            try:
                o = tlr.translate_from(e, f)
                r = o.type
            except ValueError:
                r = "-"
            except Exception as ex:
                r = ex.__class__.__name__
            print(f"  {f}: {r}")
