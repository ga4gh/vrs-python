"""convert named locations into SequenceLocations ("localize") by
reference to external data
"""
from bioutils.assemblies import make_name_ac_map
from bioutils.cytobands import get_cytoband_maps

from ga4gh.vrs import models
from ga4gh.vrs.dataproxy import _DataProxy


assy_name_to_map_name = {
    "GRCh37": "ucsc-hg19",
    "GRCh38": "ucsc-hg38",
}


class Localizer:
    """provides conversion of cytoband objects to SequenceLocation objects"""

    def __init__(self, data_proxy: _DataProxy) -> None:
        self.data_proxy = data_proxy
        # _cb_maps: maps bands to chromosomal locations
        # {'ucsc-hg19': {'1': {'p11.1': [121500000, 125000000, 'acen'],
        #              'p11.2': [120600000, 121500000, 'gneg'],
        #              'p12': [117800000, 120600000, 'gpos50'],
        #              'p13.1': [116100000, 117800000, 'gneg'],
        #              'p13.2': [111800000, 116100000, 'gpos50'],

        self._cb_maps = get_cytoband_maps()

        # _ana_maps: maps names to accessions for each assembly
        # {
        #     'GRCh37': {
        #         '1': 'NC_000001.10',
        #         '2': 'NC_000002.11',
        #         '3': 'NC_000003.11',
        #     ...
        # }
        self._ana_maps = {k: make_name_ac_map(k) for k in assy_name_to_map_name}

    def localize_named_feature(
        self, chromosome: str, start: str, end: str, assembly_name: str = "GRCh38"
    ) -> models.SequenceLocation:
        """converts named features to sequence locations

        :param chromosome: Chromosome
        :param start: Start cytoband location
        :param end: End cytoband location
        :param assembly: Assembly name
        :raise KeyError: If assembly name has no mapped cytoband or chromosome doesn't
            exit in cytoband
        :raise ValueError: If ``start``, ``end``, or ``chromosome`` not in maps
        :return: Sequence Location representation
        """
        try:
            map_name = assy_name_to_map_name[assembly_name]
        except KeyError as e:
            raise KeyError(f"No cytoband map for assembly {assembly_name}") from e

        cb_map = self._cb_maps[map_name]

        try:
            chr_cb_map = cb_map[chromosome]
        except KeyError as e:
            raise KeyError(f"{chromosome}: Chromosome name doesn't exist in cytoband"
                           f" map ({assembly_name}/{map_name})") from e

        coords = []
        for cb in (start, end):
            try:
                coords += chr_cb_map[cb][0:2]
            except (KeyError, ValueError) as e:
                err_msg = f"{cb} not in map for {assembly_name}, chromosome {chromosome}"
                raise ValueError(err_msg) from e

        try:
            ac = self._ana_maps[assembly_name][chromosome]
        except KeyError as e:
            raise ValueError(f"No accession for {chromosome} in assembly {assembly_name}") from e

        return models.SequenceLocation(
            sequenceReference=models.SequenceReference(
                refgetAccession=self.data_proxy.derive_refget_accession(ac)
            ),
            start=min(coords),  # works regardless of orientation of bands and # of bands
            end=max(coords)  # works regardless of orientation of bands and # of bands
        )
