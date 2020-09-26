import copy

from bioutils.accessions import coerce_namespace
from bioutils.assemblies import make_name_ac_map
from bioutils.cytobands import get_cytoband_maps

import ga4gh.vrs


assy_name_to_map_name = {
    "GRCh37": "ucsc-hg19",
    "GRCh38": "ucsc-hg38",
    }


class Localizer:
    """provides conversion of ChromosomeLocation objects to
    SequenceLocation objects

    """

    def __init__(self):
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


    def localize_allele(self, allele):
        # copy input variant and replace location
        # N.B. deepcopy leads to recursion errors
        allele_sl = ga4gh.vrs.models.Variation(**allele.as_dict())
        del allele_sl._id
        allele_sl.location = self.localize(allele.location)
        return allele_sl


    def localize_named_feature(self, loc, assembly_name):
        """converts named features to sequence locations

        """

        assert loc.type._value == "ChromosomeLocation", "Expected a ChromosomeLocation object"

        def _get_coords(m, cb):
            """return (start,end) of band `cb` in map `m`"""
            if cb is None:
                return None
            return chr_cb_map[cb][0:2]
            
        try:
            map_name = assy_name_to_map_name[assembly_name]
        except KeyError:
            raise KeyError(f"No cytoband map for assembly {assembly_name}")
        
        cb_map = self._cb_maps[map_name]
        
        try:
            chr_cb_map = cb_map[loc.chr]
        except KeyError:
            raise KeyError(f"{loc.chr}: Chromosome name doesn't exist in cytoband map ({assembly_name}/{map_name})")
        
        coords = []
        try:
            coords += _get_coords(chr_cb_map, loc.interval.start)
        except:
            raise ValueError(f"{loc.interval.start}: ChromosomeLocation not in map for {assembly_name}, chr {loc.chr}")
        try:
            coords += _get_coords(chr_cb_map, loc.interval.end)
        except:
            raise ValueError(f"{loc.interval.end}: ChromosomeLocation not in map for {assembly_name}, chr {loc.chr}")
 
        # the following works regardless of orientation of bands and number of bands
        start, end = min(coords), max(coords)
        
        try:
            ac = self._ana_maps[assembly_name][loc.chr]
        except KeyError:
            raise ValueError(f"No accssion for {loc.chr} in assembly {assembly_name}")

        return ga4gh.vrs.models.SequenceLocation(
            sequence_id = coerce_namespace(ac),
            interval = ga4gh.vrs.models.SimpleInterval(start=start, end=end)
            )




    
if __name__ == "__main__":
    cbl = ga4gh.vrs.models.ChromosomeLocation(chr="11", start="q22.3", end="q23.1")
    lr = Localizer()
