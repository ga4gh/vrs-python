import copy

from bioutils.assemblies import make_name_ac_map
from bioutils.cytobands import get_cytoband_maps

import ga4gh.vr
from ga4gh.vr.utils import coerce_namespace


assy_name_to_map_name = {
    "GRCh37": "ucsc-hg37",
    "GRCh38": "ucsc-hg38",
    }


class Localizer:
    """convert Location object to a SequenceLocation

    """

    def __init__(self, default_assembly_name="GRCh38"):
        self.default_assembly_name = default_assembly_name
        self._cb_maps = get_cytoband_maps()
        self._ana_maps = {k: make_name_ac_map(k) for k in assy_name_to_map_name}
        
    def _localize_cytoband(self, loc, assembly_name=None):
        assert loc.type._value == "CytobandLocation", "Expected a CytobandLocation object"

        def _get_coords(m, cb):
            """return (start,end) of band `cb` in map `m`"""
            if cb is None:
                return None
            return chr_cb_map[cb][0:2]
            
        if assembly_name is None:
            assembly_name = self.default_assembly_name

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
            coords += _get_coords(chr_cb_map, loc.start)
        except:
            raise ValueError(f"{loc.start}: CytobandLocation not in map for {assembly_name}, chr {loc.chr}")
        
        if loc.end is not None:
            try:
                coords += _get_coords(chr_cb_map, loc.end)
            except:
                raise ValueError(f"{loc.end}: CytobandLocation not in map for {assembly_name}, chr {loc.chr}")
 
        # the following works regardless of orientation of bands and number of bands
        start, end = min(coords), max(coords)
        
        try:
            ac = self._ana_maps[assembly_name][loc.chr]
        except KeyError:
            raise ValueError(f"No accssion for {loc.chr} in assembly {assembly_name}")

        return ga4gh.vr.models.SequenceLocation(
            sequence_id = coerce_namespace(ac),
            region = ga4gh.vr.models.SimpleInterval(start=start, end=end)
            )


    def localize(self, v):
        loc = v.location
        typ = loc.type._value

        if typ == "CytobandLocation":
            sloc = self._localize_cytoband(loc)
        else:
            raise ValueError(f"Cannot localized variation with location type {typ}")

        # copy input variant and replace location
        # using as_dict() to copy because deepcopy led to recursion errors
        v2 = ga4gh.vr.models.Variation(**v.as_dict())
        v2.id = None
        v2.location = sloc
        return v2


    
if __name__ == "__main__":
    from ga4gh.vr.extras.localizer import Localizer
    cbl = ga4gh.vr.models.CytobandLocation(chr="11", start="q22.3", end="q23.1")
    cnvstate = ga4gh.vr.models.CNVState(min_copies=3, max_copies=5, copy_measure="RELATIVE")
    a = ga4gh.vr.models.Allele(location=cbl, state=cnvstate)
    a.id = ga4gh.vr.computed_id(a)

    lczr = Localizer()
    
