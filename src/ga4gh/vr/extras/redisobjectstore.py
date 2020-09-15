import bz2
import collections.abc
import json

class RedisObjectStore(collections.abc.MutableMapping):
    """Provides Redis-backed storage of VR objects
    
    The intention of this class is to provide a interface that is
    indistinguishable from a dictionary for the purposes of GA4GH
    objects.  In particular, that means that the the value may be a VR
    object.

    The underlying storage is compressed(utf8-encoded(json-serialization))

    """

    def __init__(self, redis, models):
        self._enc = "utf-8"
        self.redis = redis
        self.models = models


    def __contains__(self, name):
        name = name if isinstance(name, str) else str(name)
        return self.redis.exists(name)
    

    def __delitem__(self, name):
        name = name if isinstance(name, str) else str(name)
        self.redis.delete(name)


    def __getitem__(self, name):
        """gets item in k-v storage

        e.g., os[id].type == "Allele"

        """
        name = name if isinstance(name, str) else str(name)
        d = json.loads(bz2.decompress(self.redis.get(name)).decode(self._enc))
        return self.models[d["type"]](**d)


    def __iter__(self):
        yield from (bk.decode(self._enc) for bk in self.redis.keys())


    def __len__(self):
        return self.redis.dbsize()


    def __setitem__(self, name, o):
        """sets item in k-v storage

        e.g., os[ga4gh_identify(allele)] = allele

        """
        name = name if isinstance(name, str) else str(name)
        self.redis.set(name, bz2.compress(json.dumps(o.as_dict()).encode(self._enc)))






if __name__ == "__main__":
    import redis
    from ga4gh.core import ga4gh_identify
    from ga4gh.vr import models


    ros = RedisObjectStore(redis=redis.Redis(db=15), models=models)

    a0_dict = {
        'location': {
            'interval': {'end': 44908822, 'start': 44908821, 'type': 'SimpleInterval'},
            'sequence_id': 'ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl',
            'type': 'SequenceLocation'
            },
        'state': {
            'sequence': 'T',
            'type': 'SequenceState'
            },
        'type': 'Allele'
        }
    a0 = models.Allele(**a0_dict)

    a0id = ga4gh_identify(a0)
    assert a0id == "ga4gh:VA.EgHPXXhULTwoP4-ACfs-YCXaeUQJBjH_"
    a0._id = a0id

    a0l = a0.location
    a0lid = ga4gh_identify(a0l)
    assert a0lid == "ga4gh:VSL.u5fspwVbQ79QkX6GHLF8tXPCAXFJqRPx"


    ros.redis.flushdb()         # Danger!

    # db looks empty?
    assert len(ros) == 0
    assert list(ros.keys()) == []

    # add Allele
    ros[a0._id] = a0            # N.B. a0._id is a CURIE type, auto-coerced to str
    assert len(ros) == 1
    assert list(ros.keys()) == [a0id]
    assert a0id in ros

    # add Location
    ros[a0lid] = a0l
    assert len(ros) == 2
    assert sorted(list(ros.keys())) == sorted([a0id, a0lid])
    assert a0lid in ros

    # get Allele by id
    a0b = ros[a0id]
    assert a0 is not a0b
    assert a0 == a0b
    
    # get Location by id
    a0lb = ros[a0lid]
    assert a0l is not a0lb
    assert a0l == a0lb

    # delete objects
    del ros[a0id]
    assert len(ros) == 1
    del ros[a0lid]
    assert len(ros) == 0
