#!/usr/bin/env python

import os
import pprint

import icu
import tqdm

import vmc.digest

def sorted_strings(strings, locale=None):
    if locale is None:       
        return sorted(strings)                                 
    collator = icu.Collator.createInstance(icu.Locale(locale))
    return sorted(strings, key=collator.getSortKey)        


locales = sorted([l for l in icu.Locale.getAvailableLocales().keys() if "_" in l])

binaries = sorted([os.urandom(50) for _ in range(1000)])
digests = vmc.digest.vmc_digest(os.urandom(100))


for l in tqdm.tqdm(locales):
    sdigests = sorted_strings(digests, locale=l)
    if digests != sdigests:
        pprint.pprint({"digests": digests, "sdigests": sdigests})
        raise Exception("failed on locale " + l)
