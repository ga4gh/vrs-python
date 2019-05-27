# -*- coding: utf-8 -*-
"""normalize an allele on a reference sequence

A normalized Allele is:
* Minimum length: the Allele cannot be written as a shorter Allele
  that generates the same resulting sequence; and,
* Right-aligned: the Allele cannot be written as another Allele with
  numerically greater coordinates that would generate the same
  resulting sequence.

"""

debug = False


def trim_left(alleles):
    """remove common prefix from left of all alleles, returning
    (number_trimmed, [new_alleles])

    >>> trim_left(["","AA"])
    (0, ['', 'AA'])
    >>> trim_left(["A","AA"])
    (1, ['', 'A'])
    >>> trim_left(["AT","AA"])
    (1, ['T', 'A'])
    >>> trim_left(["AA","AA"])
    (2, ['', ''])
    >>> trim_left(["CAG","CG"])
    (1, ['AG', 'G'])

    """
    trimmed = 0
    while all(len(a) > 0 for a in alleles):
        a0 = alleles[0]
        for a in alleles[1:]:
            if a0[0] != a[0]:
                return trimmed, alleles
        alleles = [a[1:] for a in alleles]
        trimmed += 1
    return (trimmed, alleles)


def trim_right(alleles):
    """remove common suffix from right of all alleles, returning
    (number_trimmed, [new_alleles])

    >>> trim_right(["","AA"])
    (0, ['', 'AA'])
    >>> trim_right(["A","AA"])
    (1, ['', 'A'])
    >>> trim_right(["AT","AA"])
    (0, ['AT', 'AA'])
    >>> trim_right(["AA","AA"])
    (2, ['', ''])
    >>> trim_right(["CAG","CG"])
    (1, ['CA', 'C'])

    """
    trimmed = 0
    while all(len(a) > 0 for a in alleles):
        a0 = alleles[0]
        for a in alleles[1:]:
            if a0[-1] != a[-1]:
                return trimmed, alleles
        alleles = [a[:-1] for a in alleles]
        trimmed += 1
    return (trimmed, alleles)


def normalize(ref, pos, allele):
    """
    Normalize allele with respect to reference sequence and position

    >>> ref = "TCTCAGCAGCATCT"
    >>> normalize(ref, (3,3), "CAG")
    ((11, 11), 'GCA')
    >>> normalize(ref, (3,6), "CG")
    ((4, 5), '')
    >>> normalize(ref, (4,4), "AGC")
    ((11, 11), 'GCA')
    >>> normalize(ref, (3,6), "")
    ((8, 11), '')
    >>> normalize(ref, (3,7), "C")
    ((8, 11), '')
    >>> normalize(ref, (3,8), "CA")
    ((8, 11), '')

    >>> ref = "TCTCAGCAGCA"     # normalizes to end of ref
    >>> normalize(ref, (3,3), "CAG")
    ((11, 11), 'GCA')
    """
    start, end = pos
    ref_allele = ref[start:end]

    if debug:                   # pragma: no cover
        _print_seq(ref)
        _print_allele(pos, allele)

    # remove common suffix, prefix, and advance start
    alleles = [ref_allele, allele]
    l_trimmed, alleles = trim_left(alleles)
    start += l_trimmed
    r_trimmed, alleles = trim_right(alleles)
    end -= r_trimmed
    ref_allele, allele = alleles

    while True:
        if debug:                   # pragma: no cover
            _print_allele((start, end), allele)
        if end == len(ref):
            break
        next_residue = ref[end:end + 1]
        trimmed, alleles = trim_left([ref_allele + next_residue, allele + next_residue])
        if trimmed == 0:
            break
        ref_allele, allele = alleles
        start += trimmed
        end += trimmed
    return ((start, end), allele)


def _print_seq(seq):    # pragma: no cover
    p = "{:15s}".format("")
    print(p + " ".join(list("01234567890123456789")))
    print(p + "|" + "|".join(list(seq)) + "|")


def _print_allele(pos, allele):    # pragma: no cover
    start, end = pos
    pfx = "[{:2d},{:2d})".format(start, end)
    p = "{:15s}".format(pfx)
    pos_str = p + "  " * start + "^"
    # pos_str = p + "  "*(start) + "|"
    if end > start:
        pos_str += "—" + "——" * (end - start - 1) + "^"
    print(pos_str + " ⇒ " + allele)


if __name__ == "__main__":    # pragma: no cover
    debug = True
    ref = "TCTCAGCAGCATCT"

    # T|C|T|C|A|G|C|A|G|C|A|T|C|T
    #  |---| ⇒ A

    tests = [
        ((3, 3), "CAG"),
        ((3, 6), "CG"),
        ((4, 4), "AGC"),
        ((3, 6), ""),
        ((3, 7), "C"),
        ((3, 8), "CA"),
    ]
    for t in tests:
        print("* " + str(t))
        print(normalize(ref, *t))
