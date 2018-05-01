from vmc.digest import computed_id

# test use of l.id as computed id when set
class AttrBag:
    pass
o = AttrBag()
o.id = "VMC:bogus"
assert o.id == computed_id(o)


