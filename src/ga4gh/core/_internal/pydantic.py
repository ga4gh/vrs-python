def is_identifiable(o: any) -> bool:
    """Determine if object is identifiable. An object is considered identifiable if
    contains a `ga4gh_digest_keys` attribute

    :param o: Object
    :return: `True` if `o` has `ga4gh_digest_keys` attribute. `False` otherwise.
    """
    return hasattr(o, "ga4gh_digest_keys")


def is_literal(o: any) -> bool:
    return isinstance(o, (str, int, float, complex, bool))


def is_list(o: any) -> bool:
    return isinstance(o, list)


def is_curie_type(o: any) -> bool:
    # return isinstance(o, CURIE)
    pass

def is_pydantic_instance(o: any) -> bool:
    # return isinstance(o, BaseModel)
    pass
