"""decorators for vrs-python"""


from collections.abc import Callable


def lazy_property(fn: Callable):  # noqa: ANN201
    """Provide a decorator that makes a property lazy-evaluated.

    [mv]
    """
    attr_name = "_lazy_" + fn.__name__

    @property
    def _lazy_property(self):  # noqa: ANN001 ANN202
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fn(self))
        return getattr(self, attr_name)

    return _lazy_property
