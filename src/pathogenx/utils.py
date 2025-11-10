from typing import Union, Iterable
from itertools import groupby
from operator import attrgetter, itemgetter


def grouper(iterable: Iterable, key: Union[str, int]):
    """Shortcut for sorting and grouping"""
    getter = attrgetter(key) if isinstance(key, str) else itemgetter(key)
    yield from groupby(sorted(iterable, key=getter), key=getter)