from typing import Union, Iterable
from itertools import groupby
from operator import attrgetter, itemgetter


# Functions ------------------------------------------------------------------------------------------------------------
def grouper(iterable: Iterable, key: Union[str, int]):
    """Shortcut for sorting and grouping"""
    getter = attrgetter(key) if isinstance(key, str) else itemgetter(key)
    yield from groupby(sorted(iterable, key=getter), key=getter)


def bold(string: str) -> str:
    """Returns the string in bold"""
    return f"\033[1m{string}\033[0m"
