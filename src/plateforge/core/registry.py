"""Named handler registries. This is the 'leave room for later' mechanism.

Adding a reader, a construct format, a hit-caller, or a reagent rule means
writing a new file and registering it. It never means editing existing code.
"""
from __future__ import annotations

from typing import Callable, Generic, Iterator, TypeVar

T = TypeVar("T")


class Registry(Generic[T]):
    def __init__(self, kind: str) -> None:
        self.kind = kind
        self._items: dict[str, T] = {}
        self._meta: dict[str, dict] = {}

    def register(self, key: str, **meta) -> Callable[[T], T]:
        """Use as a decorator: @READERS.register('spectramax', ext='.txt')"""
        def deco(obj: T) -> T:
            if key in self._items:
                raise KeyError(f"{self.kind} {key!r} already registered")
            self._items[key] = obj
            self._meta[key] = meta
            return obj
        return deco

    def add(self, key: str, obj: T, **meta) -> None:
        self.register(key, **meta)(obj)

    def get(self, key: str) -> T:
        try:
            return self._items[key]
        except KeyError:
            raise KeyError(f"unknown {self.kind} {key!r}; known: {sorted(self._items)}") from None

    def meta(self, key: str) -> dict:
        self.get(key)
        return dict(self._meta[key])

    def keys(self) -> list[str]:
        return sorted(self._items)

    def __contains__(self, key: object) -> bool:
        return key in self._items

    def __iter__(self) -> Iterator[str]:
        return iter(sorted(self._items))

    def __len__(self) -> int:
        return len(self._items)

    def __repr__(self) -> str:
        return f"<Registry {self.kind}: {', '.join(self.keys()) or 'empty'}>"


READERS = Registry("plate reader adapter")
CONSTRUCTS = Registry("construct format")
CALLERS = Registry("hit-calling strategy")
CONTROLS = Registry("control resolver")
EMITTERS = Registry("output file emitter")
RULES = Registry("reagent caveat rule")
