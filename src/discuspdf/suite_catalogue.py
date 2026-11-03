"""Catalogue of DISCUS_SUITE commands, and a checker for macros written against it.

The catalogue is generated from the DiffuseCode Fortran sources by
``tools/suite_catalogue/build_catalogue.py`` and shipped as package data. It lists, for every
menu of the suite (``suite``, ``discus``, ``discus/fourier``, ``kuplot/fit``, ...), the commands
the source dispatches, their minimum abbreviation, the sub-menu a command opens and the help
text. It describes what the engine *accepts*; it says nothing about what a command does
physically.

>>> cat = load_catalogue()
>>> cat.resolve("discus", "four").name
'fourier'
>>> [step.command for step in cat.walk("discus\\nfourier\\nrun\\nexit\\n") if step.command]
['discus', 'fourier', 'run', 'exit']
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from importlib import resources

CATALOGUE_FILES = {"6.23.01": "discus_suite_6.23.01.json"}
SCHEMA = "discuspdf.suite_catalogue"

# Lines handled by the macro machinery of lib_f90 before any menu sees them.
CONTROL_WORDS = {"do", "enddo", "if", "elseif", "else", "endif", "break"}
# Menus that process exactly one command and return to their parent (verified in the source).
ONE_SHOT_MENUS = {"discus/read"}


def str_comp(typed: str, name: str, minimum: int, compare_len: int | None = None) -> bool:
    """DISCUS's ``str_comp(a, b, j, la, lb)`` (``lib_f90/str_comp_mod.f90``).

    It compares the first ``i = min(la, lb)`` characters, and only if ``i >= j``. ``la`` is
    the length of the typed word, ``lb`` the length argument of the call (``compare_len``,
    usually ``len(name)`` but sometimes shorter, which truncates the comparison). The name is
    blank-padded to ``lb``. Case-sensitive; a word longer than the name still matches
    (``fourierx`` matches ``fourier``)."""
    lb = len(name) if compare_len is None else compare_len
    if not typed or lb <= 0:
        return False
    i = min(len(typed), lb)
    return i >= minimum and typed[:i] == name.ljust(lb)[:i]


@dataclass(frozen=True)
class Command:
    menu: str
    name: str
    min_abbrev: int
    source: str
    compare_len: int | None = None
    keywords: dict = field(default_factory=dict, hash=False, compare=False)
    optional: tuple = ()
    enters: str | None = None
    standard: bool = False
    enters_when: dict | None = field(default=None, hash=False, compare=False)
    help: dict | None = field(default=None, hash=False, compare=False)

    @property
    def syntax(self) -> list[str]:
        return list((self.help or {}).get("syntax", []))


@dataclass(frozen=True)
class Menu:
    path: str
    program: str
    routine: str
    commands: tuple
    inherits: tuple = ()
    help_only: tuple = ()
    inherits_first: tuple = ()
    word_length: int | None = None

    def command(self, name: str) -> Command | None:
        return next((c for c in self.commands if c.name == name), None)


@dataclass(frozen=True)
class Step:
    """One macro line as the engine would see it."""

    line: int
    text: str
    menu: str  # menu the line is read in
    kind: str  # "command", "comment", "control", "assignment", "macro", "unknown"
    command: str | None = None  # resolved full command name
    found_in: str | None = None  # menu that supplied the command (may be an inherited one)


class Catalogue:
    def __init__(self, data: dict):
        if data.get("schema") != SCHEMA:
            raise ValueError(f"not a {SCHEMA} file")
        self.data = data
        self.engine = data["engine"]
        self.unresolved = data.get("unresolved", [])
        self.menus: dict[str, Menu] = {}
        for path, m in data["menus"].items():
            cmds = tuple(
                Command(path, c["name"], c["min_abbrev"], c["source"], c.get("compare_len"),
                        c.get("keywords", {}),
                        tuple(c.get("optional", ())), c.get("enters"), bool(c.get("standard")),
                        c.get("enters_when"), c.get("help"))
                for c in m["commands"]
            )
            self.menus[path] = Menu(path, m["program"], m["routine"], cmds, tuple(m.get("inherits", ())),
                                    tuple(m.get("help_only", ())), tuple(m.get("inherits_first", ())),
                                    m.get("word_length"))

    def __repr__(self):
        return (f"Catalogue(DISCUS_SUITE {self.engine['version']} @ {self.engine['commit'][:7]}: "
                f"{len(self.menus)} menus, {self.n_commands} commands)")

    @property
    def n_commands(self) -> int:
        return sum(len(m.commands) for m in self.menus.values())

    def menu(self, path: str) -> Menu:
        try:
            return self.menus[path]
        except KeyError:
            raise KeyError(f"no menu {path!r}; known: {', '.join(self.menus)}") from None

    def commands(self, include_standard: bool = False):
        for m in self.menus.values():
            for c in m.commands:
                if include_standard or not c.standard:
                    yield c

    def resolve(self, menu: str, word: str) -> Command | None:
        """The command the engine runs for ``word`` typed in ``menu``.

        The word is first cut to the length the menu's dispatcher keeps (KUPLOT and MIXSCAT
        keep 4 characters). Menus listed in ``inherits_first`` are asked before the menu's own
        branches; then the branches in source order (the first match wins); then the menus it
        falls through to."""
        start = self.menus.get(menu)
        if start is None:
            return None
        if start.word_length:
            word = word[: start.word_length]
        seen: set = set()

        def search(path):
            if path in seen or path not in self.menus:
                return None
            seen.add(path)
            m = self.menus[path]
            for first in m.inherits_first:
                found = search(first)
                if found is not None:
                    return found
            for c in m.commands:
                if str_comp(word, c.name, c.min_abbrev, c.compare_len):
                    return c
            for later in m.inherits:
                if later in m.inherits_first:
                    continue
                found = search(later)
                if found is not None:
                    return found
            return None

        return search(menu)

    def shortest_abbreviation(self, cmd: Command) -> str | None:
        """Shortest prefix of ``cmd.name`` that the engine resolves to ``cmd`` in its own menu
        (``None`` if no prefix does, e.g. a branch shadowed by an earlier one)."""
        for k in range(1, len(cmd.name) + 1):
            if self.resolve(cmd.menu, cmd.name[:k]) == cmd:
                return cmd.name[:k]
        return None

    def walk(self, text: str, start: str = "suite"):
        """Follow a macro line by line, tracking the menu the engine is in.

        Yields one :class:`Step` per non-blank line. Unknown commands are reported as
        ``kind="unknown"`` and do not change the menu. Control flow (``do``/``if``) is not
        executed; every line is read once, in order."""
        stack = [start]
        for number, raw in enumerate(text.splitlines(), start=1):
            line = raw.strip()
            if not line:
                continue
            menu = stack[-1]
            if line[0] in "#!":
                yield Step(number, line, menu, "comment")
                continue
            if line.startswith("@"):
                yield Step(number, line, menu, "macro")
                continue
            word = line.split()[0].split(",")[0]
            lowered = word.lower()
            if lowered in CONTROL_WORDS:
                yield Step(number, line, menu, "control", lowered)
                continue
            if "=" in line and lowered not in ("echo", "system", "help", "?"):
                yield Step(number, line, menu, "assignment")
                continue
            cmd = self.resolve(menu, word)
            if cmd is None:
                yield Step(number, line, menu, "unknown")
                continue
            yield Step(number, line, menu, "command", cmd.name, cmd.menu)
            if cmd.name == "exit":
                if len(stack) > 1:
                    stack.pop()
                continue
            if menu in ONE_SHOT_MENUS:
                stack.pop()
            if cmd.enters and cmd.enters in self.menus and self._opens(cmd, line):
                stack.append(cmd.enters)

    @staticmethod
    def _opens(cmd: Command, line: str) -> bool:
        """Whether this command line opens ``cmd.enters`` (see ``enters_when`` rules)."""
        rule = cmd.enters_when
        if not rule:
            return True
        args = line.split(None, 1)[1].strip() if len(line.split(None, 1)) > 1 else ""
        params = [p.strip() for p in args.split(",")] if args else []
        if rule["kind"] == "no_arguments":
            return not params
        if rule["kind"] == "first_word":
            first = params[0] if params else ""
            return any(str_comp(first, w, rule["min_abbrev"]) for w in rule["words"])
        if rule["kind"] == "option":
            return any(p.replace(" ", "") == f"{rule['name']}:{rule['value']}" for p in params)
        raise ValueError(f"unknown entry rule {rule['kind']!r}")

    def check(self, text: str, start: str = "suite") -> list[Step]:
        """Lines a macro uses that the catalogue does not know (empty list = all known)."""
        return [s for s in self.walk(text, start) if s.kind == "unknown"]


@lru_cache(maxsize=None)
def load_catalogue(version: str = "6.23.01") -> Catalogue:
    """The bundled catalogue for a DISCUS_SUITE version (currently 6.23.01 only; the owner's
    fork at ``a85c9de`` has the same command set, see ``docs/SUITE_CATALOGUE.md``)."""
    try:
        name = CATALOGUE_FILES[version]
    except KeyError:
        raise ValueError(f"no catalogue for DISCUS_SUITE {version}; have {sorted(CATALOGUE_FILES)}") from None
    text = resources.files("discuspdf").joinpath("data", "suite_catalogue", name).read_text(encoding="utf-8")
    return Catalogue(json.loads(text))
