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
# Menus that process exactly one line and return to their parent (structur.f90:135-145).
ONE_SHOT_MENUS = {"discus/read"}
# Where do/if blocks are recognised (discus_loop.f90, kuplot_loop.f90, ... test line(1:3)=='do ').
PROGRAM_LEVELS = {"suite", "discus", "kuplot", "diffev", "refine", "mixscat"}


def _logical_lines(text: str):
    """(first line number, text) with ``&`` continuations joined (lib_f90/macro.f90)."""
    pending, start = "", 0
    for number, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not pending and not line:
            continue
        if not pending:
            start = number
        if line.endswith("&"):
            pending += line[:-1]
            continue
        pending += line
        if pending.strip():
            yield start, pending.strip()
        pending = ""
    if pending.strip():
        yield start, pending.strip()


def _remove_comment(line: str) -> str:
    """``remove_comment`` (lib_f90/sup.f90): cut at ``#`` or ``!`` outside quotes."""
    search = True
    for i, ch in enumerate(line):
        if ch in "'\"":
            search = not search
        elif search and ch in "#!":
            return line[:i].rstrip()
    return line


def str_comp(typed: str, name: str, minimum: int, compare_len: int | None = None) -> bool:
    """DISCUS's ``str_comp(a, b, j, la, lb)`` (``lib_f90/str_comp_mod.f90``) for a command
    word, which never contains a blank (for blank-containing strings the Fortran also cuts at
    the first blank; that case is not ported).

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
    assignments: bool = True
    assign_exclusions: tuple = ()

    def command(self, name: str) -> Command | None:
        return next((c for c in self.commands if c.name == name), None)


@dataclass(frozen=True)
class Step:
    """One macro line as the engine would see it."""

    line: int
    text: str
    menu: str  # menu the line is read in
    kind: str  # "command", "comment", "control", "assignment", "macro", "unknown", "after_end"
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
                                    m.get("word_length"), bool(m.get("assignments", True)),
                                    tuple(tuple(e) for e in m.get("assign_exclusions", ())))

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

    def _is_assignment(self, menu: Menu, line: str, word: str) -> bool:
        """``x = y`` lines go to ``do_math`` unless the command is one the menu exempts
        (``echo``, ``system``, ``help``, KUPLOT's ``tit1``...) or the line contains ``==``."""
        if not menu.assignments or "=" not in line or "==" in line:
            return False
        if menu.word_length:
            word = word[: menu.word_length]
        return not any(str_comp(word, n, m, lb) for n, m, lb in menu.assign_exclusions)

    def walk(self, text: str, start: str = "suite"):
        """Follow a macro line by line, tracking the menu the engine is in.

        Yields one :class:`Step` per logical line, as ``lib_f90`` reads it: lines ending in
        ``&`` are joined, leading blanks dropped, ``#``/``!`` comments removed outside quotes,
        and the command word ends at the first blank or tab. ``do``/``if`` blocks are
        recognised only at a program's top level, as in the engine. Blocks and ``@macro``
        calls are not executed: every line is read once, in order. Unknown commands are
        ``kind="unknown"`` and do not change the menu; lines after the final ``exit`` at the
        suite level are ``kind="after_end"``."""
        stack = [start]
        depth = 0  # open do/if blocks
        ended = False
        for number, line in _logical_lines(text):
            menu = stack[-1]
            if ended:
                yield Step(number, line, menu, "after_end")
                continue
            if line[0] in "#!":
                yield Step(number, line, menu, "comment")
                continue
            line = _remove_comment(line)
            if not line:
                yield Step(number, line, menu, "comment")
                continue
            m = self.menus[menu]
            if menu in PROGRAM_LEVELS or depth:
                if line[:3] == "do " or line[:2] == "if":
                    depth += 1
                    yield Step(number, line, menu, "control", line.split("(")[0].split()[0])
                    continue
                if depth:
                    compact = line.replace(" ", "")
                    if compact.startswith(("enddo", "endif")):
                        depth -= 1
                        yield Step(number, line, menu, "control", compact[:5])
                        continue
                    if compact.startswith(("else", "break")):
                        yield Step(number, line, menu, "control", "elseif" if compact.startswith("elseif")
                                   else compact[:4] if compact.startswith("else") else "break")
                        continue
            word = line.split()[0]
            one_shot = menu in ONE_SHOT_MENUS
            if self._is_assignment(m, line, word):
                yield Step(number, line, menu, "assignment")
            elif line.startswith("@"):
                yield Step(number, line, menu, "macro")
            else:
                cmd = self.resolve(menu, word)
                if cmd is None:
                    yield Step(number, line, menu, "unknown")
                else:
                    yield Step(number, line, menu, "command", cmd.name, cmd.menu)
                    if cmd.name == "exit":
                        if len(stack) > 1:
                            stack.pop()
                        else:
                            ended = True
                        continue
                    if one_shot:
                        stack.pop()
                        one_shot = False
                    if cmd.enters and cmd.enters in self.menus and self._opens(cmd, line):
                        stack.append(cmd.enters)
            if one_shot:  # read leaves after any non-comment line, whatever happened to it
                stack.pop()

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
