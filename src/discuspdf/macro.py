"""Build DISCUS_SUITE macros line by line, checked against the command catalogue.

Every command of the suite is reachable this way, also those without a typed discuspdf API
(``docs/SUITE_COVERAGE.md`` lists which have one). The builder writes plain DISCUS macro text;
:meth:`Macro.render` walks it through :mod:`discuspdf.suite_catalogue` and refuses a macro
with an unknown command or a line that would run in a different menu than the code says.

>>> m = Macro("hk0 plane of GaAs")
>>> with m.program("discus"):
...     with m.menu("read"):
...         _ = m.command("cell", "gaas.cll", 4, 4, 4)
...     with m.menu("fourier"):
...         _ = m.command("xray").command("ll", -2, -2, 0).command("run")
>>> print(m.render())  # doctest: +NORMALIZE_WHITESPACE
# hk0 plane of GaAs
discus
read
cell gaas.cll,4,4,4
fourier
xray
ll -2,-2,0
run
exit
exit
<BLANKLINE>

The builder checks command names and menus, not argument meaning: what a command does
physically is only known from an engine run (see the verification levels in
``docs/SUITE_CATALOGUE.md``).
"""

from __future__ import annotations

import math
import numbers
from contextlib import contextmanager
from dataclasses import dataclass

from .suite_catalogue import ONE_SHOT_MENUS, Catalogue, load_catalogue


class MacroError(ValueError):
    """A macro the engine would not run as written."""


def format_value(value) -> str:
    """A DISCUS parameter: integers as written, floats with up to 12 significant digits,
    strings unchanged (file names, keywords, expressions such as ``n[1]``)."""
    if isinstance(value, bool):
        raise MacroError("booleans are not DISCUS values; write the keyword (e.g. 'on'/'off')")
    if isinstance(value, numbers.Integral):
        return str(int(value))
    if isinstance(value, numbers.Real):
        x = float(value)
        if not math.isfinite(x):
            raise MacroError(f"non-finite number {value!r}")
        return f"{x:.12g}"
    if isinstance(value, str):
        text = value.strip()
        if not text:
            raise MacroError("empty parameter")
        _check_text(text)
        return text
    raise MacroError(f"cannot write {type(value).__name__} {value!r} as a DISCUS parameter")


def _check_text(text: str) -> None:
    if "\n" in text or "\r" in text:
        raise MacroError(f"line break in {text!r}")
    if text.rstrip().endswith("&"):
        raise MacroError(f"{text!r} ends with '&', which joins it to the next line")
    quote = None
    for ch in text:
        if quote:
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
        elif ch in "#!":
            raise MacroError(f"{text!r}: '{ch}' outside quotes starts a comment in DISCUS")


@dataclass(frozen=True)
class _Line:
    text: str
    menu: str | None  # menu the builder expects the line to be read in (None: not checked)
    kind: str  # command, raw, comment, exit


class Macro:
    """A DISCUS_SUITE macro, starting at the suite level (as ``discus_suite -macro`` does)."""

    def __init__(self, title: str | None = None, *, catalogue: Catalogue | None = None, start: str = "suite"):
        self.catalogue = catalogue or load_catalogue()
        if start not in self.catalogue.menus:
            raise MacroError(f"unknown start menu {start!r}")
        self._lines: list[_Line] = []
        self._stack = [start]
        if title:
            self.comment(title)

    # -- building --------------------------------------------------------------------------

    @property
    def current_menu(self) -> str:
        return self._stack[-1]

    def comment(self, text: str) -> "Macro":
        for part in str(text).splitlines() or [""]:
            self._lines.append(_Line(f"# {part}".rstrip(), None, "comment"))
        return self

    def command(self, word: str, *params, **options) -> "Macro":
        """One command line: ``word p1,p2,...,name:value``. ``options`` become ``name:value``
        parameters after the positional ones (a trailing ``_`` in the name is dropped, so
        ``global_=...`` writes ``global:``)."""
        cmd = self._resolve(word)
        line = self._format(word, params, options)
        self._lines.append(_Line(line, self.current_menu, "command"))
        if self.current_menu in ONE_SHOT_MENUS:
            self._stack.pop()
        if cmd.name == "exit":
            if len(self._stack) == 1:
                raise MacroError("'exit' at the top level ends the run; use the context managers")
            self._stack.pop()
        elif cmd.enters and self.catalogue._opens(cmd, line):
            raise MacroError(f"{word!r} opens the menu {cmd.enters!r}; use 'with macro.menu(...)'")
        return self

    def raw(self, line: str) -> "Macro":
        """A line written exactly as given (the escape hatch for anything the other methods do
        not cover: assignments, ``do``/``if`` blocks, ``@macro`` calls). It is still read by the
        catalogue when the macro is rendered; menu changes it causes are not tracked."""
        text = str(line).strip()
        if "\n" in text:
            raise MacroError("raw() takes one line")
        self._lines.append(_Line(text, None, "raw"))
        return self

    @contextmanager
    def menu(self, word: str, *params, **options):
        """Enter the sub-menu ``word`` opens, and leave it with ``exit`` at the end of the block
        (``read`` leaves by itself after one line)."""
        cmd = self._resolve(word)
        line = self._format(word, params, options)
        if not cmd.enters or not self.catalogue._opens(cmd, line):
            raise MacroError(f"{line!r} does not open a menu in {self.current_menu!r}")
        self._lines.append(_Line(line, self.current_menu, "command"))
        self._stack.append(cmd.enters)
        depth = len(self._stack)
        yield self
        if cmd.enters in ONE_SHOT_MENUS:
            if len(self._stack) >= depth:
                raise MacroError(f"the {cmd.enters!r} menu takes exactly one command")
            return
        if len(self._stack) != depth or self._stack[-1] != cmd.enters:
            raise MacroError(f"menu nesting broken inside {word!r}")
        self._lines.append(_Line("exit", cmd.enters, "exit"))
        self._stack.pop()

    def program(self, name: str):
        """Enter a program of the suite (``discus``, ``kuplot``, ``refine``, ``diffev``)."""
        return self.menu(name)

    # -- output ----------------------------------------------------------------------------

    def problems(self) -> list[str]:
        """Lines the engine would read differently from what the builder intended."""
        out = []
        if len(self._stack) != 1:
            out.append(f"macro ends inside {self.current_menu!r}")
        text = "\n".join(line.text for line in self._lines) + "\n"
        steps = list(self.catalogue.walk(text, self._stack[0]))
        by_line = {s.line: s for s in steps}
        for number, line in enumerate(self._lines, start=1):
            step = by_line.get(number)
            if step is None:
                continue
            if step.kind in ("unknown", "after_end"):
                out.append(f"line {number}: {line.text!r} is {step.kind} in {step.menu!r}")
            elif line.menu is not None and step.menu != line.menu:
                out.append(f"line {number}: {line.text!r} would run in {step.menu!r}, not {line.menu!r}")
        return out

    def render(self) -> str:
        problems = self.problems()
        if problems:
            raise MacroError("; ".join(problems))
        return "\n".join(line.text for line in self._lines) + "\n"

    def __str__(self) -> str:
        return self.render()

    # -- helpers ---------------------------------------------------------------------------

    def _resolve(self, word: str):
        if not isinstance(word, str) or not word or any(c.isspace() for c in word) or "," in word:
            raise MacroError(f"invalid command word {word!r}")
        cmd = self.catalogue.resolve(self.current_menu, word)
        if cmd is None:
            raise MacroError(f"{word!r} is not a command of {self.current_menu!r} "
                             f"(DISCUS_SUITE {self.catalogue.engine['version']})")
        return cmd

    @staticmethod
    def _format(word: str, params, options) -> str:
        parts = [format_value(p) for p in params]
        for name, value in options.items():
            key = name[:-1] if name.endswith("_") else name
            parts.append(f"{key}:{format_value(value)}")
        line = word if not parts else f"{word} {','.join(parts)}"
        _check_text(line)
        return line
