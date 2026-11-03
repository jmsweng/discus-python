"""Static scan of the DISCUS_SUITE Fortran sources for interactive commands.

DISCUS_SUITE dispatches commands with ``str_comp(<command variable>, '<name>', <min>, ...)``
inside IF/ELSEIF chains (``lib_f90/str_comp_mod.f90``). A sub-menu is a routine that appends a
suffix to the prompt (``prompt = prompt(...) // '/fourier'``) and then reads commands until
``exit``. This module finds, per routine:

* the command branches (name, minimum abbreviation, source line, line span),
* the keyword comparisons inside each branch (``str_comp(cpara(n), 'word', ...)``),
* the optional ``name:value`` parameter names (``DATA oname / ... /`` with ``get_optional``),
* the routines each branch calls,
* the prompt suffix, if the routine is a menu.

It is a heuristic reader of regular code, not a Fortran parser; ``build_catalogue.py`` records
what it could not resolve instead of guessing.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

# Variables that hold the command word in the dispatchers of the suite's programs.
COMMAND_VARIABLES = {"befehl", "bef", "cbef", "cbefehl", "befehl_s"}

_STR_COMP = re.compile(
    r"str_comp\s*\(\s*([A-Za-z_][\w%]*(?:\s*\(\s*[^()]*\))?)\s*,\s*'([^']*)'\s*,\s*(\d+)\s*,"
    r"\s*[^,()]+\s*,\s*(\d+)?", re.I
)
# Length caps on the command word: CHARACTER(LEN=n) :: befehl, and lbef = min(indxb - 1, n)
_CMD_DECL = re.compile(r"^\s*character\s*\(\s*(?:len\s*=\s*)?(\d+)\s*\)(?:\s*,[^:]*::|\s*::)?\s*(.*)$", re.I)
_LBEF_CAP = re.compile(r"\blbef\s*=\s*min\s*\(\s*indxb\s*-\s*1\s*,\s*(\d+)\s*\)", re.I)
_ROUTINE_START = re.compile(
    r"^\s*(?:(?:recursive|pure|elemental|module|impure)\s+)*"
    r"(?:(?:integer|logical|real|complex|character|double\s+precision|type)(?:\s*\([^)]*\))?\s+)?"
    r"(subroutine|function)\s+(\w+)", re.I
)
_ROUTINE_END = re.compile(r"^\s*end\s*(subroutine|function)\b", re.I)
_PROMPT_SUFFIX = re.compile(r"\bprompt\s*=\s*.*?//\s*'(/[^']+)'", re.I)
_CALL = re.compile(r"\bcall\s+(\w+)", re.I)
_DATA_ONAME = re.compile(r"\bdata\s+oname\s*/(.*?)/", re.I)
_ONAME_ASSIGN = re.compile(r"\boname\s*=\s*\(/(.*?)/\)", re.I)
_ONAME_ASSIGN2 = re.compile(r"\boname\s*=\s*\[(.*?)\]", re.I)
_IF_THEN = re.compile(r"^\s*(\w+\s*:\s*)?if\b.*\bthen\b", re.I)
_END_IF = re.compile(r"^\s*end\s*if\b", re.I)
_RETURN = re.compile(r"(^\s*|\)\s*)return\s*$", re.I)
_GET_OPTIONAL = re.compile(r"\bget_optional\s*\(", re.I)
_DO_HEL = re.compile(r"\bdo_hel\s*\(\s*'([^']+)'", re.I)
_QUOTED = re.compile(r"'([^']*)'")


@dataclass
class Logical:
    """One logical Fortran line (continuations joined), with its first physical line number."""

    line: int
    text: str


@dataclass
class Branch:
    name: str
    min_abbrev: int
    line: int
    end_line: int = 0
    keywords: dict = field(default_factory=dict)  # word -> min_abbrev
    calls: set = field(default_factory=set)
    optional: list = field(default_factory=list)
    uses_get_optional: bool = False
    compare_len: int | None = None  # 5th str_comp argument (caps the compared length)
    body_depth: int = 0  # IF-nesting depth of the branch body
    nested_calls: set = field(default_factory=set)  # calls inside a further IF within the branch


@dataclass
class Routine:
    name: str
    file: str
    line: int
    end_line: int = 0
    prompts: list = field(default_factory=list)
    branches: list = field(default_factory=list)
    calls: set = field(default_factory=set)
    optional: list = field(default_factory=list)
    uses_get_optional: bool = False
    keywords: dict = field(default_factory=dict)  # keyword tests outside any command branch
    help_paths: list = field(default_factory=list)  # literal prefixes passed to do_hel
    word_length: int | None = None  # longest command word the dispatcher sees
    assign_exclusions: list = field(default_factory=list)  # (name, min, len) not treated as "x = y"
    fallback_calls: dict = field(default_factory=dict)  # callee -> first line it is called
    early_return: bool = False  # a RETURN before the prompt is changed (menu may not open)
    if_depth: int = 0


def strip_comment(text: str) -> str:
    """Remove a ``!`` comment that is not inside a character constant."""
    quote = None
    for i, ch in enumerate(text):
        if quote:
            if ch == quote:
                quote = None
        elif ch in ("'", '"'):
            quote = ch
        elif ch == "!":
            return text[:i]
    return text


def logical_lines(path: Path) -> list[Logical]:
    """Join ``&`` continuations; drop comments and blank lines."""
    out: list[Logical] = []
    pending: Logical | None = None
    for number, raw in enumerate(path.read_text(errors="replace").splitlines(), start=1):
        text = strip_comment(raw).rstrip()
        if not text.strip():
            continue
        stripped = text.lstrip()
        if pending is not None:
            if stripped.startswith("&"):
                pending.text += stripped[1:]  # a leading '&' resumes without a blank
            else:
                pending.text += " " + stripped
        else:
            pending = Logical(number, text)
        if pending.text.rstrip().endswith("&"):
            pending.text = pending.text.rstrip()[:-1]
            continue
        out.append(pending)
        pending = None
    if pending is not None:
        out.append(pending)
    return out


def _words(text: str) -> list[str]:
    return [w.strip() for w in _QUOTED.findall(text) if w.strip()]


def scan_file(path: Path, root: Path) -> list[Routine]:
    routines: list[Routine] = []
    stack: list[Routine] = []
    rel = path.relative_to(root).as_posix()
    for logical in logical_lines(path):
        text = logical.text
        start = _ROUTINE_START.match(text)
        if start and not re.match(r"^\s*end\b", text, re.I):
            routine = Routine(start.group(2).lower(), rel, logical.line)
            routines.append(routine)
            stack.append(routine)
            continue
        if _ROUTINE_END.match(text):
            if stack:
                closed = stack.pop()
                closed.end_line = logical.line
                if closed.branches:
                    closed.branches[-1].end_line = closed.branches[-1].end_line or logical.line
            continue
        if not stack:
            continue
        routine = stack[-1]
        if _END_IF.match(text):
            routine.if_depth -= 1
        opens_if = bool(_IF_THEN.match(text))
        if _RETURN.search(text) and not routine.prompts:
            routine.early_return = True
        for m in _CMD_DECL.finditer(text):
            declared = {v.strip().split("=")[0].strip().lower() for v in m.group(2).split(",")}
            if declared & COMMAND_VARIABLES:
                routine.word_length = min(routine.word_length or 10**6, int(m.group(1)))
        for m in _LBEF_CAP.finditer(text):
            routine.word_length = min(routine.word_length or 10**6, int(m.group(1)))
        for m in _DO_HEL.finditer(text):
            path = m.group(1).strip()
            if path and path not in routine.help_paths:
                routine.help_paths.append(path)
        for m in _PROMPT_SUFFIX.finditer(text):
            if m.group(1) not in routine.prompts:
                routine.prompts.append(m.group(1))
        names = []
        for pattern in (_DATA_ONAME, _ONAME_ASSIGN, _ONAME_ASSIGN2):
            for m in pattern.finditer(text):
                names.extend(_words(m.group(1)))
        current = routine.branches[-1] if routine.branches else None
        if names:
            (current.optional if current else routine.optional).extend(names)
        if _GET_OPTIONAL.search(text):
            if current:
                current.uses_get_optional = True
            routine.uses_get_optional = True
        new_branch = False
        matches = list(_STR_COMP.finditer(text))
        guard_line = any(
            re.sub(r"\s+", "", m.group(1)).lower() in COMMAND_VARIABLES
            and re.search(r"\.not\.\s*\(?\s*$", text[: m.start()], re.I)
            for m in matches
        )
        if guard_line and "indxg" in text.lower():
            for m in matches:
                if re.sub(r"\s+", "", m.group(1)).lower() in COMMAND_VARIABLES:
                    entry = (m.group(2).strip(), int(m.group(3)), int(m.group(4)) if m.group(4) else None)
                    if entry not in routine.assign_exclusions:
                        routine.assign_exclusions.append(entry)
        for m in ([] if guard_line else matches):
            var = re.sub(r"\s+", "", m.group(1)).lower()
            word, minimum = m.group(2).strip(), int(m.group(3))
            if not word:
                continue
            before = text[: m.start()]
            negated = re.search(r"\.not\.\s*\(?\s*$", before, re.I)
            if var in COMMAND_VARIABLES and negated:
                continue
            if var in COMMAND_VARIABLES and re.search(r"\b(else\s*)?if\b", before, re.I):
                if any(b.name == word for b in routine.branches):
                    continue  # nested re-test of the same command inside its branch (ll/lr/ul/tl)
                if not new_branch:
                    if current:
                        current.end_line = logical.line - 1
                    current = Branch(word, minimum, logical.line,
                                     compare_len=int(m.group(4)) if m.group(4) else None,
                                     body_depth=routine.if_depth + (1 if opens_if else 0))
                    routine.branches.append(current)
                    new_branch = True
                else:
                    # several names tested in one condition (ll/lr/ul/tl, help/?) -> separate entries
                    extra = Branch(word, minimum, logical.line, body_depth=current.body_depth,
                                   compare_len=int(m.group(4)) if m.group(4) else None)
                    extra.calls = current.calls
                    extra.nested_calls = current.nested_calls
                    routine.branches.append(extra)
            elif var not in COMMAND_VARIABLES:
                (current.keywords if current is not None else routine.keywords).setdefault(word, minimum)
        branch = routine.branches[-1] if routine.branches else None
        depth_here = routine.if_depth + (1 if opens_if else 0)
        for m in _CALL.finditer(text):
            callee = m.group(1).lower()
            routine.calls.add(callee)
            routine.fallback_calls.setdefault(callee, logical.line)
            if branch is not None:
                branch.calls.add(callee)
                inline_if = bool(re.match(r"^\s*if\s*\(", text, re.I)) and not opens_if
                if depth_here > branch.body_depth or inline_if and branch.line != logical.line:
                    branch.nested_calls.add(callee)
        if opens_if:
            routine.if_depth += 1
    return routines


def scan_tree(root: Path, directories: list[str]) -> dict[str, list[Routine]]:
    """Scan every ``*.f90`` file under the given directories; returns routines by name."""
    by_name: dict[str, list[Routine]] = {}
    for directory in directories:
        for path in sorted((root / directory).rglob("*.f90")):
            for routine in scan_file(path, root):
                by_name.setdefault(routine.name, []).append(routine)
    return by_name
