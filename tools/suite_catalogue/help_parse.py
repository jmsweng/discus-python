"""Parse DISCUS_SUITE ``.hlp`` files into a topic tree.

A topic header is a line ``<level>  <name>`` (level 1-9 in column 1, then spaces). The body
follows; lines starting with ``!p`` right after the header are the syntax lines of a command.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

_HEADER = re.compile(r"^([1-9])\s+(\S.*?)\s*$")


@dataclass
class Topic:
    level: int
    name: str
    line: int
    lines: list = field(default_factory=list)
    children: list = field(default_factory=list)

    @property
    def syntax(self) -> list[str]:
        out = []
        for text in self.lines:
            if text.startswith("!p"):
                out.append(text[2:].strip())
            elif text.strip():
                break
        return [s for s in out if s]

    @property
    def summary(self) -> str:
        body = [t.strip() for t in self.lines if not t.startswith("!p")]
        para = []
        for t in body:
            if not t:
                if para:
                    break
                continue
            para.append(t)
        return " ".join(para)

    def child(self, name: str) -> "Topic | None":
        for c in self.children:
            if c.name == name:
                return c
        return None


def parse_help(path: Path) -> Topic:
    root = Topic(0, "", 0)
    stack = [root]
    for number, raw in enumerate(path.read_text(errors="replace").splitlines(), start=1):
        m = _HEADER.match(raw)
        if m and " " not in m.group(2).strip():
            topic = Topic(int(m.group(1)), m.group(2).strip(), number)
            while stack[-1].level >= topic.level:
                stack.pop()
            stack[-1].children.append(topic)
            stack.append(topic)
        else:
            stack[-1].lines.append(raw.rstrip())
    return root


def match_topic(topics: list[Topic], command: str, min_abbrev: int) -> "Topic | None":
    """Help topic for a command: same name, or one an abbreviation of the other (>= min_abbrev)."""
    for t in topics:
        if t.name == command:
            return t
    for t in topics:
        n = t.name.lower()
        if len(n) >= max(min_abbrev, 2) and (command.startswith(n) or n.startswith(command)):
            return t
    return None
