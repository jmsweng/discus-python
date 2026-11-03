"""Build the DISCUS_SUITE command catalogue from a DiffuseCode source tree.

    python tools/suite_catalogue/build_catalogue.py <DiffuseCode checkout> <out.json> \
        --version 6.23.01 --commit <sha> [--label stock]

Every menu (program top level or sub-menu) lists its commands as found in the Fortran
dispatchers, merged with the help text. The source decides what exists; help only adds text.
Anything the scanner could not resolve is listed under ``unresolved`` rather than guessed.
"""

from __future__ import annotations

import argparse
import json
from collections import deque
from pathlib import Path

from fortran_scan import Routine, scan_tree
from help_parse import Topic, match_topic, parse_help

SCHEMA = "discuspdf.suite_catalogue"
SCHEMA_VERSION = 1

# program -> (source directory, root dispatcher routine, help file)
PROGRAMS = {
    "suite": ("suite/prog", "suite_mache_kdo", "suite/prog/appl_suite.hlp"),
    "discus": ("discus/prog", "discus_mache_kdo", "discus/prog/appl_dis.hlp"),
    "kuplot": ("kuplot/prog", "kuplot_mache_kdo", "kuplot/prog/appl_kup.hlp"),
    "diffev": ("diffev/prog", "diffev_mache_kdo", "diffev/prog/appl_dif.hlp"),
    "refine": ("refine/prog", "refine_mache_kdo", "refine/prog/appl_ref.hlp"),
    "mixscat": ("mixscat/prog", "mixscat_mache_kdo", "mixscat/prog/appl_mix.hlp"),
}
LIBRARY = ("lib_f90", "kdo_all", "lib_f90/lib_f90.hlp")
# Routines that hold commands shared by several menus (reached when no branch matches).
FALLBACKS = {"kdo_all": "common", "discus_kdo_common": "discus/common", "kuplot_kdo_common": "kuplot/common"}
# Library helpers never followed when collecting keywords/options of a command.
OPAQUE = {"get_params", "ber_params", "get_optional", "do_cap", "do_low", "del_params", "errlist",
          "no_error", "do_hel", "do_math", "file_kdo", "macro_close", "do_build_name", "get_cmd",
          "do_eval", "echo", "do_operating", "do_input", "macro_continue", "kdo_all",
          "discus_kdo_common", "kuplot_kdo_common", "get_iscat"}
CALL_DEPTH = 3
# When a command opens its sub-menu only in some cases. Verified by reading the source
# (cited); the builder flags other candidates as "enters_conditional" for review.
ENTRY_RULES = {
    ("discus", "save"): {"kind": "no_arguments",
                         "evidence": "discus/prog/discus_save.f90 save_struc: ianz>0 -> write file, RETURN"},
    ("discus", "insert"): {"kind": "first_word", "words": ["domain", "molecule", "object"], "min_abbrev": 3,
                           "evidence": "discus/prog/discus_kdo.f90 'insert' branch; otherwise do_ins (one atom)"},
    ("discus/output", "run"): {"kind": "option", "name": "trans", "value": "yes",
                               "evidence": "discus/prog/graphic.f90 'run': lib_trans_menu only if trans:yes"},
}
# Error-only early returns (the menu opens whenever the arguments are valid).
UNCONDITIONAL = {("kuplot", "fit"), ("kuplot", "ksave")}
ROOTS: dict = {}  # root dispatcher routine -> program (filled in build)
# Generic commands every menu handles the same way (documented once, at the top level / command language).
STANDARD = {"continue", "echo", "evaluate", "eval", "exit", "help", "?", "system", "syst", "wait", "reset", "@", "="}


def _routine(table: dict, name: str, directory: str) -> Routine | None:
    candidates = table.get(name, [])
    local = [r for r in candidates if r.file.startswith(directory)] or candidates
    if not local:
        return None
    # interface blocks / procedure pointers repeat the name; the real routine has the body
    return max(local, key=lambda r: (len(r.branches), len(r.calls), r.end_line - r.line))


def _reachable(table, start: Routine, directory: str, menus: set, depth=CALL_DEPTH):
    """Routines reachable from ``start`` within ``depth`` calls, not entering other menus."""
    seen, out, queue = set(), [], deque([(c, 1) for c in sorted(start.calls)])
    while queue:
        name, d = queue.popleft()
        if name in seen or name in OPAQUE:
            continue
        seen.add(name)
        r = _routine(table, name, directory)
        if r is None:
            continue
        out.append((r, d))
        if name in menus or d >= depth:
            continue
        queue.extend((c, d + 1) for c in sorted(r.calls))
    return out


def _help_topic(root: Topic, words: list[str]) -> Topic | None:
    node = root
    for w in words:
        nxt = match_topic(node.children, w.lower(), 2)
        if nxt is None:
            return None
        node = nxt
    return node


def _command_topics(menu_topic: Topic | None) -> list[Topic]:
    if menu_topic is None:
        return []
    out = []
    for c in menu_topic.children:
        if c.name.lower() in ("commands", "summary", "news"):
            continue
        out.append(c)
    return out


def build(root: Path, version: str, commit: str, label: str) -> dict:
    dirs = [LIBRARY[0]] + [p[0] for p in PROGRAMS.values()]
    table = scan_tree(root, dirs)
    menu_routines = {n for n, rs in table.items() if any(r.prompts for r in rs)}
    menu_routines |= {p[1] for p in PROGRAMS.values()} | set(FALLBACKS)
    help_roots = {prog: parse_help(root / spec[2]) for prog, spec in PROGRAMS.items()}
    help_roots["common"] = parse_help(root / LIBRARY[2])

    ROOTS.clear()
    ROOTS.update({spec[1]: prog for prog, spec in PROGRAMS.items()})
    menus: dict[str, dict] = {}
    unresolved: list[dict] = []
    visited_routines: set = set()
    canonical: dict = {}  # (routine file:line, prompt suffix) -> menu path
    queue: deque = deque()

    def command_entry(branch, routine, directory, menu_path, topics):
        callees = _reachable(table, Routine(branch.name, routine.file, branch.line,
                                            calls=set(branch.calls)), directory, menu_routines)
        keywords = dict(branch.keywords)
        optional = list(branch.optional)
        if branch.uses_get_optional:
            optional += routine.optional
        enters = []
        for r, depth in callees:
            for word, m in r.keywords.items():
                keywords.setdefault(word, m)
            for b in r.branches if r.name not in menu_routines else []:
                for word, m in b.keywords.items():
                    keywords.setdefault(word, m)
            if r.uses_get_optional:
                optional += [o for o in r.optional if o not in optional]
            # a menu is entered by a direct call from the branch; a program root is reached
            # through its loop routine (suite: 'discus' -> discus_loop -> discus_mache_kdo)
            if r.name in menu_routines and r.prompts and depth == 1:
                enters.append(r)
            elif r.name in ROOTS and ROOTS[r.name] != menu_path and depth <= 2:
                entry_program = ROOTS[r.name]
                enters.append(entry_program)
        topic = match_topic(topics, branch.name.lower(), branch.min_abbrev)
        entry = {
            "name": branch.name,
            "min_abbrev": branch.min_abbrev,
            "compare_len": branch.compare_len,
            "source": f"{routine.file}:{branch.line}",
            "keywords": dict(sorted(keywords.items())),
            "optional": list(dict.fromkeys(optional)),  # declaration order: get_optional takes the first match
        }
        if topic is not None:
            entry["help"] = {"topic": topic.name, "line": topic.line, "syntax": topic.syntax,
                             "summary": topic.summary[:400]}
            subs = [{"topic": t.name, "syntax": t.syntax} for t in topic.children
                    if t.name.lower() not in ("commands", "examples", "example")]
            if subs:
                entry["help"]["subtopics"] = subs
        return entry, enters

    def visit(routine: Routine, path: str, program: str, directory: str, help_words: list[str]):
        if path in menus:
            return
        visited_routines.add(routine.name)
        fallback = path in FALLBACKS.values()
        hw = help_words
        if routine.help_paths:
            words = [w for w in routine.help_paths[-1].split() if w]
            hw = words if words and words[0] in help_roots else [program] + words
        help_root = help_roots.get(hw[0] if hw else program)
        words = ["Command_language"] if hw == ["common"] else hw
        topic = _help_topic(help_root, words) if help_root and hw else None
        topics = _command_topics(topic)
        menu = {"program": program, "routine": f"{routine.file}:{routine.line}",
                "prompt_suffix": routine.prompts, "help_path": hw,
                "inherits": [FALLBACKS[c] for c in sorted(routine.calls) if c in FALLBACKS],
                "inherits_first": sorted(FALLBACKS[c] for c, ln in routine.fallback_calls.items()
                                         if c in FALLBACKS and routine.branches and ln < routine.branches[0].line),
                "word_length": routine.word_length,
                # commands whose lines may contain '=' without being an assignment (the indxg guard)
                "assign_exclusions": [list(e) for e in routine.assign_exclusions],
                "assignments": "do_math" in routine.calls,
                "commands": []}
        menus[path] = menu
        names = set()
        children = []
        for branch in routine.branches:
            if branch.name in names:
                continue
            names.add(branch.name)
            entry, enters = command_entry(branch, routine, directory, path, topics)
            if branch.name in STANDARD:
                entry["standard"] = True
            if enters and isinstance(enters[0], str):
                entry["enters"] = enters[0]
            elif enters:
                target = enters[0]
                rule = ENTRY_RULES.get((path, branch.name))
                if rule is not None:
                    entry["enters_when"] = rule
                elif (path, branch.name) not in UNCONDITIONAL and (
                        target.name in branch.nested_calls or target.early_return):
                    entry["enters_conditional"] = True
                suffix = next((s for s in target.prompts if s.lstrip("/").startswith(branch.name[:3])),
                              target.prompts[0])
                key = (f"{target.file}:{target.line}", suffix)
                if key in canonical:
                    entry["enters"] = canonical[key]
                else:
                    child_path = f"{program}{suffix}"
                    canonical[key] = child_path
                    entry["enters"] = child_path
                    queue.append((target, child_path, program, directory, hw + [branch.name]))
            menu["commands"].append(entry)
        documented = {c["help"]["topic"] for c in menu["commands"] if "help" in c}
        menu["help_only"] = [] if fallback else sorted(t.name for t in topics if t.name not in documented)

    for program, (directory, root_name, _hlp) in PROGRAMS.items():
        r = _routine(table, root_name, directory)
        if r is None:
            unresolved.append({"program": program, "reason": f"root routine {root_name} not found"})
            continue
        visit(r, program, program, directory, [program])
        while queue:
            visit(*queue.popleft())
    for routine_name, path in FALLBACKS.items():
        directory = LIBRARY[0] if routine_name == "kdo_all" else PROGRAMS[path.split("/")[0]][0]
        r = _routine(table, routine_name, directory)
        if r is not None:
            visit(r, path, path.split("/")[0], directory, ["common"] if path == "common" else [path.split("/")[0]])

    for name in sorted(menu_routines):
        for r in table.get(name, []):
            if r.prompts and r.name not in visited_routines:
                unresolved.append({"routine": f"{r.file}:{r.line}", "name": r.name,
                                   "prompts": r.prompts, "reason": "menu routine not reached from a program root"})

    return {
        "schema": SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "engine": {"name": "DISCUS_SUITE", "version": version, "source": "tproffen/DiffuseCode",
                   "commit": commit, "label": label},
        "generator": "tools/suite_catalogue/build_catalogue.py",
        "menus": dict(sorted(menus.items())),
        "unresolved": unresolved,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("source", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--version", required=True)
    ap.add_argument("--commit", required=True)
    ap.add_argument("--label", default="stock")
    a = ap.parse_args(argv)
    cat = build(a.source, a.version, a.commit, a.label)
    a.out.write_text(json.dumps(cat, indent=1, sort_keys=False) + "\n", newline="\n")
    n = sum(len(m["commands"]) for m in cat["menus"].values())
    print(f"{len(cat['menus'])} menus, {n} commands, {len(cat['unresolved'])} unresolved -> {a.out}")


if __name__ == "__main__":
    main()
