"""Static selection semantics (MDAnalysis-like subset) + ndx exchange.

Grammar: expr := and ('or' and)* ; and := not ('and' not | not)* ;
not := 'not' not | primary ; primary := '(' expr ')' | keyword values.
Keywords: all, type, element, name, resname, resid, mol, index, id.
Adjacent keywords imply 'and'. No interactive input.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np

from .errors import LatError
from .topology import Topology

_KEYWORDS = {"all", "type", "element", "name", "resname", "resid", "mol", "index", "id"}
_OPS = {"and", "or", "not"}


@dataclass
class Selection:
    indices: np.ndarray
    literal: str

    @property
    def n_atoms(self) -> int:
        return int(len(self.indices))


def _tokenize(text: str) -> list[str]:
    return re.findall(r"[()]|[^\s()]+", text)


class _Parser:
    def __init__(self, topo: Topology, tokens: list[str], ndx: dict | None):
        self.topo = topo
        self.tokens = tokens
        self.pos = 0
        self.ndx = ndx or {}

    def peek(self):
        return self.tokens[self.pos] if self.pos < len(self.tokens) else None

    def next(self):
        tok = self.peek()
        self.pos += 1
        return tok

    def parse(self) -> np.ndarray:
        mask = self.parse_or()
        if self.peek() is not None:
            raise LatError(f"unexpected token {self.peek()!r} in selection")
        return mask

    def parse_or(self):
        mask = self.parse_and()
        while self.peek() == "or":
            self.next()
            mask = mask | self.parse_and()
        return mask

    def parse_and(self):
        mask = self.parse_not()
        while True:
            tok = self.peek()
            if tok == "and":
                self.next()
                mask = mask & self.parse_not()
            elif tok is not None and tok not in _OPS and tok not in (")",):
                # implicit 'and' before the next keyword
                mask = mask & self.parse_not()
            else:
                break
        return mask

    def parse_not(self):
        if self.peek() == "not":
            self.next()
            return ~self.parse_not()
        return self.parse_primary()

    def parse_primary(self):
        tok = self.next()
        if tok is None:
            raise LatError("unexpected end of selection")
        if tok == "(":
            mask = self.parse_or()
            if self.next() != ")":
                raise LatError("unbalanced parentheses in selection")
            return mask
        if tok in _OPS or tok == ")":
            raise LatError(f"unexpected token {tok!r} in selection")
        if tok not in _KEYWORDS:
            raise LatError(f"unknown selection keyword {tok!r}")
        values = []
        while True:
            nxt = self.peek()
            if nxt is None or nxt in _OPS or nxt in ("(", ")"):
                break
            if nxt in _KEYWORDS:
                break
            values.append(self.next())
        return self.resolve(tok, values)

    def resolve(self, keyword: str, values: list[str]) -> np.ndarray:
        topo = self.topo
        n = topo.n_atoms
        if keyword == "all":
            return np.ones(n, dtype=bool)

        if not values:
            raise LatError(f"selection keyword {keyword!r} needs at least one value")

        # ndx group reference: keyword values may name a group
        if keyword in ("index", "id") and len(values) == 1 and values[0] in self.ndx:
            return _mask_from_indices(n, self.ndx[values[0]])

        if keyword == "type":
            ts = _parse_int_tokens(values)
            return topo.mask_of_types(ts)
        if keyword in ("index", "id", "mol", "resid"):
            if keyword == "mol" and topo.mols is None:
                raise LatError("selection 'mol' needs molecule IDs (missing in topology)")
            if keyword == "resid" and topo.resids is None:
                raise LatError("selection 'resid' unavailable in this topology; use type/mol")
            vals = set(_parse_int_tokens(values))
            if keyword == "index":
                return _mask_from_indices(n, [v for v in vals if 0 <= v < n])
            if keyword == "id":
                return np.isin(topo.atom_id, list(vals))
            if keyword == "mol":
                return np.isin(topo.mols, list(vals))
            return np.isin(topo.resids, list(vals))
        if keyword == "element":
            if not topo.type_element:
                raise LatError("selection 'element' needs type->element mapping")
            wanted = {v.capitalize() for v in values}
            types = [t for t, e in topo.type_element.items() if e in wanted]
            return topo.mask_of_types(types) if types else np.zeros(n, dtype=bool)
        if keyword == "name":
            if topo.names is None:
                raise LatError("selection 'name' unavailable in this topology; use type/mol")
            return np.isin(topo.names, values)
        if keyword == "resname":
            if topo.resnames is None:
                raise LatError("selection 'resname' unavailable in this topology; use type/mol")
            return np.isin(topo.resnames, values)
        raise LatError(f"unsupported selection keyword {keyword!r}")


def _parse_int_tokens(values: list[str]) -> list[int]:
    out: list[int] = []
    for tok in values:
        for piece in tok.split(","):
            if not piece:
                continue
            if "-" in piece[1:]:
                a, b = piece.split("-", 1)
                out.extend(range(int(a), int(b) + 1))
            else:
                out.append(int(piece))
    return out


def _mask_from_indices(n: int, indices) -> np.ndarray:
    mask = np.zeros(n, dtype=bool)
    idx = np.asarray(list(indices), dtype=int)
    idx = idx[(idx >= 0) & (idx < n)]
    if len(idx):
        mask[idx] = True
    return mask


def parse(topo: Topology, literal: str | None, ndx: dict | None = None) -> Selection:
    if literal is None:
        raise LatError("selection missing (use -sel/-sel1/-sel2)")
    text = literal.strip()
    if literal in (ndx or {}):
        return Selection(_mask_from_indices(topo.n_atoms, ndx[literal]), f"ndx:{literal}")
    tokens = _tokenize(text)
    if not tokens:
        raise LatError("empty selection")
    mask = _Parser(topo, tokens, ndx).parse()
    return Selection(np.flatnonzero(mask), text)


# ---- ndx exchange ------------------------------------------------------
def read_ndx(path) -> dict[str, list[int]]:
    groups: dict[str, list[int]] = {}
    current = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith("[") and line.endswith("]"):
                current = line[1:-1].strip()
                groups[current] = []
            elif current is not None:
                groups[current].extend(int(x) for x in line.split())
    return groups


def write_ndx(path, groups: dict[str, list[int]]) -> None:
    lines = []
    for name, idx in groups.items():
        lines.append(f"[ {name} ]")
        idx = list(idx)
        for i in range(0, len(idx), 15):
            lines.append(" ".join(str(x) for x in idx[i:i + 15]))
        lines.append("")
    with open(path, "w") as fh:
        fh.write("\n".join(lines))
