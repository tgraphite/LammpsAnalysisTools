"""lat describe — list modules or show one module's contract/help."""

from __future__ import annotations

from ..config import manual_path
from .base import Module, Param
from . import all_modules, get, register


@register
class Describe(Module):
    name = "describe"
    summary = "list available modules or describe one"
    help_doc = "manual/describe.md"
    inputs = ["-module"]
    screen_only = True
    requires_topology = False
    params = [
        Param("module", flags=("-module",), help="module name to describe"),
    ]

    def run(self, ctx, *, module):
        if module is None:
            print("lat modules:")
            for name in sorted(all_modules()):
                print(f"  {name:<12} {get(name).summary}")
            return []
        cls = get(module)
        print(f"{cls.name} — {cls.summary}")
        if cls.inputs:
            print(f"inputs:  {', '.join(cls.inputs)}")
        if cls.outputs:
            print(f"outputs: {', '.join(cls.outputs)}")
        print("parameters:")
        for p in cls.params:
            flags = ", ".join(p.flags) or p.name
            origin = f"config:{p.config}" if p.config else ("required" if p.required else "optional")
            print(f"  {flags:<16} {origin:<24} {p.help}")
        doc = manual_path(cls.help_doc)
        if doc is not None:
            print(f"\nsee: {doc}")
        return []
