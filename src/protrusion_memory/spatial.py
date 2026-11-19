"""Linear reaction–diffusion fields on a fixed periodic angular grid."""

import math

import numpy as np

from . import portable

SCHEMA = "ring-v1"


def integer(value, label, low, high):
    if type(value) is not int or not low <= value <= high:
        raise ValueError(f"{label} must be an integer from {low} to {high}")
    return value


def profile(spec, modes, label):
    if not isinstance(spec, dict) or set(spec) - {str(mode) for mode in modes}:
        raise ValueError(f"{label} must map declared mode numbers to coefficients")
    return {int(mode): portable.number(value, label) for mode, value in spec.items()}


def compile_model(config):
    portable.table(config, "ring configuration",
                   fields={"schema", "ring", "time", "parameters", "modulations", "fields",
                           "inputs", "reporters", "observations"})
    if config.get("schema") != SCHEMA:
        raise ValueError(f"expected {SCHEMA}")
    ring = portable.table(config["ring"], "ring", fields={"nodes", "modes"})
    nodes = integer(ring["nodes"], "nodes", 8, 1024)
    raw_modes = ring["modes"]
    if not isinstance(raw_modes, list) or not raw_modes:
        raise ValueError("ring.modes must be a nonempty list")
    modes = [integer(mode, "mode", 0, (nodes-1)//2) for mode in raw_modes]
    if len(set(modes)) != len(modes):
        raise ValueError("ring.modes must be distinct")
    modes.sort()
    components = [(n, kind) for n in modes for kind in ("cos", "sin") if n or kind == "cos"]
    fields = portable.table(config["fields"], "fields")
    params = config.get("parameters", {})
    if not isinstance(params, dict):
        raise ValueError("ring parameters must be a table")
    if len(components)*len(fields)+len(params) > 128:
        raise ValueError("compiled ring exceeds the 128-state-and-parameter budget")

    def coefficient(value):
        if isinstance(value, str) and value in params:
            return value
        return str(portable.number(value, "coupling coefficient"))

    inputs = {}
    source_tables = config.get("inputs", {})
    if not isinstance(source_tables, dict):
        raise ValueError("ring inputs must be a table")
    for name, spec in source_tables.items():
        portable.table(spec, name, fields={"parameter", "cos", "sin"})
        param = spec["parameter"]
        if not isinstance(param, str) or param not in params:
            raise ValueError(f"{name}: input must name a declared parameter")
        inputs[name] = param, {
            "cos": profile(spec.get("cos", {}), modes, name),
            "sin": profile(spec.get("sin", {}), set(modes)-{0}, name),
        }
    names = {(name, n, kind): f"field{ix}_{n}{kind}"
             for ix, name in enumerate(fields) for n, kind in components}
    initial, equations = {}, {}
    for field, spec in fields.items():
        portable.table(spec, field, fields={"decay_s_inv", "diffusion_rad2_per_s", "coupling", "drives",
                                            "initial_cos", "initial_sin"})
        decay = portable.number(spec.get("decay_s_inv", 0), "decay_s_inv")
        diffusion = portable.number(spec.get("diffusion_rad2_per_s", 0), "diffusion_rad2_per_s")
        if min(decay, diffusion) < 0:
            raise ValueError("decay and diffusivity must be nonnegative")
        coupled, driven = spec.get("coupling", {}), spec.get("drives", {})
        if (not isinstance(coupled, dict) or coupled.keys()-fields.keys()
                or not isinstance(driven, dict) or driven.keys()-inputs.keys()):
            raise ValueError("couplings and drives must name declared fields and inputs")
        start = {k: profile(spec.get(f"initial_{k}", {}), modes if k == "cos" else set(modes)-{0},
                            f"{field} initial_{k}") for k in ("cos", "sin")}
        for n, kind in components:
            state = names[field, n, kind]
            eigenvalue = 4*math.sin(math.pi*n/nodes)**2 / (2*math.pi/nodes)**2
            terms = [f"{-(decay+diffusion*eigenvalue)}*{state}"]
            terms.extend(f"({coefficient(value)})*{names[source, n, kind]}"
                         for source, value in coupled.items())
            for source, value in driven.items():
                param, patterns = inputs[source]
                amplitude = patterns[kind].get(n, 0)
                if amplitude:
                    terms.append(f"({coefficient(value)})*{amplitude}*{param}")
            initial[state] = start[kind].get(n, 0.)
            equations[state] = " + ".join(terms)

    reporters = {}
    for name, spec in portable.table(config["reporters"], "reporters").items():
        portable.table(spec, name, fields={"field", "mode", "component", "response_s", "initial"})
        key = spec["field"], spec["mode"], spec.get("component", "cos")
        if type(spec["mode"]) is not int or key not in names:
            raise ValueError(f"{name}: reporter requires a declared field and mode component")
        reporters[name] = {**{k: v for k, v in spec.items() if k in ("response_s", "initial")},
                           "expression": names[key]}
    ode = {"schema": portable.ACQUISITION_SCHEMA, "time": config["time"], "states": initial,
           "parameters": params, "equations": equations, "modulations": config.get("modulations", []),
           "reporters": reporters, "observations": config["observations"]}
    return ode, names, nodes


def simulate(config):
    ode, names, nodes = compile_model(config)
    result = portable.simulate(ode)
    states = result.pop("states")
    coefficients = {}
    for (name, n, kind), key in names.items():
        coefficients.setdefault(name, {}).setdefault(str(n), {})[kind] = states[key]
    result.update(schema=SCHEMA, nodes=nodes, modal_coefficients=coefficients)
    return result


def compare(left, right):
    a, _, an = compile_model(left)
    b, _, bn = compile_model(right)
    if an != bn:
        raise ValueError("ring comparison requires the same grid size")
    if left["reporters"] != right["reporters"] or left["observations"] != right["observations"]:
        raise ValueError("ring comparison requires the same reporter and acquisition definitions")
    result = portable.compare(a, b)
    result["schema"] = SCHEMA
    return result
