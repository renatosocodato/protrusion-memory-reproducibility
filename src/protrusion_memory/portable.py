"""Declarative deterministic ODEs, parameter interventions and observation maps."""

import ast
from dataclasses import dataclass
from decimal import Decimal
import keyword
import math
import operator
import re

import numpy as np
from scipy import __version__ as scipy_version
from scipy.integrate import solve_ivp

SCHEMA = "ode-v1"
RTOL, ATOL = 1e-9, 1e-11
BINARY = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
          ast.Div: operator.truediv, ast.Pow: operator.pow}
UNARY = {ast.UAdd: operator.pos, ast.USub: operator.neg}
FUNCTIONS = {"abs": abs, **{name: getattr(math, name)
                           for name in ("exp", "log", "sqrt", "sin", "cos", "tanh")}}
MODES = {"set": lambda _, value: value, "add": operator.add, "scale": operator.mul}


def number(value, context):
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError(f"{context} must be a finite real number")
    return float(value)


def table(value, context, *, fields=None):
    if not isinstance(value, dict) or not value:
        raise ValueError(f"{context} must be a nonempty table")
    if fields is not None and (unknown := value.keys() - fields):
        raise ValueError(f"{context}: unknown fields {', '.join(sorted(unknown))}")
    return value


def expression(source, symbols):
    if not isinstance(source, str) or len(source) > 1024:
        raise ValueError("expressions must be text of at most 1024 characters")
    try:
        parsed = ast.parse(source.strip(), mode="eval")
    except SyntaxError as error:
        raise ValueError(f"invalid expression: {source!r}") from error
    if sum(1 for _ in ast.walk(parsed)) > 120:
        raise ValueError("expression is too complex")

    def read(node):
        if isinstance(node, ast.Constant):
            constant = number(node.value, "expression constant")
            return lambda _: constant
        if isinstance(node, ast.Name) and node.id in symbols:
            return lambda values: values[node.id]
        if isinstance(node, ast.BinOp) and type(node.op) in BINARY:
            fn = BINARY[type(node.op)]
            left, right = read(node.left), read(node.right)
            return lambda values: fn(left(values), right(values))
        if isinstance(node, ast.UnaryOp) and type(node.op) in UNARY:
            fn, item = UNARY[type(node.op)], read(node.operand)
            return lambda values: fn(item(values))
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id in FUNCTIONS and len(node.args) == 1 and not node.keywords):
            fn, item = FUNCTIONS[node.func.id], read(node.args[0])
            return lambda values: fn(item(values))
        raise ValueError(f"unsupported expression or unknown symbol: {ast.unparse(node)}")

    function = read(parsed.body)
    return lambda values: number(function(values), f"expression {source!r}")


@dataclass(frozen=True)
class Modulation:
    parameter: str
    start_s: float
    stop_s: float
    operation: str
    value: float


def read_modulations(items, parameters, end):
    if not isinstance(items, list) or len(items) > 1000:
        raise ValueError("modulations must be a list with at most 1000 entries")
    mods = []
    for item in items:
        fields = set(Modulation.__dataclass_fields__)
        table(item, "modulation", fields=fields)
        if item.keys() != fields:
            raise ValueError("a modulation requires parameter, start_s, stop_s, operation and value")
        parameter, op = item["parameter"], item["operation"]
        if not isinstance(parameter, str) or parameter not in parameters:
            raise ValueError("modulation parameter must name a declared parameter")
        if not isinstance(op, str) or op not in MODES:
            raise ValueError("modulation operation must be set, add or scale")
        start = number(item["start_s"], "modulation start_s")
        stop = number(item["stop_s"], "modulation stop_s")
        if not 0 <= start < stop <= end:
            raise ValueError("modulation must satisfy 0 <= start_s < stop_s <= time.stop_s")
        mods.append(Modulation(parameter, start, stop, op, number(item["value"], "modulation value")))
    mods.sort(key=lambda item: (item.parameter, item.start_s))
    for left, right in zip(mods, mods[1:]):
        if left.parameter == right.parameter and right.start_s < left.stop_s:
            raise ValueError(f"overlapping modulations for {left.parameter}")
    return mods


def simulate(config):
    """Solve a portable TOML system and return observed channels in declared units."""
    table(config, "configuration",
          fields={"schema", "time", "states", "parameters", "equations", "observations", "modulations"})
    if config.get("schema") != SCHEMA:
        raise ValueError(f"expected schema {SCHEMA}")
    clock = table(config["time"], "time", fields={"stop_s", "sample_s"})
    end = number(clock["stop_s"], "stop_s")
    step = number(clock["sample_s"], "sample_s")
    if not 0 < step <= end or end / step > 100000:
        raise ValueError("require 0 < sample_s <= stop_s and at most 100001 samples")
    # Form decimal grid multiples before casting so a declared 0.9 boundary
    # matches the third step of a 0.3-second clock.
    decimal_step, decimal_end = Decimal(str(step)), Decimal(str(end))
    samples = int(decimal_end // decimal_step) + 1
    times = np.fromiter((float(index * decimal_step) for index in range(samples)), dtype=float)
    times = np.append(times[times < end], end)
    states = table(config["states"], "states")
    parameters = config.get("parameters", {})
    if not isinstance(parameters, dict) or len(states) + len(parameters) > 128:
        raise ValueError("parameters must be a table and states plus parameters must be <= 128")
    if states.keys() & parameters.keys():
        raise ValueError("states and parameter names must be distinct")
    symbols = states.keys() | parameters.keys()
    if any(not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name) or keyword.iskeyword(name)
           for name in symbols):
        raise ValueError("state and parameter names must be ASCII identifiers")
    if symbols & (FUNCTIONS.keys() | {"t"}):
        raise ValueError("t and expression function names are reserved")
    y = [number(value, name) for name, value in states.items()]
    parameters = {name: number(value, name) for name, value in parameters.items()}
    eq = table(config["equations"], "equations", fields=states.keys())
    if eq.keys() != states.keys():
        raise ValueError("each state requires exactly one derivative equation")
    symbols = symbols | {"t"}
    derivatives = [expression(eq[name], symbols) for name in states]
    modulations = read_modulations(config.get("modulations", []), parameters, end)
    channels = {}
    for name, channel in table(config["observations"], "observations").items():
        table(channel, name, fields={"expression", "unit", "operator", "noise_sd"})
        unit = channel["unit"]
        kind = channel.get("operator", "point")
        if not isinstance(unit, str) or not unit.strip() or kind not in ("point", "interval_rate"):
            raise ValueError(f"{name} requires a unit and a point or interval_rate operator")
        sd = number(channel["noise_sd"], f"{name} noise_sd") if "noise_sd" in channel else None
        if sd is not None and sd <= 0:
            raise ValueError(f"{name} noise_sd must be positive")
        channels[name] = (unit, kind, sd, expression(channel["expression"], symbols))

    def parameters_at(time):
        active = parameters.copy()
        for item in modulations:
            if item.start_s <= time < item.stop_s:
                active[item.parameter] = number(MODES[item.operation](parameters[item.parameter], item.value),
                                                f"modulated {item.parameter}")
        return active

    edges = sorted({0.0, end, *(m.start_s for m in modulations), *(m.stop_s for m in modulations)})
    traces = np.empty((len(times), len(states)))
    for left, right in zip(edges, edges[1:]):
        active = parameters_at(left)

        def rhs(time, values):
            context = {**active, "t": float(time), **dict(zip(states, map(float, values)))}
            return [equation(context) for equation in derivatives]

        solution = solve_ivp(rhs, (left, right), y, dense_output=True,
                             method="DOP853", rtol=RTOL, atol=ATOL)
        if not solution.success:
            raise ValueError(f"ODE integration failed: {solution.message}")
        selected = (times >= left) & ((times < right) if right < end else (times <= right))
        if np.any(selected):
            traces[selected] = solution.sol(times[selected]).T
        y = solution.y[:, -1]
    if not np.isfinite(traces).all():
        raise ValueError("ODE trajectory contains non-finite values")

    values = {name: [] for name in channels}
    for time, state in zip(times, traces):
        context = {**parameters_at(float(time)), "t": float(time), **dict(zip(states, map(float, state)))}
        for name, (_, _, _, measure) in channels.items():
            values[name].append(measure(context))
    observed = {}
    for name, (unit, kind, sd, _) in channels.items():
        rate = kind == "interval_rate"
        series = (np.diff(values[name]) / np.diff(times)).tolist() if rate else values[name]
        for value in series:
            number(value, f"{name} observation")
        observed[name] = {"unit": unit, "operator": kind, "noise_sd": sd,
                          "time_s": times[int(rate):].tolist(), "values": series}
    return {"schema": SCHEMA, "time_s": times.tolist(), "states": dict(zip(states, traces.T.tolist())),
            "observations": observed,
            "solver": {"method": "DOP853", "rtol": RTOL, "atol": ATOL, "scipy": scipy_version}}


def compare(left, right):
    """Compare aligned expected observations without assigning mechanism probabilities."""
    a, b = simulate(left)["observations"], simulate(right)["observations"]
    if a.keys() != b.keys():
        raise ValueError("portable comparison requires the same observed channel names")
    results = {}
    for name, x in a.items():
        y = b[name]
        if any(x[key] != y[key] for key in ("unit", "operator", "time_s", "noise_sd")):
            raise ValueError(f"{name}: observation units, operators, times and noise scales must match")
        residual = np.asarray(x["values"]) - np.asarray(y["values"])
        rmse = math.hypot(*residual) / math.sqrt(len(residual))
        sd = x["noise_sd"]
        results[name] = {"unit": x["unit"], "samples": len(residual), "rmse": number(rmse, "RMSE"),
                         "max_absolute_difference": number(float(np.max(abs(residual))), "maximum deviation"),
                         "rmse_in_noise_sd": number(rmse / sd, "scaled RMSE") if sd is not None else None}
    return {"schema": SCHEMA, "comparison": "deterministic expected observations", "channels": results}
