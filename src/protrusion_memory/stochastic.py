"""Seeded Itô Euler–Maruyama integration on an event-aligned mesh."""

from dataclasses import dataclass
import math

import numpy as np


@dataclass
class Segment:
    nodes: np.ndarray
    path: np.ndarray

    def __call__(self, times):
        index = np.searchsorted(self.nodes, times)
        if np.any(index >= len(self.nodes)) or np.any(self.nodes[index] != times):
            raise ValueError("stochastic observation time is not on the integration mesh")
        return self.path[index].T


def integrate(rhs, diffusion, initial, anchors, step, rng):
    counts = np.ceil(np.diff(anchors) / step).astype(int)
    if counts.sum() > 200_000 or counts.sum()*len(initial) > 4_000_000:
        raise ValueError("stochastic segment exceeds the step or state budget")
    nodes = np.concatenate([np.linspace(a, b, n+1)[:-1]
                            for a, b, n in zip(anchors, anchors[1:], counts)] + [anchors[-1:]])
    path = np.empty((len(nodes), len(initial)), dtype=float)
    path[0] = initial
    for i, dt in enumerate(np.diff(nodes)):
        t, y = nodes[i], path[i]
        g = diffusion(t, y)
        dz = rng.standard_normal(g.shape[1]) * math.sqrt(dt)
        path[i+1] = y + np.asarray(rhs(t, y))*dt + g @ dz
        if not np.isfinite(path[i+1]).all():
            raise ValueError("nonfinite stochastic state; check the drift, diffusion and step")
    return Segment(nodes, path)
