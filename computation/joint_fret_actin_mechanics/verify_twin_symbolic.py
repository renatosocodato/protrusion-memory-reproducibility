#!/usr/bin/env python3
"""Verify the frozen out-of-family change of coordinates for arbitrary inputs."""
from pathlib import Path
from reproduction_io import output_directory
import sympy as s
import json
q,a,u,v=s.symbols('q a u v', real=True);qp=s.Rational(8,25)*q+s.Rational(17,25)*a
qd=(u-q)/8;ad=(q-a+v)/25
res_q=s.factor(s.diff(qp,q)*qd+s.diff(qp,a)*ad-(u+s.Rational(17,25)*v-qp)/25)
res_a=s.factor(ad-(qp-a+s.Rational(8,25)*v)/8)
if res_q!=0 or res_a!=0:raise ValueError('Frozen all-input identity failed')
r={'check':'exact rational change-of-coordinates','qprime_residual':str(res_q),'actin_residual':str(res_a),'scope':'all u,v in the declared linear ODE; no physical validity asserted','sympy':s.__version__}
(output_directory(Path(__file__)) / 'results_v1').joinpath('symbolic_twin.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
