#!/usr/bin/env python3
"""Symbolic identities for the Lorenz-based benchmark, not a chaos certificate."""
from pathlib import Path
import json
from reproduction_io import output_directory
import sympy as sp
out=output_directory(Path(__file__))/'symbolic_checks.json'
if out.exists():raise FileExistsError('Preserve prior verification')
U,V,W,eta,r,r0,kappa,z,d=sp.symbols('U V W eta r r0 kappa z d',real=True)
t=sp.symbols('t',nonnegative=True,real=True);s=sp.symbols('s',nonnegative=True,real=True);x0=sp.symbols('x0',real=True);f=sp.Function('f')
field=sp.Matrix([10*(V-U)+eta*sp.tanh(U/10),U*(r-W)-V,U*V-sp.Rational(8,3)*W])
lorenz=sp.Matrix([10*(V-U),U*(28-W)-V,U*V-sp.Rational(8,3)*W])
red=[sp.simplify(v) for v in (field.subs({eta:0,r:28})-lorenz)]
candidate=sp.exp(-t/2)*(x0+sp.integrate(sp.exp(s/2)*f(s)/2,(s,0,t)))
res=sp.simplify(sp.diff(candidate,t)-(f(t)-candidate)/2)
init=sp.simplify(candidate.subs(t,0)-x0)
rshift=sp.simplify(((r0-kappa*d)+kappa*(z+d))-(r0+kappa*z))
jac=field.jacobian([U,V,W]); lower=sp.simplify(sp.diff((sp.Symbol('m')-sp.Symbol('x'))/2,sp.Symbol('x'))+sp.Rational(1,2))
summary={'model_eta_zero_lorenz_difference':[str(v) for v in red],'filter_ode_residual':str(res),'filter_initial_residual':str(init),'ensemble_parameter_symmetry_residual':str(rshift),'fiber_rate_residual':str(lower),'base_jacobian_trace':str(sp.simplify(sp.trace(jac))),'identity_checks_all_zero':bool(all(v==0 for v in red) and res==0 and init==0 and rshift==0 and lower==0),'formal_closure':'UNDECIDED','tucker_theorem':'external source, not re-proved here','modified_eta_0_5':'prior floating point simulations only'}
assert summary['identity_checks_all_zero']
out.write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
