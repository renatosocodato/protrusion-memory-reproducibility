from pathlib import Path
from reproduction_io import output_directory
import json,sympy as s
out=output_directory(Path(__file__))/'symmetry_checks.json'
if out.exists():raise FileExistsError('Preserve prior verification')
u,v,w,x,r,eta=s.symbols('u v w x r eta',real=True)
M=s.diag(-1,-1,1,1)
y=s.Matrix([u,v,w,x]);p=1/(1+s.exp(-(w-25)/30))
F=s.Matrix([10*(v-u)+eta*s.tanh(u/10),u*(r-w)-v,u*v-s.Rational(8,3)*w,(p-x)/2]);subs={u:-u,v:-v,w:w,x:x}
equiv=[s.simplify(z) for z in (F.subs(subs,simultaneous=True)-M*F)]
involution=[s.simplify(z) for z in (M*M*y-y)]
q=lambda a:1/(1+s.exp(-a/20))
qsum=s.simplify(q(u)+q(-u)-1)
res={'equivariance_residual':list(map(str,equiv)),'involution_residual':list(map(str,involution)),'gate_complement_residual':str(qsum),'all_exact_zero':all(t==0 for t in equiv+involution+[qsum]),'scope':'all real r and eta, exact declared model; no cellular claim'}
out.write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
