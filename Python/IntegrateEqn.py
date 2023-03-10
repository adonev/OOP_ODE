import numpy as np
from Integrator import Integrator, BackwardEuler, CrankNicolson
from NonStiffIntegrand import NonStiffIntegrand
from LinearOperator import LinearOperator, ScalarLinearOperator

"""
Solve the ODE 
dx/dt = x+e^t
"""

def NonStiffPart(x,t):
    return np.exp(t);

N=10;
dt = 1e-4;
tf=1;
x0=np.ones((N,1))
if (N==1):
    LinearPart = ScalarLinearOperator(N);
else:
    LinearPart = LinearOperator(N);
ExplicPart = NonStiffIntegrand(NonStiffPart);
TempInt = BackwardEuler(N,LinearPart,ExplicPart,dt);
nSteps = int(tf/dt+1e-6);
xFinal, _ =TempInt.Advance(x0,0,nSteps);
xTrue = np.exp(tf)*(1+tf)
print(np.amax(xFinal-xTrue))


