import numpy as np

class Integrator(object):

    """
    This temporal integrator solves an ODE of the form
    x'(t) = L*x(t)+f(x,t)
    """

    def __init__(self,Nunknowns,LinearPart,NonlinearPart,dt):
        self._N = Nunknowns;
        self._LinearPart = LinearPart;
        self._NonLinearPart = NonlinearPart;
        self._dt = dt;
        self._implicitCoeff = 0;
        
    def InvertImplicitPart(self,b): # Donev: Invert bad word
        return b*self._dt;
        
    def FactorImplicitPart(self): # Donev: "factor" is too specific. This stuff should go into LinearOperator's update method
        """
        The matrix to invert can be written as 
        alphaId*I+alphaL*L
        Construct and factor this matrix
        """
        
    def TArgumentForExplicit(self,t):
        return t;
    
    def XArgumentForExplicit(self,x,xprev):
        return x;
  
    def Advance(self,x0,StartNum,EndNum,xprev=None): # Donev: I would make the argument be nSteps. Internally, the method should keep track of the time step number (i.e., increment by 1 inside the loop)
        x = x0;
        if (xprev is None):
            xprev = x0;
        for iT in range(StartNum,EndNum):
            t = self._dt*iT;
            RHS = self._NonLinearPart.Evaluate(self.XArgumentForExplicit(x,xprev),self.TArgumentForExplicit(t));
            RHS += (1-self._implicitCoeff)*self._LinearPart.Apply(x);
            RHS += x/self._dt;
            xprev = x;
            x = self.InvertImplicitPart(RHS);
        return x, xprev;

class ImplicitOneStepIntegrator(Integrator):

    def FactorImplicitPart(self):
        alphaId = 1/self._dt;
        alphaL = -self._implicitCoeff;
        self._LinearPart.ComputeMatrixToInvert(alphaId,alphaL);
    
    def InvertImplicitPart(self,b):
        return self._LinearPart.InvertImplicitMatrix(b);
            

class BackwardEuler(ImplicitOneStepIntegrator):

    def __init__(self,Nunknowns,LinearPart,NonlinearPart,dt):
        super().__init__(Nunknowns,LinearPart,NonlinearPart,dt);
        self._implicitCoeff=1;
        self.FactorImplicitPart();
        
class CrankNicolson(ImplicitOneStepIntegrator):

    def __init__(self,Nunknowns,LinearPart,NonlinearPart,dt):
        super().__init__(Nunknowns,LinearPart,NonlinearPart,dt);
        self._implicitCoeff=0.5;
        self.FactorImplicitPart();
    
    def TArgumentForExplicit(self,t):
        return t+self._dt/2;
    
    def XArgumentForExplicit(self,x,xprev):
        return 1.5*x-0.5*xprev;
            
