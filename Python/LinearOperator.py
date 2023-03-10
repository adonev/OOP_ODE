import numpy as np
from scipy import sparse

class LinearOperator(object): # Donev: This needs a method to "Update" the matrix (for example, same sparsity pattern but different values. This method should compute factorizations of the matrix and then use that in Invert. The initializer should only allocate memory and set any internal counters etc to zero.
    
    def __init__(self,Nunknown):
        self._N = Nunknown;
        if (self._N==1):
            raise TypeError('You need to use the scalar linear operator class')
        self._L = sparse.csr_matrix(sparse.eye(self._N));
        
    def ComputeMatrixToInvert(self,alphaIdentity,alphaL): # Donev: alphaL is redundant (can always be taken care by scaling, remove it)
        self._Id = sparse.csr_matrix(sparse.eye(self._N));
        self._MatrixToInvert = alphaIdentity*self._Id+alphaL*self._L;
        
    def InvertImplicitMatrix(self,b): # Donev: The name "invert" seems wrong for a solve, same for routine above.
        return sparse.linalg.spsolve(self._MatrixToInvert,b)
    
    def Apply(self,RHS):
        return self._L.dot(RHS);
        
class ScalarLinearOperator(LinearOperator): # Donev: Remove this to not crowd the code. Assume N>1 if that is needed (not sure why???)
    def __init__(self,N):
        self._N = 1;
        if (N>1):
            raise TypeError('Cannot use scalar linear operator with > 1 unknown!')
        self._L = 1;
        
    def ComputeMatrixToInvert(self,alphaIdentity,alphaL):
        self._Id = 1;
        self._MatrixToInvert = alphaIdentity*self._Id+alphaL*self._L;
        
    def InvertImplicitMatrix(self,b):
        return b/self._MatrixToInvert;
    
    def Apply(self,RHS):
        return self._L*RHS;
