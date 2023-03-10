import numpy as np

class NonStiffIntegrand(object):
    
    def __init__(self,function):
        self._function = function;
    
    def Evaluate(self,x,t):
        return self._function(x,t);
