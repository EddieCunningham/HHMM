from __future__ import division
import numpy as np
ZERO = 'ZERO'
UNDERFLOW = 500

def logAdd(log_a,log_b):

    if(log_a == ZERO):
        return log_b
    if(log_b == ZERO):
        return log_a

    if(log_a - log_b > UNDERFLOW):
        return log_a
    elif(log_b - log_a > UNDERFLOW):
        return log_b

    if(log_b < log_a):
        return log_a+np.log1p(np.exp(log_b-log_a))
    else:
        return log_b+np.log1p(np.exp(log_a-log_b))

def logMul(log_a,log_b):
    if(log_a == ZERO or log_b == ZERO):
        return ZERO
    return log_a + log_b

def logDiv(log_a,log_b):
    if(log_a == ZERO):
        return ZERO
    if(log_b == ZERO):
        assert 0
    return log_a - log_b

def extractVal(val):
    if('LogVar' in str(type(val)) or 'instance' in str(type(val))):
        return val.logVal
    if(val == 0):
        return ZERO
    return np.log(val)

class LogVar():
    def __init__(self,val=None,isLog=False):
        if(val == None):
            self.logVal = ZERO
        elif(isLog):
            self.logVal = val
        else:
            self.logVal = extractVal(val)

    def __add__(self,other):
        otherVal = extractVal(other)
        return LogVar(logAdd(self.logVal,otherVal),isLog=True)

    def __radd__(self,other):
        return self+other

    def __iadd__(self,other):
        otherVal = extractVal(other)
        self.logVal = logAdd(self.logVal,otherVal)
        return self

    def __mul__(self,other):
        otherVal = extractVal(other)
        return LogVar(logMul(self.logVal,otherVal),isLog=True)

    def __rmul__(self,other):
        return self*other

    def __imul__(self,other):
        otherVal = extractVal(other)
        self.logVal = logMul(self.logVal,otherVal)
        return self

    def __truediv__(self,other):
        otherVal = extractVal(other)
        return LogVar(logDiv(self.logVal,otherVal),isLog=True)

    def __idiv__(self,other):
        otherVal = extractVal(other)
        self.logVal = logDiv(self.logVal,otherVal)
        return self

    def toFloat(self):
        if(self.logVal == ZERO):
            return 0
        return np.exp(self.logVal)

    def isZero(self):
        return self.logVal == 'ZERO'

    def __float__(self):
        if(self.logVal == ZERO):
            return 0.
        return np.exp(self.logVal)

    def __int__(self):
        return int(float(self))

    def __str__(self):
        return str(float(self))

    def __repr__(self):
        return str(self)

def test():
    val1 = LogVar(val=-100,isLog=True)
    val2 = LogVar(val=-1000,isLog=True)
    for i in range(100000):
        val1 += val2
        print(val1.logVal)
