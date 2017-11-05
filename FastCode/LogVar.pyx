from __future__ import division
from libc.math cimport exp,log
cimport scipy.special.cython_special as csc

ZERO = 'ZERO'
UNDERFLOW_CUTOFF = 500

cdef _logAdd(float log_a,float log_b):
    if(log_a - log_b > UNDERFLOW_CUTOFF):
        return log_a
    elif(log_b - log_a > UNDERFLOW_CUTOFF):
        return log_b

    if(log_b < log_a):
        return log_a+csc.log1p(exp(log_b-log_a))
    else:
        return log_b+csc.log1p(exp(log_a-log_b))

def logAdd(log_a, log_b):
    if(log_a == ZERO):
        return log_b
    if(log_b == ZERO):
        return log_a

    return _logAdd(<float>log_a,<float>log_b)

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

cpdef extractVal(val):
    if(isinstance(val,LogVar)):
        return val.logVal
    if(val == 0):
        return ZERO
    return log(val)

cdef class LogVar:

    cdef public object logVal

    def __init__(self,val=None,isLog=False):
        if(isLog):
            self.logVal = val
        else:
            if(isinstance(val,type(None))):
                self.logVal = ZERO
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

    def __richcmp__(self,other,op):
        otherVal = extractVal(other)
        otherValCMP = -9999999999999.9999 if otherVal == 'ZERO' else otherVal
        logValCMP = -9999999999999.9999 if self.logVal == 'ZERO' else self.logVal
        if(op == 0):
            return logValCMP < otherValCMP
        elif(op == 1):
            return logValCMP <= otherValCMP
        elif(op == 2):
            return logValCMP == otherValCMP
        elif(op == 3):
            return logValCMP != otherValCMP
        elif(op == 4):
            return logValCMP > otherValCMP
        elif(op == 5):
            return logValCMP >= otherValCMP

    def __float__(self):
        if(self.logVal == ZERO):
            return 0.
        return exp(self.logVal)

    def __int__(self):
        return int(float(self))

    def __str__(self):
        return str(float(self))

    def __repr__(self):
        return str(self)
