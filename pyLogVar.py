from __future__ import division
import numpy as np
ZERO = 'ZERO'
UNDERFLOW = 500

def _logAdd( logA, logB ):
    # compute log( A + B ) using log( A ) and log( B )
    return logA + np.log1p( np.exp( logB - logA ) )

def _logSub( logA, logB ):
    # compute log( A - B ) using log( A ) and log( B )
    return logA + np.log1p( -np.exp( logB - logA ) )

def logAdd( logA, sgnA, logB, sgnB ):
    if( logA == ZERO ):
        return ( logB, sgnB )
    if( logB == ZERO ):
        return ( logA, sgnA )

    if( logA - logB > UNDERFLOW ):
        # B is wayyyy smaller than A and is basically 0
        return ( logA, sgnA )
    elif( logB - logA > UNDERFLOW ):
        # A is wayyyy smaller than B and is basically 0
        return ( logB, sgnB )

    if( sgnA == sgnB ):
        # either computing ( a + b ) or -( a + b )

        if( logB < logA ):
            # make sure that the # that goes into np.exp isn't too big
            return _logAdd( logA, logB ), sgnA
        else:
            return _logAdd( logB, logA ), sgnA

    else:
        # either computing ( a - b ) or -( a - b )
        if( sgnA == 1 ):
            # computing a - b
            if( logB > logA ):
                return ( _logSub( logB, logA ), -1 )
            else:
                return ( _logSub( logA, logB ), 1 )
        else:
            # computing b - a
            if( logB < logA ):
                return ( _logSub( logB, logA ), -1 )
            else:
                return ( _logSub( logA, logB ), 1 )

def logMul( logA, sgnA, logB, sgnB ):
    if( logA == ZERO or logB == ZERO ):
        return ( ZERO, 1 )
    if( sgnA == sgnB ):
        return ( logA + logB, 1 )
    else:
        return ( logA + logB, -1 )

def logDiv( logA, sgnA, logB, sgnB ):
    if( logA == ZERO ):
        return ( ZERO, 1 )
    if( logB == ZERO ):
        assert 0
    if( sgnA == sgnB ):
        return ( logA - logB, 1 )
    else:
        return ( logA - logB, -1 )

def extractVal( val, sgn=None ):

    if( isinstance( val, LogVar ) ):
        return ( val.logVal, val.sgn )

    if( val == 0 ):
        return ( ZERO, 1 )

    if( sgn is not None ):
        return ( np.log( val ), sgn )

    if( val > 0 ):
        return ( np.log( val ), 1 )

    return ( np.log( val ), -1 )

class LogVar():
    def __init__( self, val=None, sgn=None, isLog=False ):
        if( val == None ):
            self.logVal = ZERO
            self.sgn = 1
        elif( isLog ):
            self.logVal = val
            if( sgn is not None ):
                self.sgn = sgn
            else:
                self.sgn = 1
        else:
            self.logVal, self.sgn = extractVal( val, sgn )

    def genericOperator( self, other, func, otherNeg=False ):
        otherVal, sgn = extractVal( other )
        if( otherNeg ):
            sgn *= -1
        val, valSgn = func( self.logVal, self.sgn, otherVal, sgn )
        return LogVar( val=val, sgn=valSgn, isLog=True )

    def iGenericOperator( self, other, func, otherNeg=False ):
        otherVal, sgn = extractVal( other )
        if( otherNeg ):
            sgn *= -1
        self.logVal, self.sgn = func( self.logVal, self.sgn, otherVal, sgn )
        return self

    """ Addition """

    def __add__( self, other ):
        return self.genericOperator( other, logAdd )

    def __radd__( self, other ):
        return self+other

    def __iadd__( self, other ):
        return self.iGenericOperator( other, logAdd )

    """ Subtraction """

    def __sub__( self, other ):
        return self.genericOperator( other, logAdd, otherNeg=True )

    def __rsub__( self, other ):
        return self+other

    def __isub__( self, other ):
        return self.iGenericOperator( other, logAdd, otherNeg=True )

    """ Multiplication """

    def __mul__( self, other ):
        return self.genericOperator( other, logMul )

    def __rmul__( self, other ):
        return self*other

    def __imul__( self, other ):
        return self.iGenericOperator( other, logMul )

    """ Division """

    def __truediv__( self, other ):
        return self.genericOperator( other, logDiv )

    def __idiv__( self, other ):
        return self.iGenericOperator( other, logDiv )

    """ Power """

    def __pow__( self, exponent ):
        return LogVar( val=self.logVal * float( exponent ), isLog=True )

    def toFloat( self ):
        if( self.logVal == ZERO ):
            return 0
        return np.exp( self.logVal )

    def isZero( self ):
        return self.logVal == 'ZERO'

    def __float__( self ):
        if( self.logVal == ZERO ):
            return 0.
        return self.sgn * np.exp( self.logVal )

    def __int__( self ):
        return int( float( self ) )

    def __str__( self ):
        return str( float( self ) )

    def __repr__( self ):
        return str( self )

def test():
    val1 = LogVar( val=-100, isLog=True )
    val2 = LogVar( val=-1000, isLog=True )
    for i in range( 100000 ):
        val1 += val2
        print( val1.logVal )

# a = LogVar( 0 )
# a += 4
# a -= 2
# a *= 1.3
# a /= 59
# a += 3
# a -= 1.2
# b = LogVar( 0.12 )
# a -= b
# print( a )