import numpy as np
from pyLogVar import LogVar

def prettyPrint( d, indent=0 ):
    for key, value in d.items():
        print( '\t' * indent + str( key ) )
        if isinstance( value, dict ):
            prettyPrint( value, indent+1 )
        else:
            print( '\t' * ( indent+1 ) + str( value ) )


class RunningStats():

    def __init__( self, checkpoint=-1, useLogVar=False ):

        self.useLV = useLogVar
        self.checkpoint = checkpoint
        self.savedVals = []

        self.n = 0
        if( useLogVar ):
            self.M1 = LogVar( 0 )
            self.M2 = LogVar( 0 )
            self.M3 = LogVar( 0 )
            self.M4 = LogVar( 0 )
        else:
            self.M1 = 0
            self.M2 = 0
            self.M3 = 0
            self.M4 = 0

    # taken from https://www.johndcook.com/blog/skewness_kurtosis/
    def pushVal( self, x, isLog=False ):

        if( isLog ):
            if( self.useLV ):
                x = LogVar( val=x, isLog=True )
            else:
                x = np.exp( x )

        n1 = self.n
        self.n += 1
        delta = x - self.M1
        delta_n = delta / self.n
        delta_n2 = delta_n * delta_n
        term1 = delta * delta_n * n1
        self.M1 += delta_n
        self.M4 += term1 * delta_n2 * ( self.n**2 - 3*self.n + 3 ) + 6 * delta_n2 * self.M2 - 4 * delta_n * self.M3
        self.M3 += term1 * delta_n * ( self.n - 2 ) - 3 * delta_n * self.M2
        self.M2 += term1

        if( self.checkpoint != -1 and self.n % self.checkpoint == 0 ):
            if( self.useLV ):
                self.saveLogVals()
            else:
                self.saveVals()

    def _mean( self ):
        return self.M1

    def mean( self ):
        return float( self._mean() )

    def log_mean( self ):
        if( self.useLV ):
            return self._mean().logValue()
        return np.log( self._mean() )

    def _variance( self ):
        if( self.n > 1 ):
            return self.M2 / ( self.n - 1.0 )
        else:
            if( self.useLV ):
                return LogVar( 0 )
            return 0.0

    def variance( self ):
        return float( self._variance() )

    def log_variance( self ):
        if( self.useLV ):
            return self._variance().logValue()
        return np.log( self._variance() )

    def _skewness( self ):
        return np.sqrt( self.n ) * self.M3 / self.M2**1.5

    def skewness( self ):
        return float( self._skewness() )

    def log_skewness( self ):
        if( self.useLV ):
            return self._skewness().logValue()
        return np.log( self._skewness() )

    def _kurtosis( self ):
        return self.n * self.M4 / self.M2**2 - 3.0

    def kurtosis( self ):
        return float( self._kurtosis() )

    def log_kurtosis( self ):
        if( self.useLV ):
            return self._kurtosis().logValue()
        return np.log( self._kurtosis() )

    def saveVals( self ):
        self.savedVals.append( [ mean(), variance(), skewness(), kurtosis() ] )

    def saveLogVals( self ):
        self.savedVals.append( [ log_mean(), log_variance(), log_skewness(), log_kurtosis() ] )

    def getSavedValues( self ):
        times = self.n * np.arange( len( self.savedVals ) )
        return np.hstack( ( times, np.array( self.savedVals ) ) )


def test():
    a = RunningStats(useLogVar=True)
    for i in range( 100 ):
        a.pushVal( i )

    print( a.mean() )
    print( np.sqrt( a.variance() ) )
    print( a.variance() )
    print( a.skewness() )
    print( a.kurtosis() )

test()