import os
import time
import json
import itertools
from scipy.stats import dirichlet
from pyLogVar import LogVar
from functools import reduce
import numpy as np

class Dirichlet():

    def __init__( self, alpha ):
        self.alpha = alpha
        self._dir = dirichlet( np.array( alpha ) )

    def sample( self, newAlpha=None ):
        if( newAlpha ):
            return Dirichlet.sample( newAlpha )
        return self._dir.rvs( 1 )[ 0 ]

    def pdf( self, probs ):
        return self._dir.pdf( probs )

    def logpdf( self, probs ):
        return self._dir.logpdf( probs )

class Categorical():

    def __init__( self, params=None, alpha=None ):

        if( alpha is not None ):
            self.alpha = alpha
            self._prior = Dirichlet( alpha )
            self._probs = self.resample()
        else:
            assert params is not None
            self._probs = params

        self.N = self._probs.shape[ 0 ]

    def probabilities( self ):
        return self._probs

    def resample( self, newAlpha=None, observations=None ):

        if( observations is not None ):
            # print( self.alpha )
            # print( observations )
            newAlpha = self.alpha + observations

        if( newAlpha is not None ):
            self._probs = Dirichlet( newAlpha ).sample()
        else:
            self._probs = self._prior.sample()
        return self._probs

    def sample( self, newParams=None ):
        if( newParams ):
            return Categorical.sample( newParams )
        return np.random.choice( self.N, 1, p=self._probs )[ 0 ]

    def pdf( self, i ):
        return self._probs[ i ]

    def log_pdf( self, i ):
        return np.log( self.pdf( i ) )

    def log_likelihood( self ):
        return self._prior.logpdf( self._probs )

    def log_posterior( self, observations ):
        newAlpha = self.alpha + observations
        return Dirichlet( newAlpha ).logpdf( self._probs )

    @staticmethod
    def sample( probs, normalized=False ):

        N = len( probs )
        # print('probs: '+str(probs))

        if( normalized == False ):

            total = reduce( lambda x,y: x+y, probs )
            for i in range( N ):
                probs[ i ] /= total

            # in case of passing in LogVar
            if( isinstance( probs, list ) or isinstance( probs, tuple ) ):
                probs = np.array( [ float( p ) for p in probs ] )

            if( not np.isclose( reduce( lambda x,y: x+y, probs ), 1.0 ) ):
                print( probs )
                print( reduce( lambda x,y: x+y, probs ) )
                print( reduce( lambda x,y: x+y, probs ) - 1.0 )
                assert 0

        else:
            if( not np.isclose( reduce( lambda x,y: x+y, probs ), 1.0 ) ):
                print( probs )
                print( reduce( lambda x,y: x+y, probs ) )
                print( reduce( lambda x,y: x+y, probs ) - 1.0 )
                assert 0

            if( isinstance( probs, list ) or isinstance( probs, tuple ) ):
                probs = np.array( probs )

        return np.random.choice( N, 1, p=probs )[ 0 ]

    # @staticmethod
    # def pdf( probs, i ):
    #     return probs[ i ]