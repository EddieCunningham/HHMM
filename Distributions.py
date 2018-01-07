import os
import time
import json
import itertools
import scipy.stats
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

    def log_pdf( self, probs ):
        return self._dir.logpdf( probs )

    @staticmethod
    def sample( alpha ):
        return dirichlet( np.array( alpha ) ).rvs( 1 )[ 0 ]

    @staticmethod
    def pdf( alpha, probs ):
        return dirichlet( np.array( alpha ) ).pdf( probs )

class Categorical():

    def __init__( self, params, alpha=None ):

        if( alpha ):
            self.alpha = alpha
            self._prior = Dirichlet( alpha )
            self._probs = self.resample()
        else:
            self._probs = params

        self.N = _probs.shape[ 0 ]

    def resample( self, newAlpha=None, observations=None ):

        if( observations ):
            newAlpha = self.alpha + observations

        if( newAlpha ):
            self._probs = Categorical.resample( newAlpha )
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
        return self._dir.logpdf( self._probs )

    @staticmethod
    def sample( probs, normalized=False ):

        N = len( probs )

        if( normalized == False ):

            total = reduce( lambda x,y: x+y, probs )
            for i in range( N ):
                probs[ i ] /= total

            # in case of passing in LogVar
            if( isinstance( probs, list ) or isinstance( probs, tuple ) ):
                probs = np.array( [ float( p ) for p in probs ] )

        else:
            assert reduce( lambda x,y: x+y, probs ) == 1.0

            if( isinstance( probs, list ) or isinstance( probs, tuple ) ):
                probs = np.array( probs )

        return np.random.choice( N, 1, p=probs )[ 0 ]

    @staticmethod
    def pdf( probs, i ):
        return probs[ i ]