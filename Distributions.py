import os
import time
import json
import itertools
from scipy.stats import dirichlet
from pyLogVar import LogVar
from functools import reduce
import numpy as np

class HDP():
    def __init__( self, gamma, alpha_0, truncation ):
        self.gamma = gamma
        self.alpha_0 = alpha_0
        self.K = truncation
        alpha = Dirchlet.sample( np.ones( self.K ) * self.gamma / self.K )
        self.beta = Dirichlet( alpha, prior=self )

    def sampleTableDishCounts( self, custDishCounts ):
        # ( Algorithm 9.3.a in Emily Fox's PhD Thesis and in pyhsmm )

        # Samples the number of tables over all of the restaurants
        # serving the different dishes ( m_j. )

        # custDishCounts is the number of customers in each
        # restaurant eating a dish ( n_j.k ).  A transition from
        # state j to state i is equivalent to having a customer
        # in restaurant j eating dish i

        # m_jk is the number of tables in rest j serving dish k
        m = np.zeros_like( custDishCounts )

        # generate a random number for each customer
        tot = custDishCounts.sum()
        randseq = np.random.random( tot )

        starts = np.empty_like( custDishCounts )
        starts[ 0, 0 ] = 0
        # for indexing into the correct customers random number
        starts.flat[ 1: ] = np.cumsum( np.ravel( custDishCounts )[ :custDishCounts.size - 1 ] )

        for ( i, j ), n in np.ndenumerate( custDishCounts ):
            w = self.beta[ j ]
            for k in range( n ):
                # sample for each customer at table j in restaurant i
                m[ i, j ] += randseq[ starts[ i, j ] + k ] \
                        < ( self.alpha_0 * w ) / ( k + self.alpha_0 * w )

        return m

    def posteriorSampleBeta( self, transitionCounts ):
        dishCounts = self.sampleTableDishCounts( transitionCounts )
        alpha = np.hstack( ( dishCounts, np.array( [ self.gamma ] ) ) )
        return Dirichlet.sample( alpha )

    def sampleProbs( self ):
        return Dirichlet.sample( self.beta * self.alpha_0 )

    def resample( self, transitionCounts ):
        self.beta.resample( transitionCounts )

class Dirichlet():

    def __init__( self, alpha, prior=None ):
        self.alpha = alpha
        self.prior = prior

    @classmethod
    def sample( cls, alpha ):
        return dirichlet.rvs( alpha, size=1 )[ 0 ]

    @classmethod
    def pdf( cls, x, alpha ):
        return dirichlet.pdf( x, alpha )

    @classmethod
    def log_pdf( cls, x, alpha ):
        return dirichlet.logpdf( x, alpha )

    def ilog_pdf( self, x ):
        return dirichlet.logpdf( x, self.alpha )

    def resample( self, x ):
        if( isinstance( self.prior, HDP ) ):
            # x should be transitionCounts
            self.alpha = self.prior.posteriorSampleBeta( x )
        else:
            assert 0, 'Not using any other prior'

class Categorical():

    def __init__( self, params=None, alpha=None ):

        if( alpha is not None ):
            self.alpha = alpha
            self.prior = Dirichlet( alpha )
            self.probs = Dirichlet.sample( alpha )
        else:
            assert params is not None
            self.probs = params

        self.N = self.probs.shape[ 0 ]

    def probabilities( self ):
        return self.probs

    def resample( self, x ):
        alpha = self.alpha + x
        self.probs = Dirichlet.sample( alpha )
        return self.probs

    @classmethod
    def pdf( cls, i, probs ):
        return probs[ i ]

    def ipdf( self, i ):
        return self.probs[ i ]

    @classmethod
    def log_pdf( cls, i, probs ):
        return np.log( cls.pdf( i ) )

    def log_likelihood( self ):
        return self.prior.ilog_pdf( self.probs )

    def log_posterior( self, observations ):
        newAlpha = self.alpha + observations
        return Dirichlet( newAlpha ).log_pdf( self.probs )

    @classmethod
    def sample( cls, probs, normalized=True ):

        N = len( probs )

        if( normalized == False ):

            total = reduce( lambda x,y: x+y, probs )
            for i in range( N ):
                probs[ i ] /= total

            # in case of passing in LogVar
            if( isinstance( probs, list ) or isinstance( probs, tuple ) ):
                probs = np.array( [ float( p ) for p in probs ] )

            if( not np.isclose( reduce( lambda x,y: x+y, probs ), 1.0 ) ):
                assert 0

        else:
            if( not np.isclose( reduce( lambda x,y: x+y, probs ), 1.0 ) ):
                assert 0

            if( isinstance( probs, list ) or isinstance( probs, tuple ) ):
                probs = np.array( probs )

        return np.random.choice( N, 1, p=probs )[ 0 ]
