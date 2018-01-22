from HHMMMessagePasser import HiddenMarkovModelMessagePasser
from HHMMNode import NodeForHMM
from HHMMHG import MessagePassingHG
from Distributions import Categorical
import numpy as np
from cycleDetector import identifyCycles
from util import RunningStats
import sys

class HHMMModelBase():

    def __init__( self, graphs, transHypers, emissionHypers, rootHypers ):

        self.graphs = graphs

        self._transHypers = transHypers
        self.transDists = []

        self._emissionHypers = emissionHypers
        self.emissionDists = []

        # assume all roots have same hyper parameters
        self._rootHypers = rootHypers
        self.rootDists = []

        self.initializeGraphParameters()

        # preprocess graphs
        for graph in self.graphs:

            feedbackSet, blockManager = identifyCycles( graph )
            feedbackSetIds = [ node._id for node in feedbackSet ]
            graph.preprocess( feedbackSetIds )
            graph.draw()

    def transitionPosterior( self ):
        # Return P( transitionDist | X=x )
        # The posterior is a dirichlet distribution

        obs = self.countTransitions()

        tensor = []
        for i, matDist in enumerate( self.transDists ):
            matrix = []
            for j, transDis in enumerate( matDist ):
                matrix.append( self.transDists[ i ][ j ].alpha )
            tensor.append( matrix )

        return obs + np.array( tensor )

    def emissionPosterior( self ):
        # Return P( emissionDist | X=x )
        # The posterior is a dirichlet distribution

        obs = self.countEmissions()

        matrix = []
        for emissionDist in self.emissionDists:
            matrix.append( self.emissionDists[ i ].alpha )

        return obs + np.array( matrix )

    def transitionTensor( self ):

        tensor = []
        for i, matDist in enumerate( self.transDists ):
            matrix = []
            for j, transDis in enumerate( matDist ):
                matrix.append( self.transDists[ i ][ j ].probabilities() )
            tensor.append( matrix )
        return np.array( tensor )

    def emissionMatrix( self ):

        matrix = []
        for emissionDist in self.emissionDists:
            matrix.append( self.emissionDists[ i ].probabilities() )
        return np.array( matrix )

    # Closure so that to the graph the root distributions
    # look the same way as the transition and emission
    # pdf functions
    def _transPdfClosure( self, graph ):
        assert 0, 'Implement this in base class'

        def transPdf( parents, child, X, i ):
            assert parents[ 0 ].sex == 'female'
            return self.transDists[ X[ 0 ] ][ X[ 1 ] ].pdf( i )

        return transPdf

    def _emissionPdfClosure( self, graph ):
        assert 0, 'Implement this in base class'
        def emissionPdf( person, i ):
            return self.emissionDists[ i ].pdf( person.y )

        return emissionPdf

    def _rootPdfClosure( self, graph ):
        assert 0, 'Implement this in base class'
        roots = sorted( graph.roots )
        graphIndex = self.graphs.index( graph )

        def rootPdf( person, i ):
            personIndex = roots.index( person )
            return self.rootDists[ graphIndex ][ personIndex ].pdf( i )

        return rootPdf

    def initializeGraphParameters( self ):

        # make the transition parameters
        for matHyper in self._transHypers:
            transDist = []
            for transHyper in matHyper:
                transDist.append( Categorical( alpha=transHyper ) )
            self.transDists.append( transDist )

        # make the emission parameters
        for emissionHyper in self._emissionHypers:
            self.emissionDists.append( Categorical( alpha=emissionHyper ) )

        for graph in self.graphs:

            # make the root distribution parameters
            rootDist = [ Categorical( alpha=self._rootHypers ) for root in sorted( graph.roots ) ]
            self.rootDists.append( rootDist )

            transitionPdf = self._transPdfClosure( graph )
            emissionPdf   = self._emissionPdfClosure( graph )
            rootPdf       = self._rootPdfClosure( graph )
            graph.setParameters( transitionPdf, emissionPdf, rootPdf )

    def useFakeY( self ):
        # use each person's fakeY so we don't override y
        def emissionPdf( person, i ):
            return self.emissionDists[ i ].pdf( person.fakeY )

        for graph in self.graphs:
            graph.setEmissionDist( emissionPdf )

            # initialize the fake emissions
            graph.resampleEmissions()

    def useRealY( self ):
        def emissionPdf( person, i ):
            return self.emissionDists[ i ].pdf( person.y )

        for graph in self.graphs:
            graph.setEmissionDist( emissionPdf )

    def countTransitions( self ):

        N = len( self.transDists )
        obs = np.zeros( ( N, N, N ) )

        def accumulator( parents, child, X, i ):
            obs[ X[ 0 ] ][ X[ 1 ] ][ i ] += 1

        for graph in self.graphs:
            graph.countTransitions( accumulator )

        return obs

    def countEmissions( self, trueEmissions=True ):

        N = len( self.emissionDists )
        for emissionDist in self.emissionDists:
            M = emissionDist.N
            break

        obs = np.zeros( ( N, M ) )

        if( trueEmissions ):
            def accumulator( node, i ):
                obs[ i ][ node.y ] += 1
        else:
            def accumulator( node, i ):
                obs[ i ][ node.fakeY ] += 1

        for graph in self.graphs:
            graph.countEmissions( accumulator )

        return obs

    def resampleTransitionDistribution( self ):
        # Sample from P( A | X, Y, L, π )

        obs = self.countTransitions()

        for i, matDist in enumerate( self.transDists ):
            for j, transDist in enumerate( matDist ):
                self.transDists[ i ][ j ].resample( x=obs[ i ][ j ] )

        return obs

    def resampleEmissionDistribution( self, trueEmissions=True ):
        # Sample from P( L | X, Y, A, π )

        obs = self.countEmissions( trueEmissions )

        for i, emissionDist in enumerate( self.emissionDists ):
            self.emissionDists[ i ].resample( x=obs[ i ] )

        return obs

    def resampleRootDistributions( self ):
        # Sample from P( π | X, Y, A, L )

        allObs = []
        for i, graph in enumerate( self.graphs ):
            graphObs = []
            for j, root in enumerate( graph.roots ):
                obs = np.zeros( root.N )
                obs[ root.x ] = 1
                self.rootDists[ i ][ j ].resample( x=obs )
                graphObs.append( root.x )
            allObs.append( graphObs )

        return allObs


    def resampleGraphs( self ):
        # Sample from P( X | A, L, π, Y )

        nGraphs = len( self.graphs )
        for i, graph in enumerate( self.graphs ):
            graph.resampleGraphStates()

    def resample( self ):
        # Sample from P( X, A, L, π | Y )

        self.resampleGraphs()
        transCounts = self.resampleTransitionDistribution()
        emissionCounts = self.resampleEmissionDistribution()
        # because there is a different prior on each root, count is going to be 1
        rootStates = self.resampleRootDistributions()
        return ( transCounts, emissionCounts, rootStates )

    def log_joint( self ):
        # P( Y, X, A, L, π )

        val = 0

        for i, graph in enumerate( self.graphs ):

            # P( Y, X | A, L, π )
            val = graph.log_joint()

            # P( π )
            for j, root in enumerate( graph.roots ):
                val = self.rootDists[ i ][ j ].log_likelihood()

        # P( A )
        for i, matDist in enumerate( self.transDists ):
            for j, transDist in enumerate( matDist ):
                val = self.transDists[ i ][ j ].log_likelihood()

        # P( L )
        for i, emissionDist in enumerate( self.emissionDists ):
            val = self.emissionDists[ i ].log_likelihood()

        return val

    def modelEvidence( self, mixin=20, samples=50 ):
        # Use monte carlo integration to approximate P( Y )

        # Let the chain mix
        for i in range( mixin ):
            self.resample()
            joint = self.log_joint()
            # print( 'joint: %f'%joint )
            # sys.stdout.flush()

        # Accumulate the integral
        integral = RunningStats( useLogVar=True )
        for i in range( samples ):
            self.resample()
            joint = self.log_joint()
            # print( 'joint: %f'%joint )
            integral.pushVal( joint, isLog=True )

            if( i < 5 ):
                continue
            # print( integral.log_mean() )
            # print( integral.log_variance() )
            print('%d, %f, %f'%( i, integral.log_mean(), integral.log_variance() ))
            sys.stdout.flush()

        return ( integral.log_mean(), integral.log_variance() )


class WeakLimitHDPHHMM( HHMMModelBase ):

    # Put the same DP prior on all of the transition, emission and root vectors

    def __init__( self, graphs, transHDP, emissionHDP ):

        self.graphs = graphs

        self.transHDP = transHDP
        self.transDists = []

        self.emissionHDP = emissionHDP
        self.emissionDists = []

        N = transHDP.K
        self.transShape = ( N, N, N )

        M = emissionHDP.K
        self.emissionShape = ( N, M )

        # assume all roots have same hyper parameters
        # and make the weakest assumption about them
        self._rootHypers = np.ones( N )
        self.rootDists = []

        self.initializeGraphParameters()

        # preprocess graphs
        for graph in self.graphs:

            feedbackSet, blockManager = identifyCycles( graph )
            feedbackSetIds = [ node._id for node in feedbackSet ]
            graph.preprocess( feedbackSetIds )
            graph.draw()

    def _initializeTransDists( self ):
        self.transDists = []

        # make the transition parameters
        for i in range( self.transShape[ 0 ] ):
            transDist = []
            for j in range( self.transShape[ 1 ] ):
                transDist.append( Categorical( params=self.transHDP.sampleProbs(), alpha=None ) )
            self.transDists.append( transDist )

    def _initializeEmissionDists( self ):
        self.emissionDists = []

        # make the emission parameters
        for i in self.range( self.transShape[ 0 ] ):
            self.emissionDists.append( Categorical( params=self.emissionHDP.sampleProbs(), alpha=None ) )

    def _initializeRootDists( self ):

        self.rootDists = []

        for graph in self.graphs:

            # make the root distribution parameters
            rootDist = [ Categorical( alpha=self._rootHypers ) for root in sorted( graph.roots ) ]
            self.rootDists.append( rootDist )


    def initializeGraphParameters( self ):

        self._initializeTransDists()
        self._initializeEmissionDists()
        self._initializeRootDists()

        for graph in self.graphs:

            transitionPdf = self._transPdfClosure( graph )
            emissionPdf   = self._emissionPdfClosure( graph )
            rootPdf       = self._rootPdfClosure( graph )
            graph.setParameters( transitionPdf, emissionPdf, rootPdf )

    def resampleHDPs( self, transCounts, emissionCounts ):

        self.transHDP.resample( transMs )
        self.emissionHDP.resample( emissionMs )

        # Update the hyper parameters for each categorical distribution
        self._initializeTransDists()
        self._initializeEmissionDists()


    def resample( self ):
        # Sample from P( X, A, L, π, β, m | Y )

        transCounts, emissionCounts, rootStates = super( WeakLimitHDPHHMM, self ).resample()

        self.resampleTransitionTableCounts( transCounts )
        self.resampleEmissionTableCounts( emissionCounts )
        self.resampleHDPs()