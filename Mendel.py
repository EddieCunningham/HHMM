from HHMMMessagePasser import HiddenMarkovModelMessagePasser
from HHMMNode import NodeForHMM
from HHMMHG import MessagePassingHG
from Distributions import Categorical
import numpy as np
from cycleDetector import identifyCycles


class AutosomalMendelModel():

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

        self.resample()

    # Closure so that to the graph the root distributions
    # look the same way as the transition and emission
    # pdf functions
    def _transPdfClosure( self, graph ):

        def transPdf( parents, child, X, i ):
            assert parents[ 0 ].sex == 'female'
            return self.transDists[ X[ 0 ] ][ X[ 1 ] ].pdf( i )

        return transPdf

    def _emissionPdfClosure( self, graph ):
        def emissionPdf( person, i ):
            return self.emissionDists[ i ].pdf( person.y )

        return emissionPdf

    def _rootPdfClosure( self, graph ):
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

    def resampleTransitionDistribution( self ):
        # Sample from P( A | X, Y, L, π )

        N = len( self._transHypers )
        obs = np.zeros( ( N, N, N ) )

        def accumulator( parents, child, X, i ):
            obs[ X[ 0 ] ][ X[ 1 ] ][ i ] += 1

        for graph in self.graphs:
            graph.countTransitions( accumulator )

        for i, matHyper in enumerate( self._transHypers ):
            for j, transHyper in enumerate( matHyper ):
                self.transDists[ i ][ j ].resample( observations=obs[ i, j, : ] )

    def resampleEmissionDistribution( self, trueEmissions=True ):
        # Sample from P( L | X, Y, A, π )

        N = len( self._emissionHypers )
        for emissionHyper in self._emissionHypers:
            M = len( emissionHyper )
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

        for i,emissionHyper in enumerate( self._emissionHypers ):
            self.emissionDists[ i ].resample( observations=obs[ i, : ] )

    def resampleRootDistributions( self ):
        # Sample from P( π | X, Y, A, L )
        for i, graph in enumerate( self.graphs ):
            for j, root in enumerate( graph.roots ):
                obs = np.zeros( root.N )
                obs[ root.x ] = 1
                self.rootDists[ i ][ j ].resample( observations=obs )

    def resampleGraphs( self ):
        # Sample from P( X | θ, Y )

        nGraphs = len( self.graphs )
        for i, graph in enumerate( self.graphs ):
            graph.resampleGraphStates()

    def resample( self ):
        # Sample from P( X, θ | Y )

        self.resampleGraphs()
        self.resampleTransitionDistribution()
        self.resampleEmissionDistribution()
        self.resampleRootDistributions()

    def log_joint( self ):
        # P( Y, X, A, L, π )

        val = 0

        for graph in self.graphs:

            # P( Y, X | A, L, π )
            val += graph.log_joint()

            # P( π )
            for j, root in enumerate( graph.roots ):
                val += self.rootDists[ i ][ j ].log_likelihood()

        # P( A )
        for i, matHyper in enumerate( self._transHypers ):
            for j, transHyper in enumerate( matHyper ):
                val += self.transDists[ i ][ j ].log_likelihood()

        # P( L )
        for i, emissionHyper in enumerate( self._emissionHypers ):
            val += self.emissionDists[ i ].log_likelihood()

        return val

