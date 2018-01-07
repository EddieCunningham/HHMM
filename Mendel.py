from HHMMMessagePasser import HiddenMarkovModelMessagePasser
from HHMMNode import NodeForHMM
from HHMMHG import MessagePassingHG
from Distributions import Categorical



# This model assumes that the transition and emission distributions
# are known.  The things that are resampled are the root distributions
# and the latent states.
class PerfectMendelModel():

    def __init__( self, graphs, parameters, rootHypers ):

        self.graphs = graphs

        # assume all roots have same hyper parameters
        self._rootHypers = rootHypers
        self.rootDists = []

        if( not self.checkParameters( parameters ) ):
            assert 0

        self.initializeGraphParameters( **parameters )

        self.resample()

    # Closure so that to the graph the root distributions
    # look the same way as the transition and emission
    # pdf functions
    def _rootPdfClosure( self, graph ):
        roots = sorted( graph.roots )
        graphIndex = self.graphs.index( graph )

        def rootPdf( person, i ):
            personIndex = roots.index( person )
            return self.rootDists[ graphIndex ][ personIndex ].pdf( i )

        return rootPdf

    def initializeGraphParameters( self, A, L ):

        def transitionPdf( parents, child, X, i ):
            return A[ X[ 0 ] ][ X[ 1 ] ][ i ]

        def emissionPdf( person, i ):
            return L[ i ][ person.y ]

        for k, graph in enumerate( self.graphs ):

            rootDist = []
            for root in sorted( graph.roots ):
                rootDist.append( Categorical( alpha=self._rootHypers ) )
            self.rootDists.append( rootDist )

            rootPdfs = self._rootPdfClosure( graph )

            graph.setParameters( transitionPdf, emissionPdf, rootPdfs )

    def resampleRootDistributions( self ):
        # Sample from P( π | X, Y )
        for i, graph in enumerate( self.graphs ):
            for j, root in enumerate( graph.roots ):
                obs = np.zeros( root.N )
                obs[ root.x ] = 1
                self.rootDists[ i ][ j ].resample( observations=obs )

    def resampleGraphs( self ):
        # Sample from P( X | π, Y )
        for graph in self.graphs:
            graph.resampleGraphStates()

    def resample( self ):
        # Sample from P( X, π | Y )
        self.resampleRootDistributions()
        self.resampleGraphs()

    def ratioTest( self ):
        # P( X, π | Y ) / P( X', π' | Y ) =
        # P( X, π, Y ) / P( X', π', Y )
        self.resample()

        log_cond1 = 0
        log_joint1 = 0
        for i, graph in enumerate( self.graphs ):
            log_cond1  += graph.log_posterior()
            log_joint1 += graph.log_joint()

            for j, root in enumerate( graph.roots ):
                rootLogLikelihood = self.rootDists[ i ][ j ].log_likelihood()
                log_cond1  += rootLogLikelihood
                log_joint1 += rootLogLikelihood

        self.resample()

        log_cond2 = 0
        log_joint2 = 0
        for i, graph in enumerate( self.graphs ):
            log_cond2  += graph.log_posterior()
            log_joint2 += graph.log_joint()

            for j, root in enumerate( graph.roots ):
                rootLogLikelihood = self.rootDists[ i ][ j ].log_likelihood()
                log_cond2  += rootLogLikelihood
                log_joint2 += rootLogLikelihood

        return log_cond1 - log_cond2 == log_joint1 - log_joint2

    def _forwardSample( self ):
        # Sample from P( π )
        for i, graph in enumerate( self.graphs ):
            for j, root in enumerate( graph.roots ):
                self.rootDist[ i ][ j ].resample()

        # Sample from P( X, Y | π )
        for graph in self.graphs:
            graph.resampleGraphStatesAndEmissions()

    def _fullGibbsSample( self ):
        # Sample from P( X, π | Y )
        self.resample()

        # Sample from P( Y | X, π )
        for graph in self.graphs:
            graph.resampleEmissions()

    def _collectStats( self ):
        assert 0

    def _compareStats( self, stats1, stats2 ):
        assert 0

    def gewekeTest( self ):
        # Sample from P( X, π, Y ) using two methods:
        #  1. P( π ), then P( X, Y | π )
        #  2. P( X, π | Y ), then P( Y | X, π )
        # And make sure that they agree

        # use each person's fakeY so we don't override y
        def emissionPdf( person, i ):
            return L[ i ][ person.fakeY ]

        for graph in self.graphs:
            graph.setEmissionDist( emissionPdf )

        nIters = 1000

        stats1 = []
        for i in range( nIters ):
            self._forwardSample()
            stats1.append( self._collectStats() )

        stats2 = []
        for i in range( nIters ):
            self._fullGibbsSample()
            stats2.append( self._collectStats() )

        self._compareStats( stats1, stats2 )