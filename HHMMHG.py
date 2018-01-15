from HypergraphBase import BaseHyperGraph
from HHMMNode import NodeForHMM
from HHMMMessagePasser import HiddenMarkovModelMessagePasser
import graphviz
import numpy as np
from pyLogVar import LogVar
from Distributions import Categorical

class MessagePassingHG( BaseHyperGraph ):
    def __init__( self, N=None, NObs=None, NodeType=None ):
        self.N = N
        self.NObs = NObs
        super( MessagePassingHG, self ).__init__()
        if( NodeType ):
            self.setNodeType( NodeType )
        else:
            self.setNodeType( NodeForHMM )
        self._msg = HiddenMarkovModelMessagePasser( self )

        # dummy variable so that we can use the graph
        # iterator easily
        self.workspace = LogVar( 0 )

    def setTransitionDist( self, transFunc ):
        self._msg._trans = transFunc

    def setEmissionDist( self, emissionFunc ):
        self._msg._L = emissionFunc

    def setRootDist( self, rootDist ):
        self._msg._pi = rootDist

    def setParameters( self, transFunc, emissionFunc, rootFunc ):
        self._msg.setParameters( transFunc, \
                                 emissionFunc, \
                                 rootFunc )

    def addNode( self, ID, y=0, N=None, NObs=None ):
        if( N == None ):
            N = self.N
        if( NObs == None ):
            NObs = self.NObs
        return super( MessagePassingHG, self ).addNode( ID, y, N, NObs )

    def draw( self, render=True ):

        d = super( MessagePassingHG, self ).draw( render=False )

        for node in self._msg._feedbackSet:
            d.node('n( '+str( node._id )+' )', **{
                'fontcolor': 'white',
                'style': 'filled',
                'fillcolor': 'green3',
            })

        for node in self._msg._sortaRootDeps:
            d.node('n( '+str( node._id )+' )', **{
                'fontcolor': 'white',
                'style': 'filled',
                'fillcolor': 'blue1',
            })

        if( render ):
            d.render()

        return d

    def genCode( self ):
        assert self._initialized, 'call the function \'hypergraph.initialize()\''

        code = 'hg = MessagePassingHG( 2 )\n'

        for node in self._nodes:
            code += 'n%d = hg.addNode( %d )\n'%( node._id, node._id )

        for edge in self._edges:
            code += 'e%d = hg.addEdge( set( ['%( edge._id )

            for node in edge._parents:
                code += ' n%d,'%( node._id )

            code = code[:-1]+' '
            code += ' ] ), %d )\n'%( edge._id )


            for node in edge._children:
                code += 'e%d.addChild( n%d )\n'%( edge._id, node._id )

        return code

    def preprocess( self, feedbackSetIds ):
        self._msg.preprocess( feedbackSetIds )

    def resample( self ):
        self.resampleGraphStates()

    def log_joint( self ):
        # P( Y, X | θ )
        self.workspace = LogVar( 1 )

        def nodeWork( node ):
            if( len( node._parents ) == 0 ):
                self.workspace *= LogVar( self._msg._pi( node, node.x ) )
            else:
                X = tuple( [ n.x for n in node._parents ] )
                self.workspace *= LogVar( self._msg._trans( node._parents, node, X, node.x ) )
            self.workspace *= LogVar( self._msg._L( node, node.x ) )

        self.graphIterate( nodeWork )

        return self.workspace.logVal

    def log_likelihood( self ):
        # P( Y | X, θ )
        self.workspace = LogVar( 1 )

        def nodeWork( node ):
            self.workspace *= LogVar( self._msg._L( node, node.x ) )

        self.graphIterate( nodeWork )

        return self.workspace.logVal

    def log_obsProb( self ):
        # P( Y | θ )
        return self._msg.probOfAllNodeObservations().logVal

    def log_posterior( self ):
        # P( X | Y, θ )
        return self.log_joint() - self.log_obsProb()

    def resampleGraphStates( self ):
        # Sample from P( X | Y, θ )

        # Run the message passing algorithm
        self._msg.getStats()

        def nodeWork( node ):

            if( len( node._parents ) == 0 ):
                probs = [ self._msg._pi( node, i ) for i in range( node.N ) ]
                node.x = Categorical.sample( probs, normalized=True )
            else:
                X = tuple( [ n.x for n in node._parents ] )
                probs = [ self._msg.conditionalParentChild( node, X, i ) for i in range( node.N ) ]
                node.x = Categorical.sample( probs, normalized=False )

        self.graphIterate( nodeWork )

    def resampleGraphStatesAndEmissions( self ):
        # Sample from P( X, Y | θ )
        def nodeWork( node ):

            if( len( node._parents ) == 0 ):
                stateProbs = [ self._msg._pi( node, i ) for i in range( node.N ) ]
                node.x = Categorical.sample( stateProbs, normalized=True )
            else:
                X = tuple( [ n.x for n in node._parents ] )
                stateProbs = [ self._msg._trans( node._parents, node, X, i ) for i in range( node.N ) ]
                node.x = Categorical.sample( stateProbs, normalized=True )

            emissionProbs = []
            for i in range( node.NObs ):
                node.fakeY = i
                emissionProbs.append( self._msg._L( node, node.x ) )
            node.fakeY = Categorical.sample( emissionProbs, normalized=True )

        self.graphIterate( nodeWork )

    def resampleEmissions( self ):
        # Sample from P( Y | X, θ )
        def nodeWork( node ):
            emissionProbs = []
            for i in range( node.NObs ):
                node.fakeY = i
                emissionProbs.append( self._msg._L( node, node.x ) )
            node.fakeY = Categorical.sample( emissionProbs, normalized=True )

        self.graphIterate( nodeWork )

    def countTransitions( self, accumulator ):
        # Count the # of times each transition is used.
        # Accumulator should be some function that defines
        # its own count for each parent child combo

        def nodeWork( node ):
            if( len( node._parents ) > 0 ):
                X = tuple( [ n.x for n in node._parents ] )
                accumulator( node._parents, node, X, node.x )

        self.graphIterate( nodeWork )

    def countEmissions( self, accumulator ):
        def nodeWork( node ):
            accumulator( node, node.x )
        self.graphIterate( nodeWork )
