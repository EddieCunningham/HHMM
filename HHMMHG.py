from HypergraphBase import BaseHyperGraph
from HHMMNode import NodeForHMM
from HHMMMessagePasser import HiddenMarkovModelMessagePasser
import graphviz
from HGTest import marginalizeTest
import numpy as np


class MessagePassingHG( BaseHyperGraph ):
    def __init__( self, N=None ):
        self.N = N
        super( MessagePassingHG, self ).__init__()
        self.setNodeType( NodeForHMM )
        self._msg = HiddenMarkovModelMessagePasser( self )

    def setParameters( self, transFunc, emissionFunc, rootFunc ):
        self._msg.setParameters( transFunc, \
                                 emissionFunc, \
                                 rootFunc )

    def addNode( self, ID, y=0, N=None ):
        if( N == None ):
            N = self.N
        return super( MessagePassingHG, self ).addNode( ID, y, N )

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

    def test( self, printStuff=False ):
        marginalizeTest( self._msg, printStuff )

    def resampleGraphStates( self ):
        current = self._roots

        for root in current:
            probs = [ self._msg._pi( root, i ) for i in range( root.N ) ]
            root.state = np.random.choice( root.N, 1, p=probs )[0]

        while( len( current ) > 0 ):

            nextCurrent = []

            for node in current:

                if( np.all( np.array( [ n.x for n in node._parents ] ) ) ):

                    probs = np.array( [ self._msg.conditionalParentChild( node, X, i ) for i in range( node.N ) ] )
                    total = LogVar( 0 )
                    for p in probs:
                        total += p

                    for i in range( node.N ):
                        probs[ i ] /= total

                    node.state = np.random.choice( node.N, 1, p=probs )[ 0 ]

                    for edge in node._downEdges:

                        for child in node._children:
                            child.x = None

                        nextCurrent.extend( node._children )
                else:
                    nextCurrent.append( node )

            current = list( set( nextCurrent ) )