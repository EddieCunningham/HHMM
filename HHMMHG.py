from HypergraphBase import BaseHyperGraph
from HHMMNode import NodeForHMM
import graphviz

class MessagePassingHG( BaseHyperGraph ):
    def __init__( self, N=None ):
        self.N = N
        super( MessagePassingHG, self ).__init__()
        self.setNodeType( NodeForHMM )

    def addNode( self, ID, y=0, N=None ):
        if( N == None ):
            N = self.N
        return super( MessagePassingHG, self ).addNode( ID, y, N )

    def draw( self ):

        assert self._initialized, 'call the function \'hypergraph.initialize()\''

        """ Draws the hypergraph using graphviz """
        d = graphviz.Digraph()
        for e in self._edges:
            eId = e._id
            for p in e._parents:
                pId = p._id
                d.edge( 'n( '+str( pId )+' )', 'E( '+str( eId )+' )' )
            for c in e._children:
                cId = c._id
                d.edge( 'E( '+str( eId )+' )', 'n( '+str( cId )+' )' )
        d.render()
        return d
