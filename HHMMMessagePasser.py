from pyLogVar import LogVar
from cycleDetector import identifyCycles
from HHMMNode import NodeForHMM

import graphviz
import itertools
from util import prettyPrint

class HiddenMarkovModelMessagePasser():

    def __init__( self, hg ):

        self._hyperGraph = hg

        self.nodes = hg._nodes
        self.edges = hg._edges

        self._conditioning = {}

        self._srp = {}

    def setParameters( self, transFunc, emissionFunc, rootFunc ):

        self._trans = transFunc
        self._L = emissionFunc
        self._pi = rootFunc

    def preprocess( self, feedbackSetIds ):

        # for preprocessing, need to find a fvs, sortaRootProb

        self._feedbackSet = sorted( [ self._hyperGraph.getNode( n ) for n in feedbackSetIds ], key=lambda x: x._id )

        print( 'Feedback is: '+str( self._feedbackSet ) )

        for node in self.nodes: node.setMsg( self )

        # This is in case a person is basically a root because its parents are all in the fbs
        self._sortaRootDeps = [ node for node in self._feedbackSet if len( node._parents ) == 0 or \
                            ( len( [ x for x in node._parents if x not in self._feedbackSet ] ) == 0 and \
                             len( [ x for x in node._upEdge._children if x not in self._feedbackSet ] ) == 0 ) ]

        for node in self.nodes:
            if node in self._feedbackSet:
                node.inFBS = True

    """ -------------------------------------------------------------------------------------- """

    def sortaRootProb( self, conditioning ):

        # This accounts for individual disjoint nodes.
        # This same sort of logic can probably extend this
        # algorithm to handle disjoint graphs

        key = tuple( conditioning.items() )
        if( key in self._srp ):
            return self._srp[ key ]

        prod = LogVar( 1 )

        for sortaRoot in self._sortaRootDeps:
            j = conditioning[ sortaRoot ]
            _prod = LogVar( self._L( sortaRoot, j ) )

            if( len( sortaRoot._parents ) == 0 ):
                _prod *= self._pi( sortaRoot, j )
            else:
                parents = sortaRoot._parents
                X_ = [ conditioning[ p ] for p in parents ]
                _prod *= self._trans( parents, sortaRoot, X_, j )
            prod *= _prod

        self._srp[ key ] = prod
        return prod

    """ -------------------------------------------------------------------------------------- """

    def recursiveMessagePasser( self ):

        for node in self.nodes:
            if( node in self._feedbackSet ): continue

            for X in itertools.product( *[ range( n.N ) for n in self._feedbackSet ] ):
                conditioning = { node: i for node, i in zip( self._feedbackSet, X ) }
                for i in range( node.N ):
                    node.getU( i, conditioning )
                for i in range( node.N ):
                    for edge in node._downEdges:
                        node.getV( edge, i, conditioning )

        for node in self.nodes:
            if( node in self._feedbackSet ): continue

            node.accumulateFullJoint( self._feedbackSet )

        # compute the probs for the fbs
        aLeaf = self._hyperGraph._leaves.__iter__().__next__()
        self._srp = {}

        for X in itertools.product( *[ range( n.N ) for n in self._feedbackSet ] ):

            conditioning = { node: i for node, i in zip( self._feedbackSet, X ) }

            for i in range( aLeaf.N ):
                val = aLeaf.getU( i, conditioning )
                sr = self.sortaRootProb( conditioning )

                for node, x in zip( self._feedbackSet, X ):
                    node.updateFullJoint( x, val*sr )

        for node in self._feedbackSet:
            total = LogVar( 0 )
            for i in range( node.N ):
                total += node._fullJoint[ i ]

    """ -------------------------------------------------------------------------------------- """

    def isolatedParentJoint( self, node, X, i, totalProb=None ):

        # compute the probs for the fbs
        aLeaf = list( self._hyperGraph._leaves )[ 0 ]

        parents = node._parents
        jointProb = LogVar( 0 )

        parentCond = { p: j for p, j in zip( parents, X ) }

        margOut = [ n for n in self._feedbackSet if n!=node and n not in parents ]

        for _X in itertools.product( *[ range( n.N ) for n in margOut ] ):

            conditioning = { _node: _i for _node, _i in zip( margOut, _X ) }
            conditioning.update( parentCond )
            conditioning[ node ] = i

            for j in range( aLeaf.N ):
                val = aLeaf.getU( j, conditioning )
                sr = self.sortaRootProb( conditioning )
                jointProb += val*sr

        if( totalProb ):
            jointProb /= totalProb
        return jointProb

    def jointParentChild( self, node, X, i, totalProb=None ):

        # P( x_c, x_p, x_q | Y )

        if( node in self._sortaRootDeps ):
            return self.isolatedParentJoint( node, X, i, totalProb )

        parents = node._parents

        # make sure we sum over all of the possible nodes in the fbs
        latentRanges = [ [ X[ parents.index( n ) ] ] if n in parents else range( n.N ) for n in self._feedbackSet ]

        if( node.inFBS ):
            latentRanges[ self._feedbackSet.index( node ) ] = [ i ]

        """ Prob of this node """
        prob = LogVar( self._trans( parents, node, X, i ) ) * self._L( node, i )

        total = LogVar( 0 )

        # print( '\ni is: '+str( i ) )
        # print( 'node is: '+str( node ) )
        # print( 'X is: '+str( X ) )
        # print( 'parents is: '+str( parents ) )
        # print( 'self._feedbackSet is: '+str( self._feedbackSet ) )
        # print( 'latentRanges: '+str( latentRanges ) )
        for S in itertools.product( *latentRanges ):

            # print( 'S IS '+str( S ) )
            familyInFBS = { n: s for n, s in zip( self._feedbackSet, S ) }

            """ Down this node """
            _prob = LogVar( 1 )
            for e in node._downEdges:
                v = node.getMarginalizedV( e, i, familyInFBS, self._feedbackSet )
                # print( 'e: '+str( e )+' v: '+str( v ) )
                _prob *= v

            """ Out from each sibling """
            for sibling in node._upEdge._children:
                if( sibling is node ): continue
                b = sibling.getMarginalizedB( X, familyInFBS, self._feedbackSet )
                # print( 'sibling: '+str( sibling )+' b: '+str( b ) )
                _prob *= b

            """ Out from each parent """
            for parent, j in zip( parents, X ):
                a = parent.getMarginalizedA( node._upEdge, j, familyInFBS, self._feedbackSet )
                # print( 'parent: '+str( parent )+' edge: '+str( node._upEdge )+' j: '+str( j )+' a: '+str( a ) )
                _prob *= a

            srp = self.sortaRootProb( familyInFBS )
            # print( 'srp: '+str( srp ) )
            _prob *= srp
            total += _prob

        prob *= total

        """ Normalize """
        if( totalProb ):
            prob /= totalProb

        return prob

    def conditionalParentChild( self, node, X, i, totalProb=None ):

        # P( x_c | x_p, x_q, Y )

        if( node in self._sortaRootDeps ):
            # need to do something special in this case probably
            assert 0
            return self.isolatedParentJoint( node, X, i, totalProb )

        parents = node._parents

        # make sure we sum over all of the possible nodes in the fbs
        latentRanges = [ [ X[ parents.index( n ) ] ] if n in parents else range( n.N ) for n in self._feedbackSet ]

        if( node.inFBS ):
            latentRanges[ self._feedbackSet.index( node ) ] = [ i ]

        """ Prob of this node """
        prob = LogVar( self._trans( parents, node, X, i ) ) * self._L( node, i )

        total = LogVar( 0 )

        # print( '\ni is: '+str( i ) )
        # print( 'node is: '+str( node ) )
        # print( 'X is: '+str( X ) )
        # print( 'parents is: '+str( parents ) )
        # print( 'self._feedbackSet is: '+str( self._feedbackSet ) )
        # print( 'latentRanges: '+str( latentRanges ) )
        for S in itertools.product( *latentRanges ):

            # print( 'S IS '+str( S ) )
            familyInFBS = { n: s for n, s in zip( self._feedbackSet, S ) }

            """ Down this node """
            _prob = LogVar( 1 )
            for e in node._downEdges:
                v = node.getMarginalizedV( e, i, familyInFBS, self._feedbackSet )
                # print( 'e: '+str( e )+' v: '+str( v ) )
                _prob *= v

            _prob /= node.getMarginalizedB( X, self._feedbackSet )

            srp = self.sortaRootProb( familyInFBS )
            # print( 'srp: '+str( srp ) )
            _prob *= srp
            total += _prob

        prob *= total

        """ Normalize """
        if( totalProb ):
            prob /= totalProb

        return prob

    def messagePasser( self ):

        for X in itertools.product( *[ range( n.N ) for n in self._feedbackSet ] ):
            conditioning = { node: i for node, i in zip( self._feedbackSet, X ) }

            restart = True
            while( restart ):
                restart = False

                lastVList = []
                lastUList = []

                newRoots = []
                for fbNode in self._feedbackSet:
                    for edge in fbNode._downEdges:
                        newRoots.extend( edge._children )

                newLeaves = []
                for fbNode in self._feedbackSet:
                    newLeaves.extend( fbNode._parents )

                # start at leaves / roots and work in
                currentVList = set( list( self._hyperGraph._leaves ) + newLeaves ) - set( self._feedbackSet )
                currentUList = set( list( self._hyperGraph._roots ) + newRoots ) - set( self._feedbackSet )

                # print( '\n============\n============\n============\n' )
                # print( 'STARTING WITH GOING UP: %s AND GOING DOWN: %s'%( str( currentVList ), str( currentUList ) ) )

                while( len( currentVList ) + len( currentUList ) > 0 ):

                    # print( '\n=========\n=========\n' )
                    # print( '\nGOING UP: %s AND GOING DOWN: %s'%( str( currentVList ), str( currentUList ) ) )
                    nextVList = set()
                    nextUList = set()

                    # Nodes going down
                    for node in currentUList:
                        if( node.inFBS ): continue

                        addChildren = False
                        for i in range( node.N ):
                            if( node.UReady( i, conditioning ) ):
                                node.getU( i, conditioning, 0 )
                                addChildren = True

                        if( addChildren ):
                            # print( 'Computed U for node %s'%( node ) )

                            for edge in node._downEdges:

                                for mate in edge._parents:
                                    if( mate.inFBS or mate == node ):
                                        continue
                                    # print( 'ADDING %s TO GOING UP'%mate )
                                    nextVList.add( mate )

                                for child in edge._children:
                                    if( child.inFBS ):
                                        continue
                                    # print( 'ADDING %s TO GOING DOWN'%child )
                                    nextUList.add( child )

                        else:
                            # print( 'READDING %s TO GOING DOWN'%node )
                            nextUList.add( node )

                    #####################################################################################################

                    # Nodes going up
                    for node in currentVList:
                        if( node.inFBS ): continue

                        addParents = True
                        if( len( node._downEdges ) > 0 ):
                            for edge in node._downEdges:

                                for i in range( node.N ):

                                    if( node.VReady( edge, i, conditioning ) ):
                                        node.getV( edge, i, conditioning, 0 )
                                    else:
                                        addParents = False

                        if( addParents ):
                            # print( 'Computed V for node %s'%( node ) )

                            for parent in node._parents:
                                if( parent.inFBS ):
                                    continue
                                # print( 'ADDING %s TO GOING UP'%parent )
                                nextVList.add( parent )

                            if( node._upEdge ):
                                for sibling in node._upEdge._children:
                                    if( sibling.inFBS or sibling == node ):
                                        continue
                                    # print( 'ADDING %s TO GOING DOWN'%sibling )
                                    nextUList.add( sibling )
                        else:
                            # print( 'READDING %s TO GOING UP'%node )
                            nextVList.add( node )

                    if( lastVList == currentVList and lastUList == currentUList ):
                        # This hack is for the event where we don't gather
                        # all of the node dependencies the first time.  If we don't
                        # do that, then the keys get messed up and we can't find the
                        # a,b,U,V values that we need.  If we hit this, then we
                        # have most likely gathered the relevant dependencies
                        # and the algorithm should work the next time around.
                        print( '\n\n\nFAILED WITH GOING UP: %s AND GOING DOWN: %s\n'%( str( currentVList ), str( currentUList ) ) )
                        restart = True
                        break

                    lastVList = currentVList
                    lastUList = currentUList

                    currentVList = nextVList
                    currentUList = nextUList


        for node in self.nodes:
            if( node in self._feedbackSet ): continue

            node.accumulateFullJoint( self._feedbackSet )

        # compute the probs for the fbs
        aLeaf = self._hyperGraph._leaves.__iter__().__next__()
        self._srp = {}

        for X in itertools.product( *[ range( n.N ) for n in self._feedbackSet ] ):

            conditioning = { node: i for node, i in zip( self._feedbackSet, X ) }

            for i in range( aLeaf.N ):
                val = aLeaf.getU( i, conditioning )
                sr = self.sortaRootProb( conditioning )

                for node, x in zip( self._feedbackSet, X ):
                    node.updateFullJoint( x, val*sr )

        for node in self._feedbackSet:
            total = LogVar( 0 )
            for i in range( node.N ):
                total += node._fullJoint[ i ]

    def getStats( self ):

        for node in self.nodes:
            node.reset()

        self.messagePasser()

        for node in self.nodes:
            for i in range( node.N ):
                uVal = node.getFullJoint( i )

        for n in self.nodes:
            for i in range( n.N ):
                if( len( n._parents ) > 0 ):
                    for X in itertools.product( *[ range( p.N ) for p in n._parents ] ):
                        genProb = self.jointParentChild( n, X, i )

    """ -------------------------------------------------------------------------------------- """

    def probOfAllNodeObservations( self ):

        """ P( Y ) = sum_i( U_l_e( i ) ) for any leaf l """
        aLeaf = list( self._hyperGraph._leaves )[ 0 ]

        total = LogVar( 0 )
        for i in range( aLeaf.N ):
            _u = aLeaf.getFullJoint( i )
            total += _u
        return total
