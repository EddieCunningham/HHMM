from HypergraphBase import NodeBase, EdgeBase
from pyLogVar import LogVar
import itertools
from util import prettyPrint


class NodeForHMM( NodeBase ):
    def __init__( self, y, N, NObs ):
        super( NodeForHMM, self ).__init__()

        self.N = N
        self.NObs = NObs
        self.y = y

        self.fakeY = None

        self.x = None

        self.inFBS = False

        self._UDeps = set( [ self ] )
        self._VDeps = set( [ self ] )
        self._aDeps = set( [ self ] )
        self._bDeps = set( [ self ] )

        self.reset()

        self.doneUpEdge = {}
        self.doneDownEdges = {}

    def reset( self ):
        self._U = {}
        self._V = {}
        self._a = {}
        self._b = {}
        self._fullJoint = {}

    def setMsg( self, msg ):
        self._msg = msg
        self._pi = msg._pi
        self._L = msg._L
        self._trans = msg._trans

    def _getN( self, node, conditioning ):
        if( node in conditioning ):
            return ( conditioning[ node ], )
        return range( node.N )

    def keyFromCond( self, conditioning ):
        return tuple( conditioning.items() )

    """ ------------------------------------------------------------------------------------ """

    def aKey( self, conditioning ):
        actualCond = { k:v for k, v in conditioning.items() if k in self._aDeps }
        key = self.keyFromCond( actualCond )
        return key

    def needToComputeA( self, edge, i, key ):
        return edge not in self._a or i not in self._a[ edge ] or key not in self._a[ edge ][ i ]

    def setAVal( self, edge, i, key, aVal ):
        if( edge not in self._a ): self._a[ edge ] = {}
        if( i not in self._a[ edge ] ): self._a[ edge ][ i ] = {}
        assert key not in self._a[ edge ][ i ]
        self._a[ edge ][ i ][ key ] = aVal

    def getAVal( self, edge, i, key ):
        try:
            return LogVar( self._a[ edge ][ i ][ key ] )
        except:
            print( 'YOU DONE FUCKED UP' )
            prettyPrint( self._a )
            self._a[ edge ]
            self._a[ edge ][ i ]
            self._a[ edge ][ i ][ key ]
            assert 0

    def getA( self, edge, i, conditioning, depth=0 ):

        key = self.aKey( conditioning )

        """ a_n_e( i ) = P( Y \ !( e, n ), n_x = i ) """
        if( self.needToComputeA( edge, i, key ) ):

            # if we don't have the right conditioning, add it
            aVal = self._computeA( edge, i, conditioning, depth )

            # revise the key
            key = self.aKey( conditioning )

            self.setAVal( edge, i, key, aVal )

        else:
            aVal = self.getAVal( edge, i, key )
        return aVal

    def _computeA( self, edge, i, conditioning, depth=0 ):

        a_ = LogVar( 1 )

        """ All nodes up from this node """
        a_ *= self.getU( i, conditioning, depth )
        self._aDeps |= self._UDeps

        """ All nodes down from this node but not down edge """
        if( len( self._downEdges ) > 1 ):

            for e in self._downEdges:
                if( e == edge ): continue

                """ Down each mate branch """
                a_ *= self.getV( e, i, conditioning, depth )
            self._aDeps |= self._VDeps

        return a_

    def getMarginalizedA( self, edge, i, nodesToKeep, feedbackSet ):

        margOut = filter( lambda n:n not in nodesToKeep and n in self._aDeps, feedbackSet )
        a = LogVar( 0 )
        for X in itertools.product( *[ range( n.N ) for n in margOut ] ):
            conditioning = { node:_i for node, _i in zip( margOut, X ) }
            conditioning.update( nodesToKeep )
            a += self.getA( edge, i, conditioning )
        return a

    def aReady( self, edge, i, conditioning ):

        if( self.inFBS ):
            return True

        # only ready if U and V are ready
        uKey = self.UKey( conditioning )
        vKey = self.VKey( conditioning )

        if( self.needToComputeU( i, uKey ) ):
            return False

        for e in self._downEdges:
            if( e == edge ): continue

            if( self.needToComputeV( e, i, vKey ) ):
                return False

        return True

    """ ------------------------------------------------------------------------------------ """

    def bKey( self, conditioning ):
        actualCond = { k:v for k, v in conditioning.items() if k in self._bDeps }
        key = self.keyFromCond( actualCond )
        return key

    def needToComputeB( self, X, key ):
        return X not in self._b or key not in self._b[ X ]

    def setBVal( self, X, key, bVal ):
        if( X not in self._b ): self._b[ X ] = {}
        self._b[ X ][ key ] = bVal

    def getBVal( self, X, key ):
        try:
            return LogVar( self._b[ X ][ key ] )
        except:
            print( 'YOU DONE FUCKED UP' )
            prettyPrint( self._b )
            self._b[ X ][ key ]
            assert 0

    def getB( self, X, conditioning, depth=0 ):

        key = self.bKey( conditioning )

        """ B_n( X ) = P( n_y, Y \ ^( n )_y | ^( n )_x = X ) """
        if( self.needToComputeB( X, key ) ):

            bVal = self._computeB( X, conditioning, depth )

            # revise the key
            key = self.bKey( conditioning )

            self.setBVal( X, key, bVal )

        else:
            bVal = self.getBVal( X, key )
        return bVal

    def _computeB( self, X, conditioning, depth=0 ):

        b_ = LogVar( 0 )

        for k in self._getN( self, conditioning ):

            """ Prob of this sibling """
            _prod = LogVar( self._trans( self._parents, self, X, k ) ) * self._L( self, k )

            """ Branch down from sibling """
            for e in self._downEdges:
                _prod *= self.getV( e, k, conditioning, depth )

            b_ += _prod

        self._bDeps |= self._VDeps

        return b_

    def getMarginalizedB( self, X, nodesToKeep, feedbackSet ):

        margOut = filter( lambda n:n not in nodesToKeep and n in self._bDeps, feedbackSet )
        b = LogVar( 0 )
        for _X in itertools.product( *[ range( n.N ) for n in margOut ] ):
            conditioning = { node:i for node, i in zip( margOut, _X ) }
            conditioning.update( nodesToKeep )
            b += self.getB( X, conditioning )
        return b

    def bReady( self, X, conditioning ):

        if( self.inFBS or len( self._downEdges ) == 0 ):
            return True

        # only ready if V is ready
        vKey = self.VKey( conditioning )

        for k in self._getN( self, conditioning ):
            for e in self._downEdges:

                if( self.needToComputeV( e, k, vKey ) ):
                    return False

        return True

    """ ------------------------------------------------------------------------------------ """

    def UKey( self, conditioning ):
        actualCond = { k:v for k, v in conditioning.items() if k in self._UDeps }
        key = self.keyFromCond( actualCond )
        return key

    def needToComputeU( self, i, key ):
        return i not in self._U or key not in self._U[ i ]

    def setUVal( self, i, key, uVal ):
        if( i not in self._U ): self._U[ i ] = {}
        self._U[ i ][ key ] = uVal

    def getUVal( self, i, key ):
        return LogVar( self._U[ i ][ key ] )

    def getU( self, i, conditioning, depth=0 ):

        key = self.UKey( conditioning )

        """ U_n( i ) = P( n_y, ^( n )_y, n_x = i ) """
        if( self.needToComputeU( i, key ) ):

            uVal = self._computeU( i, conditioning, depth )

            # revise the key
            key = self.UKey( conditioning )

            self.setUVal( i, key, uVal )

        else:
            uVal = self.getUVal( i, key )
        return uVal

    def _computeU( self, i, conditioning, depth=0 ):

        parents = self._parents
        if( len( parents ) == 0 or self.inFBS ):

            if( self.inFBS ):
                # Prob if we conditioned on this
                u = LogVar( int( i==conditioning[ self ] ) )
            else:
                # Root dist prob
                u = LogVar( self._pi( self, i ) ) * self._L( self, i )
        else:

            u = LogVar( 0 )
            for X in itertools.product( *[ self._getN( p, conditioning ) for p in parents ] ):

                """ Prob of this node """
                prod = LogVar( self._trans( parents, self, X, i ) )

                """ Branch out from each parent """
                for parent, j in zip( parents, X ):
                    prod *= parent.getA( self._upEdge, j, conditioning, depth+1 )
                    self._UDeps |= parent._aDeps

                """ Branch out from each sibling """
                for sibling in self._upEdge._children:
                    if( sibling==self ):continue

                    prod *= sibling.getB( X, conditioning, depth+1 )
                    self._UDeps |= sibling._bDeps

                u += prod

            u *= self._L( self, i )

        return u

    def getMarginalizedU( self, i, nodesToKeep, feedbackSet ):

        margOut = filter( lambda n:n not in nodesToKeep and n in self._UDeps, feedbackSet )
        u = LogVar( 0 )
        for X in itertools.product( *[ range( n.N ) for n in margOut ] ):
            conditioning = { node:i for node, i in zip( margOut, X ) }
            conditioning.update( nodesToKeep )
            u += self.getU( i, conditioning )
        return u

    def UReady( self, i, conditioning ):

        parents = self._parents
        if( len( parents ) == 0 or self.inFBS ):
            return True
        else:
            for X in itertools.product( *[ self._getN( p, conditioning ) for p in parents ] ):

                for parent, j in zip( parents, X ):

                    if( parent.aReady( self._upEdge, j, conditioning ) == False ):
                        return False

                for sibling in self._upEdge._children:
                    if( sibling == self ):continue

                    if( sibling.bReady( X, conditioning ) == False ):
                        return False
        return True

    """ ------------------------------------------------------------------------------------ """

    def VKey( self, conditioning ):
        actualCond = { k:v for k, v in conditioning.items() if k in self._VDeps }
        key = self.keyFromCond( actualCond )
        return key

    def needToComputeV( self, edge, i, key ):
        return edge not in self._V or i not in self._V[ edge ] or key not in self._V[ edge ][ i ]

    def setVVal( self, edge, i, key, vVal ):
        if( edge not in self._V ): self._V[ edge ] = {}
        if( i not in self._V[ edge ] ): self._V[ edge ][ i ] = {}
        assert key not in self._V[ edge ][ i ]
        self._V[ edge ][ i ][ key ] = vVal

    def getVVal( self, edge, i, key ):
        return LogVar( self._V[ edge ][ i ][ key ] )

    def getV( self, edge, i, conditioning, depth=0 ):

        key = self.VKey( conditioning )

        """ V_n_e( i ) = P( !( e, n ) | n_x = i ) """
        if( self.needToComputeV( edge, i, key ) ):

            # if we don't have the right conditioning, add it
            vVal = self._computeV( edge, i, conditioning, depth )

            # revise the key
            key = self.VKey( conditioning )

            self.setVVal( edge, i, key, vVal )

        else:
            vVal = self.getVVal( edge, i, key )
        return vVal

    def _computeV( self, edge, i, conditioning, depth=0 ):

        if( self.inFBS or len( self._downEdges ) == 0 ):
            return LogVar( 1 )

        mates = [ x for x in edge._parents if x != self ]
        selfIndex = edge._parents.index( self )

        v = LogVar( 0 )
        for X_ in itertools.product( *[ self._getN( m, conditioning ) for m in mates ] ):

            prod = LogVar( 1 )

            """ Branch out from each mate """
            for mate, j in zip( mates, X_ ):

                prod *= mate.getA( edge, j, conditioning, depth+1 )
                self._VDeps |= mate._aDeps

            X = tuple( list( X_[ :selfIndex ] ) + [ i ] + list( X_[ selfIndex: ] ) )

            """ Branch out from each child """
            for child in edge._children:
                prod *= child.getB( X, conditioning, depth+1 )
                self._VDeps |= child._bDeps

            v += prod

        return v

    def getMarginalizedV( self, edge, i, nodesToKeep, feedbackSet ):

        margOut = filter( lambda n:n not in nodesToKeep and n in self._VDeps, feedbackSet )
        v = LogVar( 0 )
        for X in itertools.product( *[ range( n.N ) for n in margOut ] ):
            conditioning = { node:i for node, i in zip( margOut, X ) }
            conditioning.update( nodesToKeep )
            v += self.getV( edge, i, conditioning )
        return v

    def VReady( self, edge, i, conditioning ):

        if( self.inFBS or len( self._downEdges ) == 0 ):
            return True

        mates = [ x for x in edge._parents if x != self ]
        selfIndex = edge._parents.index( self )

        v = LogVar( 0 )
        for X_ in itertools.product( *[ self._getN( m, conditioning ) for m in mates ] ):

            for mate, j in zip( mates, X_ ):

                if( mate.aReady( edge, j, conditioning ) == False ):
                    return False

            X = tuple( list( X_[ :selfIndex ] ) + [ i ] + list( X_[ selfIndex: ] ) )

            for child in edge._children:
                if( child.bReady( X, conditioning ) == False ):
                    return False

        return True

    """ ------------------------------------------------------------------------------------ """

    def sortaRootProb( self, conditioning ):
        return self._msg.sortaRootProb( conditioning )

    def accumulateFullJoint( self, feedbackSet ):

        if( len( self._downEdges ) > 0 ):

            for i in range( self.N ):

                total = LogVar( 0 )
                for X in itertools.product( *[ range( n.N ) for n in feedbackSet ] ):
                    conditioning = { node:_i for node, _i in zip( feedbackSet, X ) }

                    prod = self.getU( i, conditioning )

                    for edge in self._downEdges:
                        prod *= self.getV( edge, i, conditioning )

                    prod *= self.sortaRootProb( conditioning )

                    total += prod

                self._fullJoint[ i ] = total
        else:
            for i in range( self.N ):

                total = LogVar( 0 )
                for X in itertools.product( *[ range( n.N ) for n in feedbackSet ] ):
                    conditioning = { node:_i for node, _i in zip( feedbackSet, X ) }

                    total += self.getU( i, conditioning ) * self.sortaRootProb( conditioning )

                self._fullJoint[ i ] = total

    def getFullJoint( self, i ):
        return self._fullJoint[ i ]

    def updateFullJoint( self, i, val ):
        if( i not in self._fullJoint ): self._fullJoint[ i ] = LogVar( 0 )
        self._fullJoint[ i ] += val
