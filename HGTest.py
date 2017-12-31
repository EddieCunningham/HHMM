from HHMMMessagePasser import HiddenMarkovModelMessagePasser
from pyLogVar import LogVar
from HHMMNode import NodeForHMM
from util import prettyPrint
import itertools


def aTest( messagePasser, printStuff=False ):

    messagePasser.getStats()

    if( printStuff ):
        for node in messagePasser.nodes:
            print( '\n----------------\nNode: '+str( node ) )
            print( '_a:\n' )
            prettyPrint( node._a )
            print( '_b:\n' )
            prettyPrint( node._b )
            print( '_U:\n' )
            prettyPrint( node._U )
            print( '_V:\n' )
            prettyPrint( node._V )
            print( 'U deps: '+str( node._UDeps ) )
            print( 'V deps: '+str( node._VDeps ) )
            print( 'a deps: '+str( node._aDeps ) )
            print( 'b deps: '+str( node._bDeps ) )

    correct = messagePasser.probOfAllNodeObservations()
    if( printStuff ):
        print( '\n\nP( Y ) is: '+str( correct ) )
        print( 'Feedback set: '+str( messagePasser._feedbackSet ) )

    for node in messagePasser.nodes:
        total = LogVar( 0 )
        for i in range( node.N ):

            uVal = node.getFullJoint( i )
            uCopy = LogVar( uVal )
            uCopy /= correct
            total += uCopy

    def isclose( a, b, rel_tol=1e-05, abs_tol=0.0 ):
        return abs( a-b ) <= max( rel_tol * max( abs( a ), abs( b ) ), abs_tol )

    """ Make sure that P( Y ) can be computed everywhere correctly """
    for node in messagePasser.nodes:
        total = LogVar( 0 )
        for i in range( node.N ):

            uVal = node.getFullJoint( i )
            uCopy = LogVar( uVal )
            total += uCopy

        if( not isclose( total.logVal, correct.logVal ) ):
            if( printStuff ):
                print( '\nFailed the P( Y ) test!' )
                print( 'Node: '+str( node ) )
                print( 'total: '+str( float( total ) ) )
                print( 'correct: '+str( float( correct ) ) )
            print( 'up edge: '+str( node._upEdge ) )
            print( 'down edges: '+str( node._downEdges ) )
            assert 0
        else:
            if( printStuff ):
                print( 'Passed P( Y ) test for node '+str( node )+' with prob '+str( total ) )
    if( printStuff ):
        print( 'Passed the P( Y ) tests!' )

    """ Make sure that the statistics we care about sum to 1 """
    for n in messagePasser.nodes:

        total1 = LogVar( 0 )
        total2 = LogVar( 0 )
        for i in range( n.N ):

            prob = n.getFullJoint( i )
            prob /= correct
            total1 += prob

            if( len( n._parents ) > 0 ):

                shouldEqualProb = LogVar( 0 )
                for X in itertools.product( *[ range( p.N ) for p in n._parents ] ):

                    prob2 = messagePasser.probOfParentsProducingNode( n, X, i, correct )
                    total2 += prob2
                    shouldEqualProb += prob2

                if( not isclose( float( shouldEqualProb ), float( prob ) ) ):
                    if( printStuff ):
                        print( '\n\n\n\nFailed this test for node: '+str( n ) )
                        print( 'This is wrong: '+str( float( shouldEqualProb*correct ) ) )
                        print( 'This is right: '+str( float( prob*correct ) ) )
                    print( 'up edge: '+str( n._upEdge ) )
                    print( 'down edges: '+str( n._downEdges ) )
                    assert 0
                else:
                    if( printStuff ):
                        print( 'Passed the marg out test P( x_c|y )=sumP( x_c, {x_p}|y ) for node '+str( n )+' with prob '+str( prob ) )
                        print( '==========================================================' )

        if( not isclose( float( total1 ), 1. ) ):
            if( printStuff ):
                print( 'The sum of all P( x_i|Y ) doesn\'t equal 1 for node: '+str( n )+', total1: '+str( total1.logVal ) )
            print( 'up edge: '+str( node._upEdge ) )
            print( 'down edges: '+str( node._downEdges ) )
            assert 0

        if( len( n._parents ) > 0 ):
            if( not isclose( float( total2 ), 1. ) ):
                if( printStuff ):
                    'The sum of all P( {x_p}|Y ) isn\'t 1 for node: '+str( n )+', total2: '+str( total2.logVal )
                print( 'up edge: '+str( node._upEdge ) )
                print( 'down edges: '+str( node._downEdges ) )
                assert 0
