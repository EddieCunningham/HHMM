from HHMMMessagePasser import HiddenMarkovModelMessagePasser
from pyLogVar import LogVar
from HHMMNode import NodeForHMM
from util import prettyPrint
import itertools

def isclose( a, b, rel_tol=1e-15, abs_tol=0.0 ):
    if( a == 'ZERO' or b == 'ZERO' ):
        if( a == 'ZERO' and b == 'ZERO' ):
            return True
        return False
    return abs( a-b ) < 1e-8

def marginalizeTest( graph, printStuff=False ):

    messagePasser = graph._msg

    messagePasser.getStats()
    allFailures = []

    # Calculates P( Y ) from a leaf
    observationProb = messagePasser.probOfAllNodeObservations()

    # Calculate sum_{ i } P( x = i, Y ) at each node and compare to P( Y )
    for node in messagePasser.nodes:
        total = LogVar( 0 )
        for i in range( node.N ):

            uVal = node.getFullJoint( i )
            uCopy = LogVar( uVal )
            total += uCopy

        if( not isclose( total.logVal, observationProb.logVal ) ):
            failString = '\n\nFailed the P( Y ) test!' + '\n'
            failString += 'Node: '+str( node ) + '\n'
            failString += 'total: '+str( total.logVal ) + '\n'
            failString += 'observationProb: '+str( observationProb.logVal ) + '\n'
            failString += 'up edge: '+str( node._upEdge ) + '\n'
            failString += 'down edges: '+str( node._downEdges ) + '\n'
            allFailures.append( failString )

    # Calculate sum_{ X } P( child_x = i, parents_x = X, Y )
    # and compare to P( child_x = i, Y )
    for n in messagePasser.nodes:

        for i in range( n.N ):

            prob = n.getFullJoint( i )

            if( len( n._parents ) > 0 ):

                shouldEqualProb = LogVar( 0 )
                for X in itertools.product( *[ range( p.N ) for p in n._parents ] ):
                    shouldEqualProb += messagePasser.jointParentChild( n, X, i )

                if( not isclose( shouldEqualProb.logVal, prob.logVal ) ):
                    failString = '\n\nFailed the joint parent child marginalization test for node: '+str( n ) + '\n'
                    failString += 'This is wrong: '+str( shouldEqualProb.logVal ) + '\n'
                    failString += 'This is right: '+str( prob.logVal ) + '\n'
                    failString += 'up edge: '+str( n._upEdge ) + '\n'
                    failString += 'down edges: '+str( n._downEdges ) + '\n'
                    allFailures.append( failString )

    # Calculate sum_{ i } P( child_x = i, parents_x = X, Y )
    # and compare to P( parents_x = X, Y )
    for n in messagePasser.nodes:

        if( len( n._parents ) > 0 ):

            for X in itertools.product( *[ range( p.N ) for p in n._parents ] ):

                prob = messagePasser.jointParents( n, X )

                total = LogVar( 0 )
                for i in range( n.N ):
                    total += messagePasser.jointParentChild( n, X, i )

                if( not isclose( total.logVal, prob.logVal ) ):
                    failString = '\n\nFailed the joint parent marginalization test for node: '+str( n ) + '\n'
                    failString += 'This is wrong: '+str( total.logVal ) + '\n'
                    failString += 'This is right: '+str( prob.logVal ) + '\n'
                    failString += 'up edge: '+str( n._upEdge ) + '\n'
                    failString += 'down edges: '+str( n._downEdges ) + '\n'
                    allFailures.append( failString )

    # Calculate P( child_x = i | parents_x = X, Y ) and see if we
    # summing over i gets us 1
    for n in messagePasser.nodes:

        if( len( n._parents ) > 0 ):

            for X in itertools.product( *[ range( p.N ) for p in n._parents ] ):

                total = LogVar( 0 )
                for i in range( n.N ):

                    total += messagePasser.conditionalParentChild( n, X, i )

                prob = LogVar( 1 )
                if( not isclose( total.logVal, prob.logVal ) ):
                    failString = '\n\nFailed the conditional parent child test for node: '+str( n ) + '\n'
                    failString += 'This is wrong: '+str( total ) + '\n'
                    failString += 'This is right: '+str( prob ) + '\n'
                    failString += 'up edge: '+str( n._upEdge ) + '\n'
                    failString += 'down edges: '+str( n._downEdges ) + '\n'
                    allFailures.append( failString )


    if( len( allFailures ) == 0 ):
        print( 'Passed all of the marginalization tests!' )
    else:
        for failure in allFailures:
            print(failure)
