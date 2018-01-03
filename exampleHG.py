from AutosomalDistribution import generic2DParameters, sampleHyperGraphParameters
from HHMMMessagePasser import HiddenMarkovModelMessagePasser
from HHMMHG import MessagePassingHG
from MMUpDown import MarkovModelMessagePasser

def cycleExample1( hiddenStateSize=2, isHidden=True ):
    hg = MessagePassingHG( hiddenStateSize )
    n0 = hg.addNode( 0 )
    n1 = hg.addNode( 1 )
    n2 = hg.addNode( 2 )
    n3 = hg.addNode( 3 )
    n7 = hg.addNode( 7 )

    e0 = hg.addEdge( set( [ n1 ] ), 0 )
    e0.addChild( n0 )

    e1 = hg.addEdge( set( [ n1 ] ), 1 )
    e1.addChild( n2 )
    e1.addChild( n3 )

    e5 = hg.addEdge( set( [ n2, n3, n1, n0 ] ), 5 )
    e5.addChild( n7 )

    hg.initialize()
    if( isHidden ):
        parameters   = generic2DParameters( hg )
        transFunc    = parameters[ 'transDist' ]
        emissionFunc = parameters[ 'emissionDist' ]
        rootFunc     = parameters[ 'rootDist' ]
        hg.setParameters( transFunc, emissionFunc, rootFunc )
    else:
        msg = MarkovModelMessagePasser( hg, generic2DParameters )
    return hg

def cycleExample2( hiddenStateSize=2, isHidden=True ):
    hg = MessagePassingHG( hiddenStateSize )
    n1 = hg.addNode( 1 )
    n3 = hg.addNode( 3 )
    n4 = hg.addNode( 4 )
    n5 = hg.addNode( 5 )
    n7 = hg.addNode( 7 )

    e1 = hg.addEdge( set( [ n1 ] ), 1 )
    e1.addChild( n3 )

    e2 = hg.addEdge( set( [ n3 ] ), 2 )
    e2.addChild( n4 )

    e3 = hg.addEdge( set( [ n4, n3 ] ), 3 )
    e3.addChild( n5 )

    e4 = hg.addEdge( set( [ n5, n1, n4 ] ), 4 )
    e4.addChild( n7 )

    hg.initialize()

    if( isHidden ):
        parameters   = generic2DParameters( hg )
        transFunc    = parameters[ 'transDist' ]
        emissionFunc = parameters[ 'emissionDist' ]
        rootFunc     = parameters[ 'rootDist' ]
        hg.setParameters( transFunc, emissionFunc, rootFunc )
    else:
        msg = MarkovModelMessagePasser( hg, generic2DParameters )
    return hg

def cycleExample3( hiddenStateSize=2, isHidden=True ):
    hg = MessagePassingHG( hiddenStateSize )
    n0 = hg.addNode( 0 )
    n1 = hg.addNode( 1 )
    n2 = hg.addNode( 2 )
    n3 = hg.addNode( 3 )
    n4 = hg.addNode( 4 )
    n5 = hg.addNode( 5 )
    n6 = hg.addNode( 6 )
    n7 = hg.addNode( 7 )
    n8 = hg.addNode( 8 )

    e1 = hg.addEdge( set( [ n1, n0 ] ), 1 )
    e1.addChild( n2 )
    e1.addChild( n3 )

    e2 = hg.addEdge( set( [ n2, n3 ] ), 2 )
    e2.addChild( n4 )

    e3 = hg.addEdge( set( [ n4, n8, n0, n2, n3 ] ), 3 )
    e3.addChild( n5 )
    e3.addChild( n6 )

    e4 = hg.addEdge( set( [ n5, n6, n1, n0, n2, n4, n6 ] ), 4 )
    e4.addChild( n7 )

    hg.initialize()

    if( isHidden ):
        parameters   = generic2DParameters( hg )
        transFunc    = parameters[ 'transDist' ]
        emissionFunc = parameters[ 'emissionDist' ]
        rootFunc     = parameters[ 'rootDist' ]
        hg.setParameters( transFunc, emissionFunc, rootFunc )
    else:
        msg = MarkovModelMessagePasser( hg, generic2DParameters )
    return hg

def cycleExample4( hiddenStateSize=2, isHidden=True ):
    hg = MessagePassingHG( hiddenStateSize )
    n1 = hg.addNode( 1 )
    n2 = hg.addNode( 2 )
    n4 = hg.addNode( 4 )
    n7 = hg.addNode( 7 )
    n8 = hg.addNode( 8 )

    e1 = hg.addEdge( set( [ n1 ] ), 1 )
    e1.addChild( n2 )

    e2 = hg.addEdge( set( [ n2, n1 ] ), 2 )
    e2.addChild( n4 )

    e3 = hg.addEdge( set( [ n4 ] ), 3 )
    e3.addChild( n7 )

    e4 = hg.addEdge( set( [ n7, n2 ] ), 4 )
    e4.addChild( n8 )

    hg.initialize()

    if( isHidden ):
        parameters   = generic2DParameters( hg )
        transFunc    = parameters[ 'transDist' ]
        emissionFunc = parameters[ 'emissionDist' ]
        rootFunc     = parameters[ 'rootDist' ]
        hg.setParameters( transFunc, emissionFunc, rootFunc )
    else:
        msg = MarkovModelMessagePasser( hg, generic2DParameters )
    return hg

def cycleExample5( hiddenStateSize=2, isHidden=True ):
    hg = MessagePassingHG( hiddenStateSize )
    n0 = hg.addNode( 0 )
    n1 = hg.addNode( 1 )
    n2 = hg.addNode( 2 )
    n3 = hg.addNode( 3 )
    n4 = hg.addNode( 4 )
    n5 = hg.addNode( 5 )
    n6 = hg.addNode( 6 )
    n7 = hg.addNode( 7 )
    n8 = hg.addNode( 8 )
    n9 = hg.addNode( 9 )
    n10 = hg.addNode( 10 )

    e1 = hg.addEdge( set( [ n1, n0 ] ), 1 )
    e1.addChild( n2 )
    e1.addChild( n3 )
    e1.addChild( n9 )

    e2 = hg.addEdge( set( [ n2, n3 ] ), 2 )
    e2.addChild( n4 )

    e3 = hg.addEdge( set( [ n4, n8, n2, n3, n1 ] ), 3 )
    e3.addChild( n5 )
    e3.addChild( n6 )

    e4 = hg.addEdge( set( [ n5, n6, n1, n2, n4, n6, n9 ] ), 4 )
    e4.addChild( n7 )
    e4.addChild( n10 )

    hg.initialize()

    if( isHidden ):
        parameters   = generic2DParameters( hg )
        transFunc    = parameters[ 'transDist' ]
        emissionFunc = parameters[ 'emissionDist' ]
        rootFunc     = parameters[ 'rootDist' ]
        hg.setParameters( transFunc, emissionFunc, rootFunc )
    else:
        msg = MarkovModelMessagePasser( hg, generic2DParameters )
    return hg


def cycleExample5_1( hiddenStateSize=2, isHidden=True ):
    hg = MessagePassingHG( hiddenStateSize )
    n1 = hg.addNode( 1 )
    n2 = hg.addNode( 2 )
    n3 = hg.addNode( 3 )
    n4 = hg.addNode( 4 )
    n5 = hg.addNode( 5 )
    n6 = hg.addNode( 6 )
    n7 = hg.addNode( 7 )
    n8 = hg.addNode( 8 )
    n9 = hg.addNode( 9 )

    e1 = hg.addEdge( set( [ n1 ] ), 1 )
    e1.addChild( n2 )
    e1.addChild( n3 )
    e1.addChild( n9 )

    e2 = hg.addEdge( set( [ n2, n3 ] ), 2 )
    e2.addChild( n4 )

    e3 = hg.addEdge( set( [ n4, n8, n2, n3, n1 ] ), 3 )
    e3.addChild( n5 )
    e3.addChild( n6 )

    e4 = hg.addEdge( set( [ n5, n6, n1, n2, n4, n6, n9 ] ), 4 )
    e4.addChild( n7 )

    hg.initialize()

    if( isHidden ):
        parameters   = generic2DParameters( hg )
        transFunc    = parameters[ 'transDist' ]
        emissionFunc = parameters[ 'emissionDist' ]
        rootFunc     = parameters[ 'rootDist' ]
        hg.setParameters( transFunc, emissionFunc, rootFunc )
    else:
        msg = MarkovModelMessagePasser( hg, generic2DParameters )
    return hg

def cycleExample6( hiddenStateSize=2, isHidden=True ):
    hg = MessagePassingHG( hiddenStateSize )
    n1 = hg.addNode( 1 )
    n2 = hg.addNode( 2 )
    n3 = hg.addNode( 3 )
    n4 = hg.addNode( 4 )

    e1 = hg.addEdge( set( [ n1 ] ), 1 )
    e1.addChild( n2 )
    e1.addChild( n3 )

    e2 = hg.addEdge( set( [ n2, n3 ] ), 2 )
    e2.addChild( n4 )

    hg.initialize()

    if( isHidden ):
        parameters   = generic2DParameters( hg )
        transFunc    = parameters[ 'transDist' ]
        emissionFunc = parameters[ 'emissionDist' ]
        rootFunc     = parameters[ 'rootDist' ]
        hg.setParameters( transFunc, emissionFunc, rootFunc )
    else:
        msg = MarkovModelMessagePasser( hg, generic2DParameters )
    return hg

def cycleExample7( hiddenStateSize=2, isHidden=True ):
    hg = MessagePassingHG( hiddenStateSize )
    n0 = hg.addNode( 0 )
    n1 = hg.addNode( 1 )
    n2 = hg.addNode( 2 )
    n3 = hg.addNode( 3 )
    n4 = hg.addNode( 4 )
    n5 = hg.addNode( 5 )

    e0 = hg.addEdge( set( [ n0 ] ), 0 )
    e0.addChild( n1 )

    e1 = hg.addEdge( set( [ n1 ] ), 1 )
    e1.addChild( n2 )
    e1.addChild( n3 )

    e2 = hg.addEdge( set( [ n2, n3 ] ), 2 )
    e2.addChild( n4 )

    e3 = hg.addEdge( set( [ n4, n0 ] ), 3 )
    e3.addChild( n5 )

    hg.initialize()

    if( isHidden ):
        parameters   = generic2DParameters( hg )
        transFunc    = parameters[ 'transDist' ]
        emissionFunc = parameters[ 'emissionDist' ]
        rootFunc     = parameters[ 'rootDist' ]
        hg.setParameters( transFunc, emissionFunc, rootFunc )
    else:
        msg = MarkovModelMessagePasser( hg, generic2DParameters )
    return hg

def cycleExample8( hiddenStateSize=2, isHidden=True ):
    hg = MessagePassingHG( hiddenStateSize )
    n0 = hg.addNode( 0 )
    n2 = hg.addNode( 2 )
    n3 = hg.addNode( 3 )
    n4 = hg.addNode( 4 )


    e1 = hg.addEdge( set( [ n0 ] ), 1 )
    e1.addChild( n2 )
    e1.addChild( n3 )

    e2 = hg.addEdge( set( [ n2, n3, n0 ] ), 2 )
    e2.addChild( n4 )

    hg.initialize()

    if( isHidden ):
        parameters   = generic2DParameters( hg )
        transFunc    = parameters[ 'transDist' ]
        emissionFunc = parameters[ 'emissionDist' ]
        rootFunc     = parameters[ 'rootDist' ]
        hg.setParameters( transFunc, emissionFunc, rootFunc )
    else:
        msg = MarkovModelMessagePasser( hg, generic2DParameters )
    return hg

def cycleExample9( hiddenStateSize=2, isHidden=True ):
    hg = MessagePassingHG( hiddenStateSize )
    n0 = hg.addNode( 0 )
    n1 = hg.addNode( 1 )
    n2 = hg.addNode( 2 )
    n3 = hg.addNode( 3 )

    e1 = hg.addEdge( set( [ n0 ] ), 1 )
    e1.addChild( n1 )
    e1.addChild( n2 )

    e2 = hg.addEdge( set( [ n0, n1, n2 ] ), 2 )
    e2.addChild( n3 )

    hg.initialize()

    if( isHidden ):
        parameters   = generic2DParameters( hg )
        transFunc    = parameters[ 'transDist' ]
        emissionFunc = parameters[ 'emissionDist' ]
        rootFunc     = parameters[ 'rootDist' ]
        hg.setParameters( transFunc, emissionFunc, rootFunc )
    else:
        msg = MarkovModelMessagePasser( hg, generic2DParameters )
    return hg


def cycleExample10( hiddenStateSize=2, isHidden=True ):
    # the current fbs algorithm makes this
    # disjoint!!!!
    hg = MessagePassingHG( hiddenStateSize )
    n0 = hg.addNode( 0 )
    n1 = hg.addNode( 1 )
    n2 = hg.addNode( 2 )
    n3 = hg.addNode( 3 )
    n4 = hg.addNode( 4 )
    n5 = hg.addNode( 5 )
    n6 = hg.addNode( 6 )
    n7 = hg.addNode( 7 )
    n8 = hg.addNode( 8 )
    n9 = hg.addNode( 9 )
    n10 = hg.addNode( 10 )
    n11 = hg.addNode( 11 )
    n12 = hg.addNode( 12 )
    n13 = hg.addNode( 13 )
    n14 = hg.addNode( 14 )
    n15 = hg.addNode( 15 )
    n16 = hg.addNode( 16 )
    n17 = hg.addNode( 17 )
    n18 = hg.addNode( 18 )
    n19 = hg.addNode( 19 )
    n20 = hg.addNode( 20 )

    e1 = hg.addEdge( set( [ n1, n0 ] ), 1 )
    e1.addChild( n2 )
    e1.addChild( n3 )
    e1.addChild( n9 )

    e2 = hg.addEdge( set( [ n2, n3 ] ), 2 )
    e2.addChild( n4 )

    e3 = hg.addEdge( set( [ n4, n8, n2, n3, n1 ] ), 3 )
    e3.addChild( n5 )
    e3.addChild( n6 )

    e4 = hg.addEdge( set( [ n5, n6, n1, n2, n4, n6, n9 ] ), 4 )
    e4.addChild( n7 )
    e4.addChild( n10 )

    e5 = hg.addEdge( set( [ n10 ] ), 5 )
    e5.addChild( n11 )

    e6 = hg.addEdge( set( [ n11 ] ), 6 )
    e6.addChild( n12 )

    e7 = hg.addEdge( set( [ n12 ] ), 7 )
    e7.addChild( n13 )

    e8 = hg.addEdge( set( [ n13 ] ), 8 )
    e8.addChild( n14 )

    e9 = hg.addEdge( set( [ n14, n0 ] ), 9 )
    e9.addChild( n15 )
    e9.addChild( n16 )
    e9.addChild( n17 )

    e10 = hg.addEdge( set( [ n1, n2, n3 ] ), 10 )
    e10.addChild( n18 )

    e11 = hg.addEdge( set( [ n11, n12, n13 ] ), 11 )
    e11.addChild( n19 )

    e12 = hg.addEdge( set( [ n9, n19 ] ), 12 )
    e12.addChild( n20 )

    hg.initialize()

    if( isHidden ):
        parameters   = generic2DParameters( hg )
        transFunc    = parameters[ 'transDist' ]
        emissionFunc = parameters[ 'emissionDist' ]
        rootFunc     = parameters[ 'rootDist' ]
        hg.setParameters( transFunc, emissionFunc, rootFunc )
    else:
        msg = MarkovModelMessagePasser( hg, generic2DParameters )
    return hg


def cycleExample11( hiddenStateSize=2, isHidden=True ):
    hg = MessagePassingHG( hiddenStateSize )
    n0 = hg.addNode( 0 )
    n1 = hg.addNode( 1 )
    n2 = hg.addNode( 2 )
    n3 = hg.addNode( 3 )
    n4 = hg.addNode( 4 )
    n5 = hg.addNode( 5 )
    n6 = hg.addNode( 6 )
    n7 = hg.addNode( 7 )
    n8 = hg.addNode( 8 )
    n9 = hg.addNode( 9 )
    n10 = hg.addNode( 10 )
    n11 = hg.addNode( 11 )
    n12 = hg.addNode( 12 )
    n13 = hg.addNode( 13 )
    n14 = hg.addNode( 14 )
    n15 = hg.addNode( 15 )
    n16 = hg.addNode( 16 )
    n17 = hg.addNode( 17 )

    e1 = hg.addEdge( set( [ n1, n0 ] ), 1 )
    e1.addChild( n2 )
    e1.addChild( n3 )
    e1.addChild( n9 )

    e2 = hg.addEdge( set( [ n2, n3 ] ), 2 )
    e2.addChild( n4 )

    e3 = hg.addEdge( set( [ n4, n8, n2, n3, n1 ] ), 3 )
    e3.addChild( n5 )
    e3.addChild( n6 )

    e4 = hg.addEdge( set( [ n5, n6, n1, n2, n4, n6, n9 ] ), 4 )
    e4.addChild( n7 )
    e4.addChild( n10 )

    e5 = hg.addEdge( set( [ n10 ] ), 5 )
    e5.addChild( n11 )

    e6 = hg.addEdge( set( [ n11 ] ), 6 )
    e6.addChild( n12 )

    e7 = hg.addEdge( set( [ n12 ] ), 7 )
    e7.addChild( n13 )

    e8 = hg.addEdge( set( [ n13 ] ), 8 )
    e8.addChild( n14 )

    e9 = hg.addEdge( set( [ n14, n0 ] ), 9 )
    e9.addChild( n15 )
    e9.addChild( n16 )
    e9.addChild( n17 )

    hg.initialize()

    if( isHidden ):
        parameters   = generic2DParameters( hg )
        transFunc    = parameters[ 'transDist' ]
        emissionFunc = parameters[ 'emissionDist' ]
        rootFunc     = parameters[ 'rootDist' ]
        hg.setParameters( transFunc, emissionFunc, rootFunc )
    else:
        msg = MarkovModelMessagePasser( hg, generic2DParameters )
    return hg

def cycleExample12( hiddenStateSize=2, isHidden=True ):
    hg = MessagePassingHG( hiddenStateSize )
    n0 = hg.addNode( 0 )
    n2 = hg.addNode( 2 )
    n3 = hg.addNode( 3 )
    n5 = hg.addNode( 5 )
    n6 = hg.addNode( 6 )

    e1 = hg.addEdge( set( [ n0 ] ), 1 )
    e1.addChild( n2 )

    e2 = hg.addEdge( set( [ n2 ] ), 2 )
    e2.addChild( n3 )

    e3 = hg.addEdge( set( [ n3 ] ), 3 )
    e3.addChild( n5 )

    e4 = hg.addEdge( set( [ n5, n0 ] ), 4 )
    e4.addChild( n6 )

    hg.initialize()

    if( isHidden ):
        parameters   = generic2DParameters( hg )
        transFunc    = parameters[ 'transDist' ]
        emissionFunc = parameters[ 'emissionDist' ]
        rootFunc     = parameters[ 'rootDist' ]
        hg.setParameters( transFunc, emissionFunc, rootFunc )
    else:
        msg = MarkovModelMessagePasser( hg, generic2DParameters )
    return hg

def disjointExample( hiddenStateSize=2, isHidden=True ):
    hg = MessagePassingHG( hiddenStateSize )

    n1 = hg.addNode( 1 )
    n2 = hg.addNode( 2 )
    n3 = hg.addNode( 3 )
    n4 = hg.addNode( 4 )
    n5 = hg.addNode( 5 )
    # n6 = hg.addNode( 6 )
    n7 = hg.addNode( 7 )
    n8 = hg.addNode( 8 )

    n10 = hg.addNode( 10 )

    e1 = hg.addEdge( set( [ n1, n2 ] ), 1 )
    e1.addChild( n3 )
    e1.addChild( n4 )

    e2 = hg.addEdge( set( [ n4 ] ), 2 )
    e2.addChild( n5 )
    # e2.addChild( n6 )

    e3 = hg.addEdge( set( [ n2, n3, n5 ] ), 3 )
    e3.addChild( n7 )

    e4 = hg.addEdge( set( [ n1, n4 ] ), 4 )
    e4.addChild( n8 )

    e5 = hg.addEdge( set( [ n4 ] ), 5 )
    e5.addChild( n10 )

    hg.initialize()

    if( isHidden ):
        parameters   = generic2DParameters( hg )
        transFunc    = parameters[ 'transDist' ]
        emissionFunc = parameters[ 'emissionDist' ]
        rootFunc     = parameters[ 'rootDist' ]
        hg.setParameters( transFunc, emissionFunc, rootFunc )
    else:
        msg = MarkovModelMessagePasser( hg, generic2DParameters )
    return hg

def nonCycle( hiddenStateSize=2, isHidden=True ):
    hg = MessagePassingHG( hiddenStateSize )
    n0 = hg.addNode( 0 )
    n1 = hg.addNode( 1 )
    n2 = hg.addNode( 2 )
    n3 = hg.addNode( 3 )

    e1 = hg.addEdge( set( [ n0 ] ), 1 )
    e1.addChild( n1 )
    e1.addChild( n2 )

    e2 = hg.addEdge( set( [ n2 ] ), 2 )
    e2.addChild( n3 )

    hg.initialize()

    if( isHidden ):
        parameters   = generic2DParameters( hg )
        transFunc    = parameters[ 'transDist' ]
        emissionFunc = parameters[ 'emissionDist' ]
        rootFunc     = parameters[ 'rootDist' ]
        hg.setParameters( transFunc, emissionFunc, rootFunc )
    else:
        msg = MarkovModelMessagePasser( hg, generic2DParameters )
    return hg
