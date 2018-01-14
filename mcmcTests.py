from Mendel import AutosomalMendelModel
from HGTest import marginalizeTest
import numpy as np

""" State marginalization test """

def marginalizeTests( model ):
    for graph in model.graphs:
        marginalizeTest( graph, printStuff=False )

""" Ratio test """

def ratioTest( model ):
    # P( Y | X, θ ) / P( Y' | X, θ ) =
    # P( Y, X, θ ) / P( Y', X, θ )

    model.resample()

    model.useFakeY()

    log_cond1 = 0
    log_joint1 = 0
    for i, graph in enumerate( model.graphs ):
        log_cond1  += graph.log_likelihood()
        log_joint1 += graph.log_joint()
        for j, root in enumerate( graph.roots ):
            log_joint1 += model.rootDists[ i ][ j ].log_likelihood()

    for graph in model.graphs:
        graph.resampleEmissions()

    log_cond2 = 0
    log_joint2 = 0
    for i, graph in enumerate( model.graphs ):
        log_cond2  += graph.log_likelihood()
        log_joint2 += graph.log_joint()
        for j, root in enumerate( graph.roots ):
            log_joint2 += model.rootDists[ i ][ j ].log_likelihood()


    if( not np.isclose( log_cond1 - log_cond2, log_joint1 - log_joint2 ) ):
        print('Failed the ratio test!')
        print('log P( X, π | Y ) = %f'%log_cond1)
        print('log P( X\', π\' | Y ) = %f'%log_cond2)
        print('\n')
        print('log P( X, π, Y ) = %f'%log_joint1)
        print('log P( X\', π\', Y ) = %f'%log_joint2)
        print('\n')
        print('log_cond1 - log_cond2 = %f'%( log_cond1 - log_cond2 ))
        print('log_joint1 - log_joint2 = %f'%( log_joint1 - log_joint2 ))
        print('\n')
        assert 0

    print('Passed a ratio test!')
    return log_cond1 - log_cond2 == log_joint1 - log_joint2

""" Gewek test """


def _forwardSample( model ):
    # Sample from P( θ )
    for i, graph in enumerate( model.graphs ):
        for j, root in enumerate( graph.roots ):
            model.rootDist[ i ][ j ].resample()

    # Sample from P( X, Y | θ )
    for graph in model.graphs:
        graph.resampleGraphStatesAndEmissions()

def _fullGibbsSample( model ):
    # Sample from P( X, θ | Y )
    model.resample()

    # Sample from P( Y | X, θ )
    for graph in model.graphs:
        graph.resampleEmissions()

def _collectStats( model ):
    assert 0

def _compareStats( model, stats1, stats2 ):
    assert 0

def gewekeTest( model ):
    # Sample from P( X, θ, Y ) using two methods:
    #  1. P( θ ), then P( X, Y | θ )
    #  2. P( X, θ | Y ), then P( Y | X, θ )
    # And make sure that they agree


    nIters = 1000

    stats1 = []
    for i in range( nIters ):
        model._forwardSample()
        stats1.append( model._collectStats() )

    stats2 = []
    for i in range( nIters ):
        model._fullGibbsSample()
        stats2.append( model._collectStats() )

    model._compareStats( stats1, stats2 )