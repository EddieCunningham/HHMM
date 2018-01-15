import numpy as np

def generic2DParameters(hyperGraph):

    transDist = np.array([ \
                             [ \
                                 [ \
                                   [[
                                    [[0.99,0.01],[0.5,0.5]], \
                                     [[0.25,0.75],[0.5,0.5]]], \
                                    [[[0.6,0.4],[0.5,0.5]], \
                                     [[0.01,0.99],[0.2,0.8]]]], \
                                   [[
                                    [[0.23,0.77],[0.1,0.9]], \
                                     [[0.8,0.2],[0.5,0.5]]], \
                                    [[[0.65,0.35],[0.5,0.5]], \
                                     [[0.11,0.89],[0.5,0.5]]]] \
                                ], \
                                [ \
                                   [[
                                    [[0.99,0.01],[0.5,0.5]], \
                                     [[0.7,0.3],[0.5,0.5]]], \
                                    [[[0.6,0.4],[0.5,0.5]], \
                                     [[0.01,0.99],[0.2,0.8]]]], \
                                   [[
                                    [[0.23,0.77],[0.1,0.9]], \
                                     [[0.8,0.2],[0.5,0.5]]], \
                                    [[[0.65,0.35],[0.5,0.5]], \
                                     [[0.11,0.89],[0.5,0.5]]] \
                                   ] \
                                ] \
                            ], \
                            [ \
                                 [ \
                                   [[
                                    [[0.99,0.01],[0.6,0.4]], \
                                     [[0.2,0.8],[0.2,0.8]]], \
                                    [[[0.6,0.4],[0.5,0.5]], \
                                     [[0.01,0.99],[0.2,0.8]]]], \
                                   [[
                                    [[0.99,0.01],[0.5,0.5]], \
                                     [[0.8,0.2],[0.5,0.5]]], \
                                    [[[0.65,0.35],[0.5,0.5]], \
                                     [[0.11,0.89],[0.5,0.5]]]] \
                                ], \
                                [ \
                                   [[
                                    [[0.23,0.77],[0.1,0.9]], \
                                     [[0.2,0.8],[0.5,0.5]]], \
                                    [[[0.6,0.4],[0.5,0.5]], \
                                     [[0.01,0.99],[0.2,0.8]]]], \
                                   [[
                                    [[0.29,0.71],[0.1,0.9]], \
                                     [[0.8,0.2],[0.5,0.5]]], \
                                    [[[0.65,0.35],[0.1,0.9]], \
                                     [[0.11,0.89],[0.25,0.75]]] \
                                   ] \
                                ] \
                            ] \
                        ])

    def rootDistFunction(person,i):
        theDist = np.array([[0.25,0.75] for _ in hyperGraph.roots])
        index = sorted(hyperGraph.roots).index(person)
        return theDist[index][i]

    def transDistFunction(parents,child,X,i):

        theDist = transDist
        for j in range(7-len(parents)-1):
            theDist = theDist[0]

        for j in range(len(parents)):
            theDist = theDist[X[j]]

        return theDist[i]

    def emissionDistFunction(person,i):
        theDist = np.array([[0.1],[0.9]])
        return theDist[i][person.y]

    return {
        'rootDist':rootDistFunction,
        'transDist':transDistFunction,
        'emissionDist':emissionDistFunction
    }


def autosomalHyperParameters():


    A_hyper = np.array([
            [
                [ 1.0, 0.5, 0.5, 0.0 ],
                [ 0.5, .25, .25, 0.0 ],
                [ 0.5, .25, .25, 0.0 ],
                [ 0.0, 0.0, 0.0, 0.0 ]
            ],
            [
                [ 0.0, 0.5, 0.5, 1.0 ],
                [ 0.0, .25, .25, 0.5 ],
                [ 0.0, .25, .25, 0.5 ],
                [ 0.0, 0.0, 0.0, 0.0 ]
            ],
            [
                [ 0.0, 0.0, 0.0, 0.0 ],
                [ 0.5, .25, .25, 0.0 ],
                [ 0.5, .25, .25, 0.0 ],
                [ 1.0, 0.5, 0.5, 0.0 ]
            ],
            [
                [ 0.0, 0.0, 0.0, 0.0 ],
                [ 0.0, .25, .25, 0.5 ],
                [ 0.0, .25, .25, 0.5 ],
                [ 0.0, 0.5, 0.5, 1.0 ]
            ]
        ])

    A_hyper = np.array( [ [ [ A_hyper[ k ][ i ][ j ] for k in range( A_hyper.shape[ 2 ] ) ] \
                                                     for j in range( A_hyper.shape[ 1 ] ) ] \
                                                     for i in range( A_hyper.shape[ 0 ] ) ] )

    L_hyper = np.array( [ [1.,0.],
                          [1.,0.],
                          [1.,0.],
                          [0.,1.] ] )

    A_hyper = 1. + ( 5 * A_hyper )
    L_hyper = 1. + ( 5 * L_hyper )
    pi_hyper = np.ones( 4 )

    return ( A_hyper, L_hyper, pi_hyper )
