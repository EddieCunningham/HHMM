
import scipy.stats
import numpy as np
from scipy.special import digamma
from numpy.random import RandomState
r = RandomState(5)

def uniformRootPrior(hyperGraph,strength=1):
    """ Generates a uniform for the dirichlet distribution """
    ans = []
    for r in sorted(hyperGraph._roots):
        ans.append([1. for _ in range(r.N)])
    return strength*np.array(ans)

def autosomalTransPrior(hyperGraph,strength=1):
    """ This is the prior for autosomal inheritance """
    """ Generates the dirichlet prior for the transition tensor.
        Because the transition tensor is (N,N,N), will return a
        Dirichlet prior over each NxN slice of the tensor """
    g = np.array([
            [
                [1.0,0.5,0.5,0.0],
                [0.5,.25,.25,0.0],
                [0.5,.25,.25,0.0],
                [0.0,0.0,0.0,0.0]
            ],
            [
                [0.0,0.5,0.5,1.0],
                [0.0,.25,.25,0.5],
                [0.0,.25,.25,0.5],
                [0.0,0.0,0.0,0.0]
            ],
            [
                [0.0,0.0,0.0,0.0],
                [0.5,.25,.25,0.0],
                [0.5,.25,.25,0.0],
                [1.0,0.5,0.5,0.0]
            ],
            [
                [0.0,0.0,0.0,0.0],
                [0.0,.25,.25,0.5],
                [0.0,.25,.25,0.5],
                [0.0,0.5,0.5,1.0]
            ]
        ])

    """ Transform matrix so that it is indexed by [parentA][parentB][child] """
    transformedG = np.array([[[g[k][i][j] for k in range(g.shape[2])] for j in range(g.shape[1])] for i in range(g.shape[0])])
    return 1.+(strength*transformedG)

def dominantEmissionPrior(hyperGraph,strength=1):
    """ Prior for dominant emission scheme """
    return strength*(1+np.array([
            [0.,0.],
            [0.,0.],
            [0.,0.],
            [0.,1.]
        ]))

# def hyperGraphBaseHyperParameters(hyperGraph,strength=1):
#     """ Will return base hyperparameters for the hypergraph """
#     alpha_r = uniformRootPrior(hyperGraph,strength=strength)
#     alpha_t = autosomalTransPrior(hyperGraph,strength=strength)
#     alpha_e = dominantEmissionPrior(hyperGraph,strength=strength)
#     return {
#         'rootPriors':alpha_r,
#         'transPriors':alpha_t,
#         'emissionPriors':alpha_e
#     }

def generateRootDist(hyperGraph,alpha_r):
    """ Returns a function that will sample root distribution parameters """
    theDist = np.array([r.dirichlet(a) for a in alpha_r])

    def rootDistFunction(person,i):
        index = sorted(hyperGraph._roots).index(person)
        return theDist[index][i]
    return rootDistFunction

def generateTransDist(hyperGraph,alpha_t):
    """ Returns a function that will sample transition distribution parameters """
    theDist = np.array([[r.dirichlet(row) for row in mat] for mat in alpha_t])

    def transDistFunction(parents,child,X,i):
        return theDist[X[0]][X[1]][i]
    return transDistFunction

def generateEmissionDist(hyperGraph,alpha_e):
    """ Returns a function that will sample emission distribution parameters """
    theDist = np.array([np.random.dirichlet(a) for a in alpha_e])

    def emissionDistFunction(person,i):
        return theDist[i][person.y]
    return emissionDistFunction

def hyperGraphBaseHyperParameters(hyperGraph,transType,emissionType,strength=1):
    """ Will return base hyperparameters for the hypergraph """
    alpha_r = uniformRootPrior(hyperGraph,strength=strength)

    if(transType == 'autosome'):
        alpha_t = autosomalTransPrior(hyperGraph,strength=strength)
    else:
        assert 0

    if(emissionType == 'dominant'):
        alpha_e = dominantEmissionPrior(hyperGraph,strength=strength)
    elif(emissionType == 'recessive'):
        alpha_e = recessiveEmissionPrior(hyperGraph,strength=strength)
    else:
        assert 0

    return {
        'rootPriors':alpha_r,
        'transPriors':alpha_t,
        'emissionPriors':alpha_e
    }

def sampleHyperGraphParameters(hyperGraph):
    """ Will sample parameters of hmm from prior distribution """
    rootDistribution = generateRootDist(hyperGraph,hyperGraph._hyperParams['rootPriors'])
    transitionDistribution = generateTransDist(hyperGraph,hyperGraph._hyperParams['transPriors'])
    emissionDistribution = generateEmissionDist(hyperGraph,hyperGraph._hyperParams['emissionPriors'])
    return {
        'rootDist':rootDistribution,
        'transDist':transitionDistribution,
        'emissionDist':emissionDistribution
    }


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
    # transDist = np.array([ \
    #                          [ \
    #                              [ \
    #                                [[
    #                                  [[1.0,0.0],[1.0,0.0]], \
    #                                  [[0.0,1.0],[0.0,1.0]]], \
    #                                 [[[1.0,0.0],[1.0,0.0]], \
    #                                  [[0.0,1.0],[0.0,1.0]]]] \
    #                             ]
    #                         ]
    #                     ])


    def rootDistFunction(person,i):
        theDist = np.array([[0.25,0.75] for _ in hyperGraph._roots])
        index = sorted(hyperGraph._roots).index(person)
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
