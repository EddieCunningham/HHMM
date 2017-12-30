import scipy.stats
import numpy as np
import os
import time
from scipy.special import digamma
import matplotlib.pyplot as plt
# from CythonCode.LogVarCode import LogVar
import json
import itertools
from model import Pedigree
from PedigreeHypergraph import PedigreeHG
from AutosomalDistribution import *
# from CythonCode.HHMMUpDownFast import HiddenMarkovModelMessagePasser,MessagePassingHG
from HHMMUpDown import HiddenMarkovModelMessagePasser,MessagePassingHG
# from MMUpDown import MarkovModelMessagePasser

def printListOfNodes(nodes):
    print('[ '),
    for node in nodes:
        print(str(node._id)+' '),
    print(']')

def pedigreeExample(name):
    pedigreeFolderName = '/Users/Eddie/kec-bot/app/pedigreeDataOLDBUTWORKS/'

    pedigreeNames = [name]

    allPedigrees = {}
    for name in pedigreeNames:

        filename = os.path.join(pedigreeFolderName,name+'.json')
        with open(filename) as data_file:
            data = json.loads(json.load(data_file))
        pedigree = Pedigree(data)
        """ This does nothing, just haven't gotten around
            to changing how model.py works """
        pedigree.setAffectedFunctions(lambda x: None)

        hg = PedigreeHG(name)
        hg.initialize(pedigree)
        allPedigrees[name] = hg

    hg = allPedigrees[name]

    hg.draw()
    hg.initHyperParams('autosome','dominant',1000)
    msg = HiddenMarkovModelMessagePasser(hg,sampleHyperGraphParameters)
    return hg,msg


def cycleExample1(isHidden):
    hg = MessagePassingHG(2)
    n0 = hg.addNode(0)
    n1 = hg.addNode(1)
    n2 = hg.addNode(2)
    n3 = hg.addNode(3)
    n7 = hg.addNode(7)

    e0 = hg.addEdge(set([n1]),0)
    e0.addChild(n0)

    e1 = hg.addEdge(set([n1]),1)
    e1.addChild(n2)
    e1.addChild(n3)

    e5 = hg.addEdge(set([n2,n3,n1,n0]),5)
    e5.addChild(n7)

    hg.initialize()

    hg.draw()
    if(isHidden):
        msg = HiddenMarkovModelMessagePasser(hg,generic2DParameters)
    else:
        msg = MarkovModelMessagePasser(hg,generic2DParameters)
    return hg,msg

def cycleExample2(isHidden):
    hg = MessagePassingHG(1)
    n1 = hg.addNode(1)
    n3 = hg.addNode(3)
    n4 = hg.addNode(4)
    n5 = hg.addNode(5)
    n7 = hg.addNode(7)

    e1 = hg.addEdge(set([n1]),1)
    e1.addChild(n3)

    e2 = hg.addEdge(set([n3]),2)
    e2.addChild(n4)

    e3 = hg.addEdge(set([n4,n3]),3)
    e3.addChild(n5)

    e4 = hg.addEdge(set([n5,n1,n4]),4)
    e4.addChild(n7)

    hg.initialize()

    hg.draw()

    if(isHidden):
        msg = HiddenMarkovModelMessagePasser(hg,generic2DParameters)
    else:
        msg = MarkovModelMessagePasser(hg,generic2DParameters)
    return hg,msg

def cycleExample3(isHidden):
    hg = MessagePassingHG(2)
    n0 = hg.addNode(0)
    n1 = hg.addNode(1)
    n2 = hg.addNode(2)
    n3 = hg.addNode(3)
    n4 = hg.addNode(4)
    n5 = hg.addNode(5)
    n6 = hg.addNode(6)
    n7 = hg.addNode(7)
    n8 = hg.addNode(8)

    e1 = hg.addEdge(set([n1,n0]),1)
    e1.addChild(n2)
    e1.addChild(n3)

    e2 = hg.addEdge(set([n2,n3]),2)
    e2.addChild(n4)

    e3 = hg.addEdge(set([n4,n8,n0,n2,n3]),3)
    e3.addChild(n5)
    e3.addChild(n6)

    e4 = hg.addEdge(set([n5,n6,n1,n0,n2,n4,n6]),4)
    e4.addChild(n7)

    hg.initialize()

    hg.draw()

    if(isHidden):
        msg = HiddenMarkovModelMessagePasser(hg,generic2DParameters)
    else:
        msg = MarkovModelMessagePasser(hg,generic2DParameters)
    return hg,msg

def cycleExample4(isHidden=True):
    hg = MessagePassingHG(2)
    n1 = hg.addNode(1)
    n2 = hg.addNode(2)
    n4 = hg.addNode(4)
    n7 = hg.addNode(7)
    n8 = hg.addNode(8)

    e1 = hg.addEdge(set([n1]),1)
    e1.addChild(n2)

    e2 = hg.addEdge(set([n2,n1]),2)
    e2.addChild(n4)

    e3 = hg.addEdge(set([n4]),3)
    e3.addChild(n7)

    e4 = hg.addEdge(set([n7,n2]),4)
    e4.addChild(n8)

    hg.initialize()

    hg.draw()

    if(isHidden):
        msg = HiddenMarkovModelMessagePasser(hg,generic2DParameters)
    else:
        msg = MarkovModelMessagePasser(hg,generic2DParameters)
    return hg,msg

def cycleExample5(isHidden=True):
    hg = MessagePassingHG(2)
    n0 = hg.addNode(0)
    n1 = hg.addNode(1)
    n2 = hg.addNode(2)
    n3 = hg.addNode(3)
    n4 = hg.addNode(4)
    n5 = hg.addNode(5)
    n6 = hg.addNode(6)
    n7 = hg.addNode(7)
    n8 = hg.addNode(8)
    n9 = hg.addNode(9)
    n10 = hg.addNode(10)

    e1 = hg.addEdge(set([n1,n0]),1)
    e1.addChild(n2)
    e1.addChild(n3)
    e1.addChild(n9)

    e2 = hg.addEdge(set([n2,n3]),2)
    e2.addChild(n4)

    e3 = hg.addEdge(set([n4,n8,n2,n3,n1]),3)
    e3.addChild(n5)
    e3.addChild(n6)

    e4 = hg.addEdge(set([n5,n6,n1,n2,n4,n6,n9]),4)
    e4.addChild(n7)
    e4.addChild(n10)

    hg.initialize()

    hg.draw()

    if(isHidden):
        msg = HiddenMarkovModelMessagePasser(hg,generic2DParameters)
    else:
        msg = MarkovModelMessagePasser(hg,generic2DParameters)
    return hg,msg


def cycleExample5_1(isHidden=True):
    hg = MessagePassingHG(2)
    n1 = hg.addNode(1)
    n2 = hg.addNode(2)
    n3 = hg.addNode(3)
    n4 = hg.addNode(4)
    n5 = hg.addNode(5)
    n6 = hg.addNode(6)
    n7 = hg.addNode(7)
    n8 = hg.addNode(8)
    n9 = hg.addNode(9)

    e1 = hg.addEdge(set([n1]),1)
    e1.addChild(n2)
    e1.addChild(n3)
    e1.addChild(n9)

    e2 = hg.addEdge(set([n2,n3]),2)
    e2.addChild(n4)

    e3 = hg.addEdge(set([n4,n8,n2,n3,n1]),3)
    e3.addChild(n5)
    e3.addChild(n6)

    e4 = hg.addEdge(set([n5,n6,n1,n2,n4,n6,n9]),4)
    e4.addChild(n7)

    hg.initialize()

    hg.draw()

    if(isHidden):
        msg = HiddenMarkovModelMessagePasser(hg,generic2DParameters)
    else:
        msg = MarkovModelMessagePasser(hg,generic2DParameters)
    return hg,msg

def cycleExample6(isHidden=True):
    hg = MessagePassingHG(2)
    n1 = hg.addNode(1)
    n2 = hg.addNode(2)
    n3 = hg.addNode(3)
    n4 = hg.addNode(4)

    e1 = hg.addEdge(set([n1]),1)
    e1.addChild(n2)
    e1.addChild(n3)

    e2 = hg.addEdge(set([n2,n3]),2)
    e2.addChild(n4)

    hg.initialize()

    hg.draw()

    if(isHidden):
        msg = HiddenMarkovModelMessagePasser(hg,generic2DParameters)
    else:
        msg = MarkovModelMessagePasser(hg,generic2DParameters)
    return hg,msg

def cycleExample7(isHidden=True):
    hg = MessagePassingHG(2)
    n0 = hg.addNode(0)
    n1 = hg.addNode(1)
    n2 = hg.addNode(2)
    n3 = hg.addNode(3)
    n4 = hg.addNode(4)
    n5 = hg.addNode(5)

    e0 = hg.addEdge(set([n0]),0)
    e0.addChild(n1)

    e1 = hg.addEdge(set([n1]),1)
    e1.addChild(n2)
    e1.addChild(n3)

    e2 = hg.addEdge(set([n2,n3]),2)
    e2.addChild(n4)

    e3 = hg.addEdge(set([n4,n0]),3)
    e3.addChild(n5)

    hg.initialize()

    hg.draw()

    if(isHidden):
        msg = HiddenMarkovModelMessagePasser(hg,generic2DParameters)
    else:
        msg = MarkovModelMessagePasser(hg,generic2DParameters)
    return hg,msg

def cycleExample8(isHidden=True):
    hg = MessagePassingHG(2)
    n0 = hg.addNode(0)
    n2 = hg.addNode(2)
    n3 = hg.addNode(3)
    n4 = hg.addNode(4)


    e1 = hg.addEdge(set([n0]),1)
    e1.addChild(n2)
    e1.addChild(n3)

    e2 = hg.addEdge(set([n2,n3,n0]),2)
    e2.addChild(n4)

    hg.initialize()

    hg.draw()

    if(isHidden):
        msg = HiddenMarkovModelMessagePasser(hg,generic2DParameters)
    else:
        msg = MarkovModelMessagePasser(hg,generic2DParameters)
    return hg,msg

def cycleExample9(isHidden=True):
    hg = MessagePassingHG(2)
    n0 = hg.addNode(0)
    n1 = hg.addNode(1)
    n2 = hg.addNode(2)
    n3 = hg.addNode(3)

    e1 = hg.addEdge(set([n0]),1)
    e1.addChild(n1)
    e1.addChild(n2)

    e2 = hg.addEdge(set([n0,n1,n2]),2)
    e2.addChild(n3)

    hg.initialize()

    hg.draw()

    if(isHidden):
        msg = HiddenMarkovModelMessagePasser(hg,generic2DParameters)
    else:
        msg = MarkovModelMessagePasser(hg,generic2DParameters)
    return hg,msg


def cycleExample10(isHidden=True):
    hg = MessagePassingHG(2)
    n0 = hg.addNode(0)
    n1 = hg.addNode(1)
    n2 = hg.addNode(2)
    n3 = hg.addNode(3)
    n4 = hg.addNode(4)
    n5 = hg.addNode(5)
    n6 = hg.addNode(6)
    n7 = hg.addNode(7)
    n8 = hg.addNode(8)
    n9 = hg.addNode(9)
    n10 = hg.addNode(10)
    n11 = hg.addNode(11)
    n12 = hg.addNode(12)
    n13 = hg.addNode(13)
    n14 = hg.addNode(14)
    n15 = hg.addNode(15)
    n16 = hg.addNode(16)
    n17 = hg.addNode(17)
    n18 = hg.addNode(18)
    n19 = hg.addNode(19)
    n20 = hg.addNode(20)

    e1 = hg.addEdge(set([n1,n0]),1)
    e1.addChild(n2)
    e1.addChild(n3)
    e1.addChild(n9)

    e2 = hg.addEdge(set([n2,n3]),2)
    e2.addChild(n4)

    e3 = hg.addEdge(set([n4,n8,n2,n3,n1]),3)
    e3.addChild(n5)
    e3.addChild(n6)

    e4 = hg.addEdge(set([n5,n6,n1,n2,n4,n6,n9]),4)
    e4.addChild(n7)
    e4.addChild(n10)

    e5 = hg.addEdge(set([n10]),5)
    e5.addChild(n11)

    e6 = hg.addEdge(set([n11]),6)
    e6.addChild(n12)

    e7 = hg.addEdge(set([n12]),7)
    e7.addChild(n13)

    e8 = hg.addEdge(set([n13]),8)
    e8.addChild(n14)

    e9 = hg.addEdge(set([n14,n0]),9)
    e9.addChild(n15)
    e9.addChild(n16)
    e9.addChild(n17)

    e10 = hg.addEdge(set([n1,n2,n3]),10)
    e10.addChild(n18)

    e11 = hg.addEdge(set([n11,n12,n13]),11)
    e11.addChild(n19)

    e12 = hg.addEdge(set([n9,n19]),12)
    e12.addChild(n20)

    hg.initialize()

    hg.draw()

    if(isHidden):
        msg = HiddenMarkovModelMessagePasser(hg,generic2DParameters)
    else:
        msg = MarkovModelMessagePasser(hg,generic2DParameters)
    return hg,msg


def cycleExample11(isHidden=True):
    hg = MessagePassingHG(2)
    n0 = hg.addNode(0)
    n1 = hg.addNode(1)
    n2 = hg.addNode(2)
    n3 = hg.addNode(3)
    n4 = hg.addNode(4)
    n5 = hg.addNode(5)
    n6 = hg.addNode(6)
    n7 = hg.addNode(7)
    n8 = hg.addNode(8)
    n9 = hg.addNode(9)
    n10 = hg.addNode(10)
    n11 = hg.addNode(11)
    n12 = hg.addNode(12)
    n13 = hg.addNode(13)
    n14 = hg.addNode(14)
    n15 = hg.addNode(15)
    n16 = hg.addNode(16)
    n17 = hg.addNode(17)

    e1 = hg.addEdge(set([n1,n0]),1)
    e1.addChild(n2)
    e1.addChild(n3)
    e1.addChild(n9)

    e2 = hg.addEdge(set([n2,n3]),2)
    e2.addChild(n4)

    e3 = hg.addEdge(set([n4,n8,n2,n3,n1]),3)
    e3.addChild(n5)
    e3.addChild(n6)

    e4 = hg.addEdge(set([n5,n6,n1,n2,n4,n6,n9]),4)
    e4.addChild(n7)
    e4.addChild(n10)

    e5 = hg.addEdge(set([n10]),5)
    e5.addChild(n11)

    e6 = hg.addEdge(set([n11]),6)
    e6.addChild(n12)

    e7 = hg.addEdge(set([n12]),7)
    e7.addChild(n13)

    e8 = hg.addEdge(set([n13]),8)
    e8.addChild(n14)

    e9 = hg.addEdge(set([n14,n0]),9)
    e9.addChild(n15)
    e9.addChild(n16)
    e9.addChild(n17)

    hg.initialize()

    hg.draw()

    if(isHidden):
        msg = HiddenMarkovModelMessagePasser(hg,generic2DParameters)
    else:
        msg = MarkovModelMessagePasser(hg,generic2DParameters)
    return hg,msg

def disjointExample(isHidden=True):
    hg = MessagePassingHG(2)

    n1 = hg.addNode(1)
    n2 = hg.addNode(2)
    n3 = hg.addNode(3)
    n4 = hg.addNode(4)
    n5 = hg.addNode(5)
    # n6 = hg.addNode(6)
    n7 = hg.addNode(7)
    n8 = hg.addNode(8)

    n10 = hg.addNode(10)

    e1 = hg.addEdge(set([n1,n2]),1)
    e1.addChild(n3)
    e1.addChild(n4)

    e2 = hg.addEdge(set([n4]),2)
    e2.addChild(n5)
    # e2.addChild(n6)

    e3 = hg.addEdge(set([n2,n3,n5]),3)
    e3.addChild(n7)

    e4 = hg.addEdge(set([n1,n4]),4)
    e4.addChild(n8)

    e5 = hg.addEdge(set([n4]),5)
    e5.addChild(n10)

    hg.initialize()

    hg.draw()

    if(isHidden):
        msg = HiddenMarkovModelMessagePasser(hg,generic2DParameters)
    else:
        msg = MarkovModelMessagePasser(hg,generic2DParameters)
    return hg,msg

def nonCycle(isHidden=True):
    hg = MessagePassingHG(1)
    n0 = hg.addNode(0)
    n1 = hg.addNode(1)
    n2 = hg.addNode(2)
    n3 = hg.addNode(3)

    e1 = hg.addEdge(set([n0]),1)
    e1.addChild(n1)
    e1.addChild(n2)

    e2 = hg.addEdge(set([n2]),2)
    e2.addChild(n3)

    hg.initialize()

    hg.draw()

    if(isHidden):
        msg = HiddenMarkovModelMessagePasser(hg,generic2DParameters)
    else:
        msg = MarkovModelMessagePasser(hg,generic2DParameters)
    return hg,msg


# hg,msg = pedigreeExample()
# hg,msg = cycleExample4(True)
# hg,msg = cycleExample5(True)
# hg,msg = cycleExample5_1(True)
# hg,msg = cycleExample6()
# hg,msg = cycleExample7(True)
# hg,msg = disjointExample()

# print('Feedback set: '+str(msg._feedbackSet))
# print('Feedback set: '+str(msg._feedbackStack))

# n4 = hg.addNode(4)
# val0,str0 = msg.fullBruteForce([n4],[0])
# val1,str1 = msg.fullBruteForce([n4],[1])
# print(val0+val1)
# print(str0)
# assert 0
def HMMTest(msg):
    start = time.time()
    msg.preprocess()
    end = time.time()
    print('\nPreprocess time: '+str(end-start))
    msg.aTest(True)
    start = time.time()
    for _ in range(1):
        msg.getStats()
    end = time.time()
    print('Traversal time: '+str(end-start))
    msg.aTest()

def MMTest():

    for node in msg.nodes:
        for i in range(node.N):
            comp = msg.getW(node,i)
            bf = LogVar(0)#msg.getWBruteForce([node],[i])
            print(comp)

# hg,msg = nonCycle(True)
# hg,msg = cycleExample2(True)
# hg,msg = cycleExample5(True)
# hg,msg = cycleExample9(True)
# hg,msg = cycleExample10(True)
hg,msg = cycleExample11(True)
hg.draw()
HMMTest(msg)
# HMMTest(cycleExample1(True)[1])
# HMMTest(cycleExample3(True)[1])
# HMMTest(cycleExample4(True)[1])
# HMMTest(cycleExample5(True)[1])
# HMMTest(cycleExample5_1(True)[1])
# HMMTest(cycleExample6(True)[1])
# HMMTest(cycleExample7(True)[1])
# HMMTest(pedigreeExample('3818J')[1])