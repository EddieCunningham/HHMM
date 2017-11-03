import scipy.stats
import numpy as np

import time
from scipy.special import digamma
import matplotlib.pyplot as plt
from pyLogVar import LogVar
import json
import itertools
from model import Pedigree
from PedigreeHypergraph import PedigreeHG
from AutosomalDistribution import *
from HHMMUpDown import HiddenMarkovModelMessagePasser,MessagePassingHG
from MMUpDown import MarkovModelMessagePasser

def printListOfNodes(nodes):
    print('[ '),
    for node in nodes:
        print(str(node._id)+' '),
    print(']')

def pedigreeExample():
    pedigreeFolderName = '/Users/Eddie/kec-bot/app/pedigreeDataOLDBUTWORKS/'


    # 3239PB, 3818J, 5092AD, 5546EL, 5596IN, 5712CS, 5713BS, 5865MH, 5992VM, 5992VM, 6050MM
    weird = ['5777AH','235TL','3239PB','5092AD','5546EL',\
    '5596IN','5712CS','5713BS','5865MH','5992VM','5992VM','6050MM']

    # pedigreeNames = []
    # for filename in os.listdir(pedigreeFolderName):
    #     if('json' in filename and 'test' not in filename):
    #         if(filename.strip('.json') in weird):
    #             continue
    #         pedigreeNames.append(filename.strip('.json'))

    pedigreeNames = ['3818J']

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

    hg = allPedigrees['3818J']

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
    hg = MessagePassingHG(2)
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


# hg,msg = pedigreeExample()
# hg,msg = cycleExample4(False)
# hg,msg = cycleExample5(False)
hg,msg = cycleExample5_1(False)
# hg,msg = cycleExample6()
# hg,msg = cycleExample7(False)
# hg,msg = disjointExample()
# hg,msg = pedigreeExample()

# print('Feedback set: '+str(msg._feedbackSet))
# print('Feedback set: '+str(msg._feedbackStack))

# n4 = hg.addNode(4)
# val0,str0 = msg.fullBruteForce([n4],[0])
# val1,str1 = msg.fullBruteForce([n4],[1])
# print(val0+val1)
# print(str0)
# assert 0
def HMMTest():
    start = time.time()
    msg.preprocess()
    end = time.time()
    print('Preprocess time: '+str(end-start))
    msg.aTest()
    start = time.time()
    msg.getStats()
    end = time.time()
    print('Traversal time: '+str(end-start))
    msg.aTest()

def MMTest():

    for node in msg.nodes:
        for i in range(node.N):
            comp = msg.getW(node,i)
            bf = msg.getWBruteForce([node],[i])
            print(comp._logVal - bf._logVal)

MMTest()