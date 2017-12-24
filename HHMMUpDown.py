from HypergraphBase import NodeBase,EdgeBase,BaseHyperGraph
# from CythonCode.HypergraphBaseFast import NodeBase,EdgeBase,BaseHyperGraph
# from CythonCode.LogVarCode import LogVar
from pyLogVar import LogVar
from cycleDetector import identifyCycles
import numpy as np
import graphviz
import itertools
from functools import reduce

def addIt(x,y):
    return x+y
def mulIt(x,y):
    return x*y

def calcN(person):
    return 2

def getY(person):
    if('diagnoses' in dir(person)):
        return int(len(person.diagnoses) > 0)
    return 0

def prettyPrint(d, indent=0):
    for key, value in d.items():
        print('\t' * indent + str(key))
        if isinstance(value, dict):
            prettyPrint(value, indent+1)
        else:
            print('\t' * (indent+1) + str(value))

class NodeForHMM(NodeBase):
    def __init__(self,y,N=None):
        super(NodeForHMM,self).__init__()
        if(N == None):
            self.N = calcN(self)
        else:
            self.N = N

        self.y = y

        self.inFBS = False

        self._UDeps = set([self])
        self._VDeps = set([self])
        self._aDeps = set([self])
        self._bDeps = set([self])

        self.reset()

    def reset(self):
        self._U = {}
        self._V = {}
        self._a = {}
        self._b = {}
        self._fullJoint = {}

    def setMsg(self,msg):
        self._msg = msg
        self._pi = msg._pi
        self._L = msg._L
        self._trans = msg._trans

    def _getN(self,node,conditioning):
        if(node in conditioning):
            return (conditioning[node],)
        return range(node.N)

    def keyFromCond(self,conditioning):
        return tuple(conditioning.items())

    """ ------------------------------------------------------------------------------------ """

    def aKey(self,conditioning):
        actualCond = {k:v for k,v in conditioning.items() if k in self._aDeps}
        key = self.keyFromCond(actualCond)
        return key

    def needToComputeA(self,edge,i,key):
        return edge not in self._a or i not in self._a[edge] or key not in self._a[edge][i]

    def setAVal(self,edge,i,key,aVal):
        if(edge not in self._a): self._a[edge] = {}
        if(i not in self._a[edge]): self._a[edge][i] = {}
        assert key not in self._a[edge][i]
        self._a[edge][i][key] = aVal

    def getAVal(self,edge,i,key):
        try:
            return LogVar(self._a[edge][i][key])
        except:
            print('YOU DONE FUCKED UP')
            prettyPrint(self._a)
            self._a[edge]
            self._a[edge][i]
            self._a[edge][i][key]
            assert 0

    def getA(self,edge,i,conditioning):

        key = self.aKey(conditioning)

        """ a_n_e(i) = P(Y\!(e,n),n_x=i) """
        if(self.needToComputeA(edge,i,key)):

            # if we don't have the right conditioning, add it
            aVal = self._computeA(edge,i,conditioning)
            self.setAVal(edge,i,key,aVal)

            if(self._msg.preprocessing):
                self._msg._traverseOrder.append(('a',self,edge,i,conditioning))
        else:
            aVal = self.getAVal(edge,i,key)
        return aVal

    def _computeA(self,edge,i,conditioning):

        a_ = LogVar(1)

        """ All nodes up from this node """
        a_ *= self.getU(i,conditioning)
        self._aDeps |= self._UDeps

        """ All nodes down from this node but not down edge """
        if(len(self._downEdges) > 1):

            for e in self._downEdges:
                if(e == edge): continue

                """ Down each mate branch """
                a_ *= self.getV(e,i,conditioning)
            self._aDeps |= self._VDeps

        return a_

    def getMarginalizedA(self,edge,i,nodesToKeep,feedbackSet):

        margOut = filter(lambda n:n not in nodesToKeep and n in self._aDeps, feedbackSet)
        a = LogVar(0)
        for X in itertools.product(*[range(n.N) for n in margOut]):
            conditioning = {node:_i for node,_i in zip(margOut,X)}
            conditioning.update(nodesToKeep)
            a += self.getA(edge,i,conditioning)
        return a

    """ ------------------------------------------------------------------------------------ """

    def bKey(self,conditioning):
        actualCond = {k:v for k,v in conditioning.items() if k in self._bDeps}
        key = self.keyFromCond(actualCond)
        return key

    def needToComputeB(self,X,key):
        return X not in self._b or key not in self._b[X]

    def setBVal(self,X,key,bVal):
        if(X not in self._b): self._b[X] = {}
        self._b[X][key] = bVal

    def getBVal(self,X,key):
        try:
            return LogVar(self._b[X][key])
        except:
            print('YOU DONE FUCKED UP')
            prettyPrint(self._b)
            self._b[X][key]
            assert 0

    def getB(self,X,conditioning):

        key = self.bKey(conditioning)

        """ B_n(X) = P(n_y,Y\^(n)_y|^(n)_x=X) """
        if(self.needToComputeB(X,key)):

            bVal = self._computeB(X,conditioning)
            self.setBVal(X,key,bVal)

            if(self._msg.preprocessing):
                self._msg._traverseOrder.append(('b',self,X,conditioning))
        else:
            bVal = self.getBVal(X,key)
        return bVal

    def _computeB(self,X,conditioning):

        b_ = LogVar(0)
        if(len(self._downEdges) > 0):

            for k in self._getN(self,conditioning):

                """ Prob of this sibling """
                _prod = LogVar(self._trans(self._parents,self,X,k)) * self._L(self,k)

                """ Branch down from sibling """
                for e in self._downEdges:
                    _prod *= self.getV(e,k,conditioning)

                b_ += _prod

            self._bDeps |= self._VDeps
        else:

            for k in self._getN(self,conditioning):

                """ Prob of this sibling """
                b_ += self._trans(self._parents,self,X,k) * self._L(self,k)

        return b_


    def getMarginalizedB(self,X,nodesToKeep,feedbackSet):

        margOut = filter(lambda n:n not in nodesToKeep and n in self._bDeps, feedbackSet)
        b = LogVar(0)
        for _X in itertools.product(*[range(n.N) for n in margOut]):
            conditioning = {node:i for node,i in zip(margOut,_X)}
            conditioning.update(nodesToKeep)
            b += self.getB(X,conditioning)
        return b

    """ ------------------------------------------------------------------------------------ """

    def UKey(self,conditioning):
        actualCond = {k:v for k,v in conditioning.items() if k in self._UDeps}
        key = self.keyFromCond(actualCond)
        return key

    def needToComputeU(self,i,key):
        return i not in self._U or key not in self._U[i]

    def setUVal(self,i,key,uVal):
        if(i not in self._U): self._U[i] = {}
        self._U[i][key] = uVal

    def getUVal(self,i,key):
        return LogVar(self._U[i][key])

    def getU(self,i,conditioning):

        key = self.UKey(conditioning)

        """ U_n(i) = P(n_y,^(n)_y,n_x=i) """
        if(self.needToComputeU(i,key)):

            uVal = self._computeU(i,conditioning)
            self.setUVal(i,key,uVal)

            if(self._msg.preprocessing):
                self._msg._traverseOrder.append(('U',self,i,conditioning))
        else:
            uVal = self.getUVal(i,key)
        return uVal

    def _computeU(self,i,conditioning):

        parents = self._parents
        if(len(parents) == 0 or self.inFBS):

            if(self.inFBS):
                # Prob if we conditioned on this
                u = LogVar(int(i==conditioning[self]))
            else:
                # Root dist prob
                u = LogVar(self._pi(self,i)) * self._L(self,i)
        else:

            u = LogVar(0)
            for X in itertools.product(*[self._getN(p,conditioning) for p in parents]):

                """ Prob of this node """
                prod = LogVar(self._trans(parents,self,X,i))

                """ Branch out from each parent """
                for parent,j in zip(parents,X):
                    prod *= parent.getA(self._upEdge,j,conditioning)
                    self._UDeps |= parent._aDeps

                """ Branch out from each sibling """
                for sibling in self._upEdge._children:
                    if(sibling==self):continue

                    prod *= sibling.getB(X,conditioning)
                    self._UDeps |= sibling._bDeps

                u += prod

            u *= self._L(self,i)

        return u

    def getMarginalizedU(self,i,nodesToKeep,feedbackSet):

        margOut = filter(lambda n:n not in nodesToKeep and n in self._UDeps, feedbackSet)
        u = LogVar(0)
        for X in itertools.product(*[range(n.N) for n in margOut]):
            conditioning = {node:i for node,i in zip(margOut,X)}
            conditioning.update(nodesToKeep)
            u += self.getU(i,conditioning)
        return u

    """ ------------------------------------------------------------------------------------ """

    def VKey(self,conditioning):
        actualCond = {k:v for k,v in conditioning.items() if k in self._VDeps}
        key = self.keyFromCond(actualCond)
        return key

    def needToComputeV(self,edge,i,key):
        return edge not in self._V or i not in self._V[edge] or key not in self._V[edge][i]

    def setVVal(self,edge,i,key,vVal):
        if(edge not in self._V): self._V[edge] = {}
        if(i not in self._V[edge]): self._V[edge][i] = {}
        assert key not in self._V[edge][i]
        self._V[edge][i][key] = vVal

    def getVVal(self,edge,i,key):
        return LogVar(self._V[edge][i][key])

    def getV(self,edge,i,conditioning):

        key = self.VKey(conditioning)

        """ V_n_e(i) = P(!(e,n)|n_x=i) """
        if(self.needToComputeV(edge,i,key)):

            # if we don't have the right conditioning, add it
            vVal = self._computeV(edge,i,conditioning)

            self.setVVal(edge,i,key,vVal)
            if(self._msg.preprocessing):

                self._msg._traverseOrder.append(('V',self,edge,i,conditioning))
        else:
            vVal = self.getVVal(edge,i,key)
        return vVal

    def _computeV(self,edge,i,conditioning):

        if(self.inFBS or len(self._downEdges) == 0):
            return LogVar(1)

        mates = [x for x in edge._parents if x != self]

        v = LogVar(0)
        for X_,X in zip(itertools.product(*[self._getN(m,conditioning) for m in mates]),\
            itertools.product(*[self._getN(m,conditioning) if m!=self else [i] for m in edge._parents])):

            prod = LogVar(1)

            """ Branch out from each mate """
            for mate,j in zip(mates,X_):

                prod *= mate.getA(edge,j,conditioning)
                self._VDeps |= mate._aDeps

            """ Branch out from each child """
            for child in edge._children:
                prod *= child.getB(X,conditioning)
                self._VDeps |= child._bDeps

            v += prod

        return v

    def getMarginalizedV(self,edge,i,nodesToKeep,feedbackSet):

        margOut = filter(lambda n:n not in nodesToKeep and n in self._VDeps, feedbackSet)
        v = LogVar(0)
        for X in itertools.product(*[range(n.N) for n in margOut]):
            conditioning = {node:i for node,i in zip(margOut,X)}
            conditioning.update(nodesToKeep)
            v += self.getV(edge,i,conditioning)
        return v

    """ ------------------------------------------------------------------------------------ """

    def sortaRootProb(self,conditioning):
        return self._msg.sortaRootProb(conditioning)

    def accumulateFullJoint(self,feedbackSet):

        if(len(self._downEdges) > 0):

            for i in range(self.N):

                total = LogVar(0)
                for X in itertools.product(*[range(n.N) for n in feedbackSet]):
                    conditioning = {node:_i for node,_i in zip(feedbackSet,X)}

                    prod = self.getU(i,conditioning)

                    for edge in self._downEdges:
                        prod *= self.getV(edge,i,conditioning)

                    prod *= self.sortaRootProb(conditioning)

                    total += prod

                self._fullJoint[i] = total
        else:
            for i in range(self.N):

                total = LogVar(0)
                for X in itertools.product(*[range(n.N) for n in feedbackSet]):
                    conditioning = {node:_i for node,_i in zip(feedbackSet,X)}

                    total += self.getU(i,conditioning) * self.sortaRootProb(conditioning)

                self._fullJoint[i] = total

    def getFullJoint(self,i):
        return self._fullJoint[i]

    def updateFullJoint(self,i,val):
        if(i not in self._fullJoint): self._fullJoint[i] = LogVar(0)
        self._fullJoint[i] += val

class MessagePassingHG(BaseHyperGraph):
    def __init__(self,N=None):
        self.N = N
        super(MessagePassingHG,self).__init__()
        self.setNodeType(NodeForHMM)

    def addNode(self,ID):
        y = getY(ID)
        N = self.N
        return super(MessagePassingHG,self).addNode(ID,y,N)

    def draw(self):

        assert self._initialized, 'call the function \'hypergraph.initialize()\''

        """ Draws the hypergraph using graphviz """
        d = graphviz.Digraph()
        for e in self._edges:
            eId = e._id
            for p in e._parents:
                pId = p._id
                d.edge('n('+str(pId)+')','E('+str(eId)+')')
            for c in e._children:
                cId = c._id
                d.edge('E('+str(eId)+')','n('+str(cId)+')')
        d.render()
        return d

class HiddenMarkovModelMessagePasser():

    def __init__(self,hg,hgParamFunction):

        self._hyperGraph = hg
        self._paramGenerator = hgParamFunction

        self.nodes = hg._nodes
        self.edges = hg._edges

        self._params = self._paramGenerator(self._hyperGraph)
        self._trans = self._params['transDist']
        self._L = self._params['emissionDist']

        self._conditioning = {}

        self._pi = self._params['rootDist']
        self._srp = {}

    def preprocess(self):

        # for preprocessing, need to find a fvs, sortaRootProb,
        # and find the path that the algorithm travels so that
        # we can just use that

        feedbackSet,blockManager = identifyCycles(self._hyperGraph)

        blockManager.accumulateAllChains()

        self._blockManager = blockManager
        self._feedbackSet = sorted(feedbackSet,key=lambda x:x._id)

        # self._feedbackSet = [self._hyperGraph.getNode(4),self._hyperGraph.getNode(5)]

        for node in self.nodes: node.setMsg(self)

        self._sortaRootDeps = [node for node in self._feedbackSet if len(node._parents) == 0 or \
                            (len([x for x in node._parents if x not in self._feedbackSet]) == 0 and \
                             len([x for x in node._upEdge._children if x not in self._feedbackSet]) == 0)]

        for node in self.nodes:
            if node in self._feedbackSet:
                node.inFBS = True

        self._traverseOrder = []

        self.preprocessing = True
        self.computeForPreprocessing()
        self.preprocessing = False

    """ -------------------------------------------------------------------------------------- """

    def sortaRootProb(self,conditioning):

        key = tuple(conditioning.items())
        if(key in self._srp):
            return self._srp[key]

        prod = LogVar(1)
        if len(self._sortaRootDeps) > 0:

            for sortaRoot in self._sortaRootDeps:
                j = conditioning[sortaRoot]
                _prod = LogVar(self._L(sortaRoot,j))

                if(len(sortaRoot._parents) == 0):
                    _prod *= self._pi(sortaRoot,j)
                else:
                    parents = sortaRoot._parents
                    X_ = [conditioning[p] for p in parents]
                    _prod *= self._trans(parents,sortaRoot,X_,j)
                prod *= _prod

        self._srp[key] = prod
        return prod

    """ -------------------------------------------------------------------------------------- """

    def computeForPreprocessing(self):

        for node in self.nodes:
            if(node in self._feedbackSet): continue

            for X in itertools.product(*[range(n.N) for n in self._feedbackSet]):
                conditioning = {node:i for node,i in zip(self._feedbackSet,X)}
                for i in range(node.N):
                    node.getU(i,conditioning)
                for i in range(node.N):
                    for edge in node._downEdges:
                        node.getV(edge,i,conditioning)

        for node in self.nodes:
            if(node in self._feedbackSet): continue

            node.accumulateFullJoint(self._feedbackSet)

        # compute the probs for the fbs
        aLeaf = self._hyperGraph._leaves.__iter__().__next__()
        self._srp = {}

        for X in itertools.product(*[range(n.N) for n in self._feedbackSet]):

            conditioning = {node:i for node,i in zip(self._feedbackSet,X)}

            for i in range(aLeaf.N):
                val = aLeaf.getU(i,conditioning)
                sr = self.sortaRootProb(conditioning)

                for node,x in zip(self._feedbackSet,X):
                    node.updateFullJoint(x,val*sr)

        for node in self._feedbackSet:
            total = LogVar(0)
            for i in range(node.N):
                total += node._fullJoint[i]

    """ -------------------------------------------------------------------------------------- """

    def getCounts(self):

        for node in self.nodes:
            node.reset()

        for thing in self._traverseOrder:

            theType = thing[0]
            args = thing[1:]

            if(theType == 'a'):
                node,edge,i,conditioning = args
                node.getA(edge,i,conditioning)

            elif(theType == 'b'):
                node,X,conditioning = args
                node.getB(X,conditioning)

            elif(theType == 'U'):
                node,i,conditioning = args
                node.getU(i,conditioning)

            elif(theType == 'V'):
                node,edge,i,conditioning = args
                node.getV(edge,i,conditioning)

        for node in self.nodes:
            if(node in self._feedbackSet): continue

            node.accumulateFullJoint(self._feedbackSet)

        # compute the probs for the fbs
        aLeaf = self._hyperGraph._leaves.__iter__().__next__()
        self._srp = {}

        for X in itertools.product(*[range(n.N) for n in self._feedbackSet]):

            conditioning = {node:i for node,i in zip(self._feedbackSet,X)}

            for i in range(aLeaf.N):
                val = aLeaf.getU(i,conditioning)
                sr = self.sortaRootProb(conditioning)
                [node.updateFullJoint(x,val*sr) for node,x in zip(self._feedbackSet,X)]

        for node in self._feedbackSet:
            total = reduce(addIt,node._fullJoint)

    def isolatedParentJoint(self,node,X,i,totalProb=None):

        # compute the probs for the fbs
        aLeaf = list(self._hyperGraph._leaves)[0]

        parents = node._parents
        jointProb = LogVar(0)

        parentCond = {p:j for p,j in zip(parents,X)}

        margOut = [n for n in self._feedbackSet if n!=node and n not in parents]

        for _X in itertools.product(*[range(n.N) for n in margOut]):

            conditioning = {_node:_i for _node,_i in zip(margOut,_X)}
            conditioning.update(parentCond)
            conditioning[node] = i

            for j in range(aLeaf.N):
                val = aLeaf.getU(j,conditioning)
                sr = self.sortaRootProb(conditioning)
                jointProb += val*sr

        if(totalProb):
            jointProb /= totalProb
        return jointProb

    def probOfParentsProducingNode(self,node,X,i,totalProb=None):

        if(node in self._sortaRootDeps):
            return self.isolatedParentJoint(node,X,i,totalProb)

        parents = node._parents

        # make sure we sum over all of the possible nodes in the fbs
        latentRanges = [[X[parents.index(n)]] if n in parents else range(n.N) for n in self._feedbackSet]

        if(node.inFBS):
            latentRanges[self._feedbackSet.index(node)] = [i]

        """ Prob of this node """
        prob = LogVar(self._trans(parents,node,X,i)) * self._L(node,i)

        total = LogVar(0)

        # print('\ni is: '+str(i))
        # print('node is: '+str(node))
        # print('X is: '+str(X))
        # print('parents is: '+str(parents))
        # print('self._feedbackSet is: '+str(self._feedbackSet))
        # print('latentRanges: '+str(latentRanges))
        for S in itertools.product(*latentRanges):

            # print('S IS '+str(S))
            familyInFBS = {n:s for n,s in zip(self._feedbackSet,S)}

            """ Down this node """
            _prob = LogVar(1)
            for e in node._downEdges:
                v = node.getMarginalizedV(e,i,familyInFBS,self._feedbackSet)
                # print('e: '+str(e)+' v: '+str(v))
                _prob *= v

            """ Out from each sibling """
            for sibling in node._upEdge._children:
                if(sibling is node): continue
                b = sibling.getMarginalizedB(X,familyInFBS,self._feedbackSet)
                # print('sibling: '+str(sibling)+' b: '+str(b))
                _prob *= b

            """ Out from each parent """
            for parent,j in zip(parents,X):
                a = parent.getMarginalizedA(node._upEdge,j,familyInFBS,self._feedbackSet)
                # print('parent: '+str(parent)+' edge: '+str(node._upEdge)+' j: '+str(j)+' a: '+str(a))
                _prob *= a

            srp = self.sortaRootProb(familyInFBS)
            # print('srp: '+str(srp))
            _prob *= srp
            total += _prob

        prob *= total

        """ Normalize """
        if(totalProb):
            prob /= totalProb

        return prob

    def getStats(self):

        self.getCounts()

        for node in self.nodes:
            for i in range(node.N):
                uVal = node.getFullJoint(i)

        for n in self.nodes:
            for i in range(n.N):
                if(len(n._parents) > 0):
                    for X in itertools.product(*[range(p.N) for p in n._parents]):
                        genProb = self.probOfParentsProducingNode(n,X,i)

    """ -------------------------------------------------------------------------------------- """

    def probOfAllNodeObservations(self):

        """ P(Y) = sum_i(U_l_e(i)) for any leaf l """
        aLeaf = list(self._hyperGraph._leaves)[0]

        total = LogVar(0)
        for i in range(aLeaf.N):
            _u = aLeaf.getFullJoint(i)
            total += _u
        return total

    def aTest(self,printStuff=False):

        if(self.preprocessing):
            self.computeForPreprocessing()
        else:
            self.getStats()

        if(printStuff):
            for node in self.nodes:
                print('\n----------------\nNode: '+str(node))
                print('_a:\n')
                prettyPrint(node._a)
                print('_b:\n')
                prettyPrint(node._b)
                print('_U:\n')
                prettyPrint(node._U)
                print('_V:\n')
                prettyPrint(node._V)
                print('U deps: '+str(node._UDeps))
                print('V deps: '+str(node._VDeps))
                print('a deps: '+str(node._aDeps))
                print('b deps: '+str(node._bDeps))

        correct = self.probOfAllNodeObservations()
        if(printStuff):
            print('\n\nP(Y) is: '+str(correct))
            print('Feedback set: '+str(self._feedbackSet))

        for node in self.nodes:
            total = LogVar(0)
            for i in range(node.N):

                uVal = node.getFullJoint(i)
                uCopy = LogVar(uVal)
                uCopy /= correct
                total += uCopy

        def isclose(a, b, rel_tol=1e-05, abs_tol=0.0):
            return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

        """ Make sure that P(Y) can be computed everywhere correctly """
        for node in self.nodes:
            total = LogVar(0)
            for i in range(node.N):

                uVal = node.getFullJoint(i)
                uCopy = LogVar(uVal)
                total += uCopy

            if(not isclose(total.logVal,correct.logVal)):
                if(printStuff):
                    print('\nFailed the P(Y) test!')
                    print('Node: '+str(node))
                    print('total: '+str(float(total)))
                    print('correct: '+str(float(correct)))
                print('up edge: '+str(node._upEdge))
                print('down edges: '+str(node._downEdges))
                assert 0
            else:
                if(printStuff):
                    print('Passed P(Y) test for node '+str(node)+' with prob '+str(total))
        if(printStuff):
            print('Passed the P(Y) tests!')

        """ Make sure that the statistics we care about sum to 1 """
        for n in self.nodes:

            total1 = LogVar(0)
            total2 = LogVar(0)
            for i in range(n.N):

                prob = n.getFullJoint(i)
                prob /= correct
                total1 += prob

                if(len(n._parents) > 0):

                    shouldEqualProb = LogVar(0)
                    for X in itertools.product(*[range(p.N) for p in n._parents]):

                        prob2 = self.probOfParentsProducingNode(n,X,i,correct)
                        total2 += prob2
                        shouldEqualProb += prob2

                    if(not isclose(float(shouldEqualProb),float(prob))):
                        if(printStuff):
                            print('\n\n\n\nFailed this test for node: '+str(n))
                            print('This is wrong: '+str(float(shouldEqualProb*correct)))
                            print('This is right: '+str(float(prob*correct)))
                        print('up edge: '+str(n._upEdge))
                        print('down edges: '+str(n._downEdges))
                        assert 0
                    else:
                        if(printStuff):
                            print('Passed the marg out test P(x_c|y)=sumP(x_c,{x_p}|y) for node '+str(n)+' with prob '+str(prob))
                            print('==========================================================')

            if(not isclose(float(total1),1.)):
                if(printStuff):
                    print('The sum of all P(x_i|Y) doesn\'t equal 1 for node: '+str(n)+', total1: '+str(total1.logVal))
                print('up edge: '+str(node._upEdge))
                print('down edges: '+str(node._downEdges))
                assert 0

            if(len(n._parents) > 0):
                if(not isclose(float(total2),1.)):
                    if(printStuff):
                        'The sum of all P({x_p}|Y) isn\'t 1 for node: '+str(n)+', total2: '+str(total2.logVal)
                    print('up edge: '+str(node._upEdge))
                    print('down edges: '+str(node._downEdges))
                    assert 0





    def printTraversalPath(self):
        pass
