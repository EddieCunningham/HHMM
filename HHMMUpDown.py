from HypergraphBase import *
from pyLogVar import *
from cycleDetector import *
import numpy as np

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

class NodeForHMM(Node):
    def __init__(self,y,N=None):
        super(NodeForHMM,self).__init__()
        if(N == None):
            self.N = calcN(self)
        else:
            self.N = N

        self.y = y

        self.inFBS = False
        self.gaveRootProb = [False for _ in range(self.N)]

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
        return tuple(sorted(conditioning.items()))

    """ ------------------------------------------------------------------------------------ """

    def aKey(self,conditioning):
        actualCond = {k:v for k,v in conditioning.items() if k in self._aDeps}
        key = self.keyFromCond(actualCond)
        return key

    def needToComputeA(self,edge,i,conditioning):
        key = self.aKey(conditioning)
        if(edge not in self._a or i not in self._a[edge] or key not in self._a[edge][i]):
            return True
        return False

    def setAVal(self,edge,i,conditioning,aVal):
        key = self.aKey(conditioning)
        if(edge not in self._a): self._a[edge] = {}
        if(i not in self._a[edge]): self._a[edge][i] = {}
        assert key not in self._a[edge][i]
        self._a[edge][i][key] = aVal

    def getAVal(self,edge,i,conditioning):
        key = self.aKey(conditioning)
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
        # print('a node: '+str(self)+' edge: '+str(edge)+' i: '+str(i)+' cond: '+str(conditioning))

        """ a_n_e(i) = P(Y\!(e,n),n_x=i) """
        if(self.needToComputeA(edge,i,conditioning)):

            # if we don't have the right conditioning, add it
            aVal = self._computeA(edge,i,conditioning)
            self.setAVal(edge,i,conditioning,aVal)

            if(self._msg.preprocessing):
                self._msg._traverseOrder.append(('a',self,edge,i,conditioning))
        else:
            aVal = self.getAVal(edge,i,conditioning)
        return aVal

    def _computeA(self,edge,i,conditioning):

        a_ = LogVar(1)

        """ All nodes up from this node """
        a_ *= self.getU(i,conditioning)
        self._aDeps |= self._UDeps

        """ All nodes down from this node but not down edge """
        for e in [x for x in self._downEdges if x!=edge]:
            a_ *= self.getV(e,i,conditioning)
            self._aDeps |= self._VDeps

        return a_

    def getMarginalizedA(self,edge,i,nodesToKeep,feedbackSet):
        # marginalize out the condition nodes that aren't
        # in nodesToKeep
        assert len([n for n in nodesToKeep if n not in feedbackSet]) == 0

        margOut = [n for n in feedbackSet if n not in nodesToKeep and n in self._aDeps]

        total = LogVar(0)
        for X in itertools.product(*[range(n.N) for n in margOut]):

            conditioning = {node:i for node,i in zip(margOut,X)}
            conditioning.update(nodesToKeep)

            total += self.getA(edge,i,conditioning)
            # total += self.getAVal(edge,i,conditioning)

        return total

    """ ------------------------------------------------------------------------------------ """

    def bKey(self,conditioning):
        actualCond = {k:v for k,v in conditioning.items() if k in self._bDeps}
        key = self.keyFromCond(actualCond)
        return key

    def needToComputeB(self,X,conditioning):
        key = self.bKey(conditioning)
        if(X not in self._b or key not in self._b[X]):
            return True
        return False

    def setBVal(self,X,conditioning,bVal):
        key = self.bKey(conditioning)
        if(X not in self._b): self._b[X] = {}
        self._b[X][key] = bVal

    def getBVal(self,X,conditioning):
        key = self.bKey(conditioning)
        try:
            return LogVar(self._b[X][key])
        except:
            print('YOU DONE FUCKED UP')
            prettyPrint(self._b)
            self._b[X][key]
            assert 0

    def getB(self,X,conditioning):
        # print('b node: '+str(self)+' X: '+str(X)+' cond: '+str(conditioning))

        """ B_n(X) = P(n_y,Y\^(n)_y|^(n)_x=X) """
        if(self.needToComputeB(X,conditioning)):

            bVal = self._computeB(X,conditioning)
            self.setBVal(X,conditioning,bVal)

            if(self._msg.preprocessing):
                self._msg._traverseOrder.append(('b',self,X,conditioning))
        else:
            bVal = self.getBVal(X,conditioning)
        return bVal

    def _computeB(self,X,conditioning):
        b_ = LogVar(0)
        for k in self._getN(self,conditioning):

            _prod = LogVar(1)

            """ Prob of this sibling """
            _prod *= self._trans(sorted(self._parents),self,X,k)
            _prod *= self._L(self,k)

            """ Branch down from sibling """
            for e in self._downEdges:
                _prod *= self.getV(e,k,conditioning)
                self._bDeps |= self._VDeps

            b_ += _prod
        return b_

    def getMarginalizedB(self,X,nodesToKeep,feedbackSet):
        # marginalize out the condition nodes that aren't
        # in nodesToKeep
        assert len([n for n in nodesToKeep if n not in feedbackSet]) == 0

        margOut = [n for n in feedbackSet if n not in nodesToKeep and n in self._bDeps]

        total = LogVar(0)
        for _X in itertools.product(*[range(n.N) for n in margOut]):

            conditioning = {node:i for node,i in zip(margOut,_X)}
            conditioning.update(nodesToKeep)

            total += self.getB(X,conditioning)
            # total += self.getBVal(X,conditioning)

        return total

    """ ------------------------------------------------------------------------------------ """

    def UKey(self,conditioning):
        actualCond = {k:v for k,v in conditioning.items() if k in self._UDeps}
        key = self.keyFromCond(actualCond)
        return key

    def needToComputeU(self,i,conditioning):
        key = self.UKey(conditioning)
        if(i not in self._U or key not in self._U[i]):
            return True
        return False

    def setUVal(self,i,conditioning,uVal):
        key = self.UKey(conditioning)
        if(i not in self._U): self._U[i] = {}
        self._U[i][key] = uVal

    def getUVal(self,i,conditioning):
        key = self.UKey(conditioning)
        return LogVar(self._U[i][key])

    def getU(self,i,conditioning):
        # print('U node: '+str(self)+' i: '+str(i)+' cond: '+str(conditioning))

        """ U_n(i) = P(n_y,^(n)_y,n_x=i) """
        if(self.needToComputeU(i,conditioning)):

            uVal = self._computeU(i,conditioning)
            self.setUVal(i,conditioning,uVal)

            if(self._msg.preprocessing):
                self._msg._traverseOrder.append(('U',self,i,conditioning))
        else:
            uVal = self.getUVal(i,conditioning)
        return uVal

    def _computeU(self,i,conditioning):

        u = LogVar(0)

        parents = sorted(self._parents)
        if(len(parents) == 0 or self.inFBS):

            if(self.inFBS):
                # Prob if we conditioned on this
                rootProb = LogVar(int(i==conditioning[self]))
                emissionProb = LogVar(1)
            else:
                # Root dist prob
                rootProb = LogVar(self._pi(self,i))
                emissionProb = LogVar(self._L(self,i))

            u = rootProb*emissionProb

        else:

            for X in itertools.product(*[self._getN(p,conditioning) for p in parents]):

                prod = LogVar(1)

                """ Prob of this node """
                prod *= LogVar(self._trans(parents,self,X,i))
                prod *= LogVar(self._L(self,i))

                """ Branch out from each parent """
                for parent,j in zip(parents,X):

                    prod *= parent.getA(self._upEdge,j,conditioning)
                    self._UDeps |= parent._aDeps

                """ Branch out from each sibling """
                for sibling in [x for x in self._upEdge._children if x!=self]:

                    prod *= sibling.getB(X,conditioning)
                    self._UDeps |= sibling._bDeps

                u += prod
        return u

    def getMarginalizedU(self,i,nodesToKeep,feedbackSet):
        # marginalize out the condition nodes that aren't
        # in nodesToKeep
        assert len([n for n in nodesToKeep if n not in feedbackSet]) == 0

        margOut = [n for n in feedbackSet if n not in nodesToKeep and n in self._UDeps]

        total = LogVar(0)
        for X in itertools.product(*[range(n.N) for n in margOut]):

            conditioning = {node:_i for node,_i in zip(margOut,X)}
            conditioning.update(nodesToKeep)

            total += self.getU(i,conditioning)
            # total += self.getUVal(i,conditioning)

        return total

    """ ------------------------------------------------------------------------------------ """

    def VKey(self,conditioning):
        actualCond = {k:v for k,v in conditioning.items() if k in self._VDeps}
        key = self.keyFromCond(actualCond)
        return key

    def needToComputeV(self,edge,i,conditioning):
        key = self.VKey(conditioning)
        if(edge not in self._V or i not in self._V[edge] or key not in self._V[edge][i]):
            return True
        return False

    def setVVal(self,edge,i,conditioning,vVal):
        key = self.VKey(conditioning)
        if(edge not in self._V): self._V[edge] = {}
        if(i not in self._V[edge]): self._V[edge][i] = {}
        assert key not in self._V[edge][i]
        self._V[edge][i][key] = vVal

    def getVVal(self,edge,i,conditioning):
        key = self.VKey(conditioning)
        return LogVar(self._V[edge][i][key])

    def getV(self,edge,i,conditioning):
        # print('V node: '+str(self)+' edge: '+str(edge)+' i: '+str(i)+' cond: '+str(conditioning))

        """ V_n_e(i) = P(!(e,n)|n_x=i) """
        if(self.needToComputeV(edge,i,conditioning)):

            # if we don't have the right conditioning, add it
            vVal = self._computeV(edge,i,conditioning)

            self.setVVal(edge,i,conditioning,vVal)
            if(self._msg.preprocessing):

                self._msg._traverseOrder.append(('V',self,edge,i,conditioning))
        else:
            vVal = self.getVVal(edge,i,conditioning)
        return vVal

    def _computeV(self,edge,i,conditioning):

        if(self.inFBS or len(self._downEdges) == 0):
            return LogVar(1)

        v = LogVar(0)

        assert self in edge._parents

        mates = [x for x in sorted(edge._parents) if x != self]
        for X_,X in zip(itertools.product(*[self._getN(m,conditioning) for m in mates]),\
            itertools.product(*[self._getN(m,conditioning) if m!=self else [i] for m in sorted(edge._parents)])):

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
        # marginalize out the condition nodes that aren't
        # in nodesToKeep
        assert len([n for n in nodesToKeep if n not in feedbackSet]) == 0, [n for n in nodesToKeep if n not in feedbackSet]

        margOut = [n for n in feedbackSet if n not in nodesToKeep and n in self._VDeps]

        total = LogVar(0)
        for X in itertools.product(*[range(n.N) for n in margOut]):

            conditioning = {node:i for node,i in zip(margOut,X)}
            conditioning.update(nodesToKeep)

            total += self.getV(edge,i,conditioning)
            # total += self.getVVal(edge,i,conditioning)

        return total

    """ ------------------------------------------------------------------------------------ """

    def sortaRootProb(self,conditioning):
        return self._msg.sortaRootProb(conditioning)

    def accumulateFullJoint(self,feedbackSet):

        for i in range(self.N):

            for X in itertools.product(*[range(n.N) for n in feedbackSet]):
                conditioning = {node:i for node,i in zip(feedbackSet,X)}

                prod = LogVar(1)

                u = self.getUVal(i,conditioning)
                prod *= u

                for edge in self._downEdges:
                    v = self.getVVal(edge,i,conditioning)
                    prod *= v

                sr = self.sortaRootProb(conditioning)
                prod *= sr

                self.updateFullJoint(i,prod)

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

    def preprocess(self):

        # for preprocessing, need to find a fvs, sortaRootProb,
        # and find the path that the algorithm travels so that
        # we can just use that

        feedbackSet,blockManager = identifyCycles(self._hyperGraph)

        # feedbackSet = [hg.getNode(1),hg.getNode(2),hg.getNode(3),hg.getNode(4),hg.getNode(5)]

        blockManager.accumulateAllChains()

        self._blockManager = blockManager
        self._feedbackSet = feedbackSet

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

    def getStats(self):

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
        aLeaf = self._hyperGraph._leaves.__iter__().next()

        for X in itertools.product(*[range(n.N) for n in self._feedbackSet]):

            conditioning = {node:i for node,i in zip(self._feedbackSet,X)}

            for i in range(aLeaf.N):
                val = aLeaf.getUVal(i,conditioning)
                sr = self.sortaRootProb(conditioning)

                for node,x in zip(self._feedbackSet,X):
                    node.updateFullJoint(x,val*sr)

        for node in self._feedbackSet:
            total = LogVar(0)
            for i in range(node.N):
                total += node._fullJoint[i]


    def sortaRootProb(self,conditioning):
        prod = LogVar(1)
        for sortaRoot in self._sortaRootDeps:
            j = conditioning[sortaRoot]
            prod *= self._L(sortaRoot,j)
            if(len(sortaRoot._parents) == 0):
                prod *= self._pi(sortaRoot,j)
            else:
                parents = sorted(sortaRoot._parents)
                X_ = [conditioning[p] for p in parents]
                prod *= self._trans(parents,sortaRoot,X_,j)
        return prod

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
        aLeaf = self._hyperGraph._leaves.__iter__().next()

        for X in itertools.product(*[range(n.N) for n in self._feedbackSet]):

            conditioning = {node:i for node,i in zip(self._feedbackSet,X)}

            for i in range(aLeaf.N):
                val = aLeaf.getUVal(i,conditioning)
                sr = self.sortaRootProb(conditioning)

                for node,x in zip(self._feedbackSet,X):
                    node.updateFullJoint(x,val*sr)

        for node in self._feedbackSet:
            total = LogVar(0)
            for i in range(node.N):
                total += node._fullJoint[i]

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

        # assert 0

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

        def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
            return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

        """ Make sure that P(Y) can be computed everywhere correctly """
        for node in self.nodes:
            total = LogVar(0)
            for i in range(node.N):

                uVal = node.getFullJoint(i)
                uCopy = LogVar(uVal)
                total += uCopy

            if(not isclose(total._logVal,correct._logVal)):
                if(printStuff):
                    print('\nFailed the P(Y) test!')
                    print('Node: '+str(node))
                    print('total: '+str(total.toFloat()))
                    print('correct: '+str(correct.toFloat()))
                    assert 0
            else:
                if(printStuff):
                    print('Passed P(Y) test for node '+str(node))
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
                    for X in itertools.product(*[range(p.N) for p in sorted(n._parents)]):

                        prob2 = self.probOfParentsProducingNode(n,X,i,correct)
                        total2 += prob2
                        shouldEqualProb += prob2

                    if(not isclose(shouldEqualProb.toFloat(),prob.toFloat())):
                        if(printStuff):
                            print('\n\n\n\nFailed this test for node: '+str(node))
                            print('This is wrong: '+str(shouldEqualProb.toFloat()))
                            print('This is right: '+str(prob.toFloat()))
                            assert 0
                    else:
                        if(printStuff):
                            print('Passed the marg out test P(x_c|y)=sumP(x_c,{x_p}|y)')
                            print('==========================================================')

            if(not isclose(total1.toFloat(),1.)):
                if(printStuff):
                    print('The sum of all P(x_i|Y) doesn\'t equal 1 for node: '+str(n)+', total1: '+str(total1._logVal))

            if(len(n._parents) > 0):
                if(not isclose(total2.toFloat(),1.)):
                    if(printStuff):
                        'The sum of all P({x_p}|Y) isn\'t 1 for node: '+str(n)+', total2: '+str(total2._logVal)

    def isolatedParentJoint(self,node,X,i,totalProb=None):

        # compute the probs for the fbs
        aLeaf = list(self._hyperGraph._leaves)[0]

        parents = sorted(node._parents)
        jointProb = LogVar(0)

        margOut = [n for n in self._feedbackSet if n!=node and n not in parents]

        for _X in itertools.product(*[range(n.N) for n in margOut]):

            conditioning = {_node:_i for _node,_i in zip(margOut,_X)}
            conditioning.update({p:j for p,j in zip(parents,X)})
            conditioning[node] = i

            for j in range(aLeaf.N):
                val = aLeaf.getUVal(j,conditioning)
                sr = self.sortaRootProb(conditioning)
                jointProb += val*sr

        if(totalProb):
            jointProb /= totalProb
        return jointProb


    def probOfParentsProducingNode(self,node,X,i,totalProb=None):

        if(node in self._sortaRootDeps):
            return self.isolatedParentJoint(node,X,i,totalProb)


        parents = sorted(node._parents)

        # print('\n\n-------------------\nCurrent node: '+str(node))
        # print('Parents: '+str(parents))
        # print('FBS: '+str(self._feedbackSet))
        # print('X: '+str(X))
        # print('i: '+str(i))

        # make sure we sum over all of the possible nodes in the fbs
        latentRanges = [range(n.N) for n in sorted(self._feedbackSet)]
        nodesInFBS = [n for n in sorted(self._feedbackSet)]
        for p,j in zip(parents,X):
            if(p.inFBS):
                latentRanges[nodesInFBS.index(p)] = [j]
        if(node.inFBS):
            latentRanges[nodesInFBS.index(node)] = [i]

        # print('nodes in FBS: '+str(nodesInFBS))
        prob = LogVar(1)

        """ Prob of this node """
        transProb = self._trans(parents,node,X,i)
        emissionProb = self._L(node,i)
        prob *= transProb
        prob *= emissionProb
        # print('Prob of this node: '+str(prob))

        total = LogVar(0)
        for S in itertools.product(*latentRanges):

            # print('\n\n\nS: '+str(S))

            familyInFBS = {n:s for n,s in zip(nodesInFBS,S)}

            _prob = LogVar(1)

            """ Down this node """
            for e in node._downEdges:
                if(node.inFBS): continue

                # print('\ncomputing v_ for node: '+str(node)+' e: '+str(e)+' i: '+str(i))
                v_ = node.getMarginalizedV(e,i,familyInFBS,self._feedbackSet)
                # print('vDeps: '+str(node._VDeps))
                # print('v_: '+str(v_))
                _prob *= v_
                # print('Prob is now: '+str(_prob))

            """ Out from each sibling """
            if(node._upEdge):

                for sibling in [x for x in node._upEdge._children if x!=node]:

                    # print('\ncomputing b_ for sibling: '+str(sibling)+' X: '+str(X))
                    b_ = sibling.getMarginalizedB(X,familyInFBS,self._feedbackSet)
                    # print('bDeps: '+str(sibling._bDeps))
                    # print('b_: '+str(b_))
                    _prob *= b_
                    # print('Prob is now: '+str(_prob))

                """ Out from each parent """
                for parent,j in zip(parents,X):
                    if(parent.inFBS): continue

                    # print('\ncomputing a_ for parent: '+str(parent)+' edge: '+str(node._upEdge)+' j: '+str(j))
                    a_ = parent.getMarginalizedA(node._upEdge,j,familyInFBS,self._feedbackSet)
                    # print('aDeps: '+str(parent._aDeps))
                    # print('a_: '+str(a_))
                    _prob *= a_
                    # print('Prob is now: '+str(_prob))


                spr = self.sortaRootProb(familyInFBS)
                # print('\nSorta root prob: '+str(spr))
                _prob *= spr

                total += _prob

        prob *= total



        # print('\n\nMultiplied together: '+str(prob))
        # print('Total prob: '+str(totalProb))

        """ Normalize """
        if(totalProb):
            prob /= totalProb

        return prob

    def getStats(self):

        for node in self.nodes:
            for i in range(node.N):
                uVal = node.getFullJoint(i)

        for n in self.nodes:
            for i in range(n.N):
                if(len(n._parents) > 0):
                    for X in itertools.product(*[range(p.N) for p in sorted(n._parents)]):
                        genProb = self.probOfParentsProducingNode(n,X,i)


