from HypergraphBase import *
from pyLogVar import *
from cycleDetector import *
import numpy as np

class MarkovModelMessagePasser():

    def __init__(self,hg,hgParamFunction):

        self._hyperGraph = hg
        self._paramGenerator = hgParamFunction

        self.nodes = hg._nodes
        self.edges = hg._edges

        self._params = self._paramGenerator(self._hyperGraph)
        self._trans = self._params['transDist']
        self._L = self._params['emissionDist']

        self._pi = self._params['rootDist']

        self._conditionStack = [{}]
        self._WStack = [{}]

        feedbackSet,blockManager = identifyCycles(self._hyperGraph)
        blockManager.accumulateAllChains()

        self._blockManager = blockManager



    def getWBruteForce(self,nodes,X_):

        rootGroups = []
        transGroups = []
        allAncestors = list(nodes)
        current = list(nodes)
        while(len(current) > 0):
            newCurrent = []
            for _node in current:
                good = True
                # make sure all siblings are in current
                if(len(_node._parents) > 0 and _node not in nodes):
                    for sibling in _node._upEdge._children:
                        if(sibling not in current):
                            if(sibling not in newCurrent):
                                newCurrent.append(sibling)
                            good = False
                if(good):
                    if(len(_node._parents) == 0):
                        if(_node not in rootGroups):
                            rootGroups.append(_node)
                        if(_node not in allAncestors):
                            allAncestors.append(_node)
                    else:
                        if([sorted(list(_node._parents)),_node] not in transGroups):
                            transGroups.append([sorted(list(_node._parents)),_node])
                        for p in sorted(list(_node._parents)):
                            if(p not in allAncestors):
                                allAncestors.append(p)

                    for p in _node._parents:
                        if(p not in newCurrent):
                            newCurrent.append(p)
            current = newCurrent

        val = LogVar(0)
        for X in itertools.product(*[range(_node.N) if _node not in nodes else [X_[nodes.index(_node)]] for _node in allAncestors]):
            _val = LogVar(1)
            for root in rootGroups:
                _val *= self._pi(root,X[allAncestors.index(root)])
            for _parents,_node in transGroups:
                _i = X[allAncestors.index(_node)]
                _X = [X[allAncestors.index(p)] for p in _parents]
                _val *= self._trans(_parents,_node,_X,_i)
            val += _val
        return val


    def _pushW(self,node,i):

        currentConditioning = {}
        oldHeads = self._conditionStack[-1]
        currentConditioning.update(oldHeads)

        currentConditioning[node] = i

        self._conditionStack.append(currentConditioning)
        self._WStack.append({})

    def _popW(self):

        self._conditionStack.pop()
        self._WStack.pop()

    def _getN(self,node):
        # check if we've conditioned on this:
        if(node in self._conditionStack[-1]):
            return [self._conditionStack[-1][node]]
        return range(node.N)

    def _applyWConditioning(self,condX,conditionNodes):
        for node,i in zip(conditionNodes,condX):
            self._pushW(node,i)

    def _removeWConditioning(self,conditionNodes):
        for _ in conditionNodes:
            self._popW()


    def chainRule(self,nodes,X):

        latentNodes = self._blockManager.getLatentNodes(nodes)

        # make sure that all parents are included, even parents of
        # latent nodes that we will sum over
        nodesOfInterest = nodes+latentNodes

        latentParents = []

        for node in nodesOfInterest:

            missingParents = [p for p in node._parents if p not in nodesOfInterest and p not in self._conditionStack[-1]]

            if(len(missingParents) != len(node._parents) and\
               len(missingParents) != 0):

                for parent in missingParents:

                    if(parent not in self._conditionStack[-1] and\
                       parent not in latentParents and\
                       parent not in latentNodes):
                        latentParents.append(parent)

        if(len(latentNodes+latentParents) > 0):
            # then return so that we can sum over the latent
            # nodes and then next pass we can actually perform the
            # chain rule
            return [[nodes,X,[],[],latentNodes+latentParents]]

        heldOutEdges = {}
        chainRuleGroups = []

        # here, we know that there are no latent nodes to worry about
        # and we just need to worry about how to split up the nodes
        # such that the chain rule is satisfied
        for node in nodes:

            anyParentInNodes = [p for p in node._parents if p in nodes and p not in self._conditionStack[-1]]
            parents = []
            parentsX = []
            if(len(anyParentInNodes) > 0):

                for parent in node._parents:

                    if(parent not in parents and\
                       parent not in self._conditionStack[-1]):

                        parents.append(parent)
                        parentsX.append(X[nodes.index(parent)])

            if(node._upEdge is None):
                if(node not in self._conditionStack[-1]):
                    chainRuleGroups.append([[node],[X[nodes.index(node)]],[],[],[]])
            else:
                edgeId = node._upEdge._id
                if(edgeId not in heldOutEdges):
                    heldOutEdges[edgeId] = [[],[],parents,parentsX,[]]

                if(node not in self._conditionStack[-1]):
                    heldOutEdges[edgeId][0].append(node)
                    heldOutEdges[edgeId][1].append(X[nodes.index(node)])

        if(len(heldOutEdges) > 0):
            for v in heldOutEdges.values():
                if(len(v[0]) > 0):
                    chainRuleGroups.append(v)

        for c in chainRuleGroups: assert len(c[0]) > 0
        return chainRuleGroups


    def _jointProb(self,nodes,X):

        # using the chain rule, split the nodes up so that
        # we can condition on them
        chainRuleGroups = self.chainRule(nodes,X)

        """ DO THE COMPUTATION """
        jointProb = LogVar(1)
        for _nodes,_X,conditionNodes,condX,latentNodes in chainRuleGroups:

            _jointProb = LogVar(0)

            iterRange = [self._getN(_node) if _node not in conditionNodes else [_x] for _node,_x in zip(conditionNodes,condX)]
            for _condX in itertools.product(*iterRange):

                if(len(latentNodes) > 0):

                    _iterRange = [self._getN(l)for l in latentNodes]
                    for _latentX in itertools.product(*_iterRange):

                        self._applyWConditioning(_condX,conditionNodes)

                        __jointProb = self._jointProb(_nodes+latentNodes,_X+list(_latentX))

                        self._removeWConditioning(conditionNodes)

                        _jointProb += __jointProb

                else:

                    self._applyWConditioning(_condX,conditionNodes)

                    w = self._computeAssumeIndep(_nodes,_X)

                    self._removeWConditioning(conditionNodes)
                    conditionalProb = w
                    _jointProb += conditionalProb

            jointProb *= _jointProb

        return jointProb

    def _computeAssumeIndep(self,nodes,X):

        upEdge = nodes[0]._upEdge
        assert len([_ for n in nodes if n._upEdge != upEdge]) == 0

        w = LogVar(1)
        parents = sorted(list(nodes[0]._parents))
        inCondStack = [node for node in nodes if node in self._conditionStack[-1]]

        if(len(parents) == 0):
            for node,i in zip(nodes,X):

                assert len(node._parents) == 0

                if(node in inCondStack):
                    w *= int(i==self._conditionStack[-1][node])
                else:
                    w *= self._pi(node,i)

        else:

            _w = LogVar(0)
            for _X in itertools.product(*[self._getN(p) for p in parents]):

                jointParentProb = self.getJoint(parents,_X)
                transProb = LogVar(1)

                for node,i in zip(nodes,X):
                    if(node in inCondStack):
                        continue
                    transProb *= self._trans(parents,node,_X,i)

                _w += transProb * jointParentProb

            w *= _w

        return w

    def getJoint(self,nodes,X):

        if(len(nodes) == 0):
            return LogVar(1)

        nodes = list(nodes)
        X = list(X)

        W = self._WStack[-1]
        key1 = str(nodes)
        key2 = str(X)


        """ W_N(X) = P(N_x=X) """
        if(key1 not in W or key2 not in W[key1]):
            if(key1 not in W):
                W[key1] = {}

            """ Compute W """
            W[key1][key2] = self._jointProb(nodes,X)
        return W[key1][key2]

    def getW(self,node,i):
        return self.getJoint([node],[i])

