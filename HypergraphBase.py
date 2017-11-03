import graphviz
import itertools
import random
import operator

class Node(object):
    def __init__(self):
        self._parents = set()
        self._childrenForEdge = {}
        self._upEdge = None
        self._downEdges = set()
        self._id = -1

        self._cycleHeads = set()
        self._cycleBases = set()

        self.isRoot = False
        self.isLeaf = False

    def addUpEdge(self,edge):
        """ Add edge to upEdges and add self to edge's children """
        if(self._upEdge != None):
            assert self._upEdge == edge, 'Edge before: '+str(self._upEdge)+' but tried setting: '+str(edge)
            return
        self._upEdge = edge
        edge._children.add(self)

    def addDownEdge(self,edge):
        """ Add edge to downEdges and add self to edge's parents """
        if(edge not in self._downEdges):
            self._downEdges.add(edge)
            edge._parents.add(self)
            self._childrenForEdge[edge] = set()

    def __repr__(self):
        return str(self._id)

class Edge(object):
    def __init__(self):
        self._parents = set()
        self._children = set()
        self._id = -1

    def addParent(self,node):
        for child in self._children:
            child._parents.add(node)
        node.addDownEdge(self)
        node._childrenForEdge[self] |= self._children

    def addChild(self,node):
        node._parents |= self._parents
        node.addUpEdge(self)

    def __repr__(self):
        return str(self._id)

class BaseHyperGraph(object):
    def __init__(self):
        self._nodes = set()
        self._leaves = set()
        self._roots = set()
        self._edges = set()
        self._initialized = False
        self._NodeType = Node
        self._EdgeType = Edge
        self._nodeIDs = {}
        self._edgeIDs = {}

    def setNodeType(self,NodeType):
        self._NodeType = NodeType

    def setEdgeType(self,EdgeType):
        self._EdgeType = EdgeType

    def addNode(self,ID,*args):
        if(self._initialized):
            assert 0, 'Graph already initialized'
        node = self._NodeType(*args)
        node._id = ID
        self._nodes.add(node)
        assert ID not in self._nodeIDs
        self._nodeIDs[ID] = node
        return node

    def hasNode(self,ID):
        return ID in self._nodeIDs

    def getNode(self,ID):
        return self._nodeIDs[ID]

    def addEdge(self,parents,ID):
        if(self._initialized):
            assert 0, 'Graph already initialized'
        newEdge = self._EdgeType()
        newEdge._id = ID
        self._edges.add(newEdge)
        [newEdge.addParent(p) for p in parents]
        assert newEdge not in self._edgeIDs
        self._edgeIDs[ID] = newEdge
        return newEdge

    def hasEdge(self,ID):
        return ID in self._edgeIDs

    def getEdge(self,ID):
        return self._edgeIDs[ID]

    def initialize(self):
        self._initialized = True

        for n in self._nodes:
            if(len(n._parents) == 0):
                self._roots.add(n)
                n.isRoot = True
            if(len(n._downEdges) == 0):
                self._leaves.add(n)
                n.isLeaf = True

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

