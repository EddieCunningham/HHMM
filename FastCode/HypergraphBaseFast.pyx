cdef class Node:

    cdef public set _parents
    cdef public dict _childrenForEdge
    cdef public Edge _upEdge
    cdef public set _downEdges
    cdef public int _id
    cdef public bint isRoot
    cdef public bint isLeaf

    def __init__(self):
        self._parents = set()
        self._childrenForEdge = {}
        self._upEdge = None
        self._downEdges = set()

        self.isRoot = False
        self.isLeaf = False

    def __hash__(self):
        return hash('n'+str(self._id))

    cpdef addUpEdge(self,Edge edge):
        """ Add edge to upEdges and add self to edge's children """
        if(self._upEdge != None): return
        self._upEdge = edge
        edge._children.add(self)

    cpdef addDownEdge(self,Edge edge):
        """ Add edge to downEdges and add self to edge's parents """
        if(edge not in self._downEdges):
            self._downEdges.add(edge)
            edge._parents.add(self)
            self._childrenForEdge[edge] = set()

    def __repr__(self):
        return str(self._id)

    def __richcmp__(self,other,op):
        if(op == 0):
            return hash(self._id) < hash(other._id)
        elif(op == 1):
            return hash(self._id) <= hash(other._id)
        elif(op == 2):
            return hash(self._id) == hash(other._id)
        elif(op == 3):
            return hash(self._id) != hash(other._id)
        elif(op == 4):
            return hash(self._id) > hash(other._id)
        elif(op == 5):
            return hash(self._id) >= hash(other._id)


cdef class Edge:
    cdef public set _parents
    cdef public set _children
    cdef public int _id

    def __init__(self):
        self._parents = set()
        self._children = set()

    def __hash__(self):
        return hash('e'+str(self._id))

    cpdef addParent(self,Node node):
        for child in self._children:
            child._parents.add(node)
        node.addDownEdge(self)
        node._childrenForEdge[self] |= self._children

    cpdef addChild(self,Node node):
        node._parents |= self._parents
        node.addUpEdge(self)

    def __repr__(self):
        return str(self._id)

cdef class BaseHyperGraph:

    cdef public set _nodes
    cdef public set _leaves
    cdef public set _roots
    cdef public set _edges
    cdef public bint _initialized
    cdef public object _NodeType
    cdef public object _EdgeType
    cdef public dict _nodeIDs
    cdef public dict _edgeIDs

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
        return
