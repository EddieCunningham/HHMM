# distutils: language=c++

from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "HHMMMessagePassing.h":

    cdef cppclass MessagePassingHyperGraphWrapper:

        MessagePassingHyperGraphWrapper() except +
        MessagePassingHyperGraphWrapper( unsigned int latentStateSize ) except +
        void initialize()
        void addNodeId( unsigned int _id, unsigned int y )
        void addEdgeId( vector[ unsigned int ] parents, vector[ unsigned int ] children, unsigned int _id )
        void preprocessWithInt( vector[ unsigned int ] fbs )
        float logProbOfAllNodeObservations();
        void getStats();
        float logFullJoint( unsigned int _id, unsigned int i )

cdef class MessagePasser:

    cdef MessagePassingHyperGraphWrapper mp

    def __init__( self, latentStateSize ):
        self.mp = MessagePassingHyperGraphWrapper( latentStateSize )

    def initialize( self ):
        self.mp.initialize()

    def addNode( self, _id, y=0 ):
        self.mp.addNodeId( _id, y )

    def addEdge( self, parentIds=[], childIds=[], edgeId=-1 ):
        if( edgeId == -1 ):
            assert 0
        self.mp.addEdgeId( parentIds, childIds, edgeId )

    def preprocess( self, fbs ):
        self.mp.preprocessWithInt( fbs )

    def logObservationProb( self ):
        return self.mp.logProbOfAllNodeObservations()

    def getStats( self ):
        return self.mp.getStats()

    def logFullJoint( self, _id, i ):
        return self.mp.logFullJoint( _id, i )