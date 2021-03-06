# distutils: language=c++

from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "HHMM.h":

    ctypedef vector[ unsigned int ] parentStates

    cdef cppclass HyperGraph:

        HyperGraph() except +
        HyperGraph( unsigned int latentStateSize ) except +
        void initialize()
        void addNodeId( unsigned int _id, unsigned int y )
        void addEdgeId( vector[ unsigned int ] parents, vector[ unsigned int ] children, unsigned int _id )
        void preprocessWithInt( vector[ unsigned int ] fbs )
        float logProbOfAllNodeObservations()
        void getStats()
        float logFullJoint( unsigned int _id, unsigned int i )

        unsigned int indexOfRoot( unsigned int id )

        void setPi( np.ndarray[ unsigned int, ndim=2 ] pi )
        void setL( np.ndarray[ unsigned int, ndim=1 ] L )
        void setTrans( np.ndarray[ unsigned int, ndim=3 ] trans )

cdef class MessagePasser:

    cdef HyperGraph mp

    def __init__( self, latentStateSize ):
        self.mp = HyperGraph( latentStateSize )

    def initialize( self ):
        self.mp.initialize()

    def addNode( self, _id, y=0 ):
        self.mp.addNodeId( _id, y )

    def addEdge( self, parentIds=[], childIds=[], edgeId=-1 ):
        if( edgeId == -1 ):
            assert 0
        self.mp.addEdgeId( parentIds, childIds, edgeId )

    cpdef setRootDistribution( self, vector[ vector[ unsigned int ] ] rootDist ):
        self.mp.setPi( rootDist )

    cpdef setEmissionDistribution( self, vector[ unsigned int ] emissionDist ):
        self.mp.setL( emissionDist )

    cpdef setTransitionDistribution( self, vector[ vector[ vector[ unsigned int ] ] ] transDist ):
        self.mp.setTrans( transDist )

    def preprocess( self, fbs ):
        self.mp.preprocessWithInt( fbs )

    def logObservationProb( self ):
        return self.mp.logProbOfAllNodeObservations()

    def getStats( self ):
        return self.mp.getStats()

    def logFullJoint( self, _id, i ):
        return self.mp.logFullJoint( _id, i )