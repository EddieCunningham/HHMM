# distutils: language=c++

from libcpp cimport bool
from libcpp.vector cimport vector
from cycleDetector import identifyCycles

cdef extern from "HypergraphBase.h":

    cdef cppclass MessagePassingHyperGraphWrapper[ EmissionType ]:

        MessagePassingHyperGraphWrapper() except +
        void addNodeId( unsigned int _id )
        void addEdgeId( vector[ unsigned int ]& parents, unsigned int _id )
        void preprocess()
        float logProbOfAllNodeObservations();
        bool test();
        void getStats();