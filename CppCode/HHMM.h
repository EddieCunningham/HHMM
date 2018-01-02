#ifndef __HYPERGRAPH_H__
#define __HYPERGRAPH_H__

#define notInSet( s, x ) s.find( x ) == s.end()
#define inSet( s, x )    s.find( x ) != s.end()
#define notInVector( v, x ) std::find( v.begin(), v.end(), x ) == v.end()
#define inVector( v, x )    std::find( v.begin(), v.end(), x ) != v.end()

#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
#include "LogVar.h"

#define map std::unordered_map
#define set std::set

class Node;
class Edge;
class HyperGraph;

typedef Node*       Node_ptr;
typedef Edge*       Edge_ptr;
typedef HyperGraph* HyperGraph_ptr;

typedef std::string                          condKey;
typedef std::unordered_map< Node_ptr, uint > conditioning;
typedef std::vector< uint >                  parentStates;

typedef struct {
    /* Hash function for vectors of uints taken from                                  */
    /* https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector */
    std::size_t operator()( std::vector< uint > const& vec ) const {
        std::size_t seed = vec.size();
        for( auto& i : vec ) {
            seed ^= i + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
        }
        return seed;
    }
} StateHash;

/* ================================================================================================================== */

class ConditioningProduct {
private:
    parentStates _X;
    std::vector< Node_ptr > _nodes;
public:
    ConditioningProduct( const std::vector< Node_ptr >& nodes );

    bool next();
    conditioning getConditioning();
    parentStates getParentStates();
};

/* Like Conditioning Product, but allows
   some values to be constant
*/
class CartesianNodeProduct {
private:
    parentStates        _X;
    std::vector< uint > _maxVals;
    std::vector< uint > _constIndices;
public:
    CartesianNodeProduct( const std::vector< uint >& maxVals, const std::vector< std::pair< uint, uint > >& constValues );
    bool next();

    /* Get the states that we're incrementing over */
    parentStates getNonConstStates( int index = -1 );
    parentStates getParentStates();
};

std::string keyString( const Node_ptr node, uint state );

/* Get the key for condition nodes in cond.  We only need to condition
   on nodes that are in deps
*/
condKey getKey( const conditioning& cond, const set< Node_ptr >& deps );

/* Get the key without restricting to deps */
condKey getKey( const conditioning& cond );

/* ================================================================================================================== */


class Node {
private:

    HyperGraph_ptr _msg;

    map< Edge_ptr    , map< uint   , map< condKey, LogVar > > > _a;
    map< parentStates, map< condKey, LogVar >, StateHash >      _b;
    map< uint        , map< condKey, LogVar > >                 _U;
    map< Edge_ptr    , map< uint   , map< condKey, LogVar > > > _V;

    set< Node_ptr > _aDependencies;
    set< Node_ptr > _bDependencies;
    set< Node_ptr > _UDependencies;
    set< Node_ptr > _VDependencies;

    std::vector< LogVar > _fullJoint;

    LogVar transFunc   ( const parentStates& X, uint k );
    LogVar emissionFunc( uint i                        );
    LogVar rootFunc    ( uint i                        );

    condKey _aKey            ( const conditioning& cond                                                                 );
    bool    _needToComputeA  ( Edge_ptr edge, uint i, condKey key                                                       );
    void    _setAVal         ( Edge_ptr edge, uint i, condKey key, LogVar aVal                                          );
    LogVar  _getAVal         ( Edge_ptr edge, uint i, condKey key                                                       );
    LogVar  _getA            ( Edge_ptr edge, uint i, const conditioning& cond                                          );
    LogVar  _computeA        ( Edge_ptr edge, uint i, const conditioning& cond                                          );
    LogVar  _getMarginalizedA( Edge_ptr edge, uint i, const conditioning& nodesToKeep, const set<Node_ptr>& feedbackSet );
    bool    _aReady          ( Edge_ptr edge, uint i, const conditioning& cond                                          );

    condKey _bKey            ( const conditioning& cond                                                          );
    bool    _needToComputeB  ( parentStates X, condKey key                                                       );
    void    _setBVal         ( parentStates X, condKey key, LogVar bVal                                          );
    LogVar  _getBVal         ( parentStates X, condKey key                                                       );
    LogVar  _getB            ( parentStates X, const conditioning& cond                                          );
    LogVar  _computeB        ( parentStates X, const conditioning& cond                                          );
    LogVar  _getMarginalizedB( parentStates X, const conditioning& nodesToKeep, const set<Node_ptr>& feedbackSet );
    bool    _bReady          ( parentStates X, const conditioning& cond                                          );

    condKey _UKey            ( const conditioning& cond                                                  );
    bool    _needToComputeU  ( uint i, condKey key                                                       );
    void    _setUVal         ( uint i, condKey key, LogVar uVal                                          );
    LogVar  _getUVal         ( uint i, condKey key                                                       );
    LogVar  _getU            ( uint i, const conditioning& cond                                          );
    LogVar  _computeU        ( uint i, const conditioning& cond                                          );
    LogVar  _getMarginalizedU( uint i, const conditioning& nodesToKeep, const set<Node_ptr>& feedbackSet );
    bool    _UReady          ( uint i, const conditioning& cond                                          );

    condKey _VKey            ( const conditioning& cond                                                                 );
    bool    _needToComputeV  ( Edge_ptr edge, uint i, condKey key                                                       );
    void    _setVVal         ( Edge_ptr edge, uint i, condKey key, LogVar vVal                                          );
    LogVar  _getVVal         ( Edge_ptr edge, uint i, condKey key                                                       );
    LogVar  _getV            ( Edge_ptr edge, uint i, const conditioning& cond                                          );
    LogVar  _computeV        ( Edge_ptr edge, uint i, const conditioning& cond                                          );
    LogVar  _getMarginalizedV( Edge_ptr edge, uint i, const conditioning& nodesToKeep, const set<Node_ptr>& feedbackSet );
    bool    _VReady          ( Edge_ptr edge, uint i, const conditioning& cond                                          );

    LogVar _sortaRootProb      ( const conditioning& cond         );
    void   _accumulateFullJoint( const set<Node_ptr>& feedbackSet );
    void   _updateFullJoint    ( uint i, LogVar val               );

    std::pair< uint, uint > _getN( Node_ptr node , const conditioning& cond );

public:

    set< Node_ptr >                  parents;
    map< Edge_ptr, set< Node_ptr > > childrenForEdge;
    Edge_ptr                         upEdge;
    set< Edge_ptr >                  downEdges;

    uint  id    = -1;
    bool isRoot = false;
    bool isLeaf = false;

    Node():
        _aDependencies(),
        _bDependencies(),
        _UDependencies(),
        _VDependencies(),
        parents( set< Node_ptr >() ),
        childrenForEdge( map< Edge_ptr, set< Node_ptr > >() ),
        upEdge( nullptr ),
        downEdges( set< Edge_ptr >() ) {
            reset();
        }

    uint N;
    uint y;
    bool inFeedbackSet;

    void    reset            (        );
    LogVar  getFullJoint     ( uint i );

    void addUpEdge  ( Edge_ptr edge );
    void addDownEdge( Edge_ptr edge );

    bool operator == (const Node& other) { return id == other.id; }
    bool operator != (const Node& other) { return id != other.id; }

    friend class Edge;
    friend class HyperGraph;
};

class Edge {
public:

    set< Node_ptr > parents;
    set< Node_ptr > children;

    uint id = -1;

    Edge():
        parents ( set< Node_ptr >() ),
        children( set< Node_ptr >() ) {}

    void addParent( Node_ptr node );
    void addChild ( Node_ptr node );

    bool operator == (const Node& other) { return id == other.id; }
    bool operator != (const Node& other) { return id != other.id; }

    friend class Node;
    friend class HyperGraph;
};

class HyperGraph {
private:

    conditioning           _conditioning;
    map< condKey, LogVar > _sortaRootProbs;
    set< Node_ptr >        _sortaRootDeps;
    bool                   _preprocessing;

    std::vector< std::vector< float > > pi;
    std::vector< std::vector< float > > L;
    std::vector< std::vector< std::vector < float > > > trans;

    void   _computeForPreprocessing(                          );
    LogVar _sortaRootProb          ( const conditioning& cond );

public:

    map< uint, Node > nodeIds;
    map< uint, Edge > edgeIds;
    set< Node_ptr >  leaves;
    set< Node_ptr >  roots;
    set< Node_ptr >  nodes;
    set< Edge_ptr >  edges;

    set< Node_ptr > feedbackSet;
    uint N;

    bool initialized = false;

    HyperGraph():
        _conditioning(),
        _sortaRootProbs(),
        _sortaRootDeps(),
        nodeIds( map< uint, Node >() ),
        edgeIds( map< uint, Edge >() ),
        leaves( set< Node_ptr >() ),
        roots( set< Node_ptr >() ),
        nodes( set< Node_ptr >() ),
        edges( set< Edge_ptr >() ),
        feedbackSet(),
        N( -1 ) {}

    HyperGraph( uint latentStateSize ):
        _conditioning(),
        _sortaRootProbs(),
        _sortaRootDeps(),
        nodeIds( map< uint, Node >() ),
        edgeIds( map< uint, Edge >() ),
        leaves( set< Node_ptr >() ),
        roots( set< Node_ptr >() ),
        nodes( set< Node_ptr >() ),
        edges( set< Edge_ptr >() ),
        feedbackSet(),
        N( latentStateSize ) {}

    uint indexOfNode( uint id );
    uint indexOfRoot( uint id );

    Node_ptr addNode( uint id );
    bool     hasNode( uint id );
    Node_ptr getNode( uint id );

    Edge_ptr addEdge( const std::vector< Node_ptr > & parents, uint id );
    bool     hasEdge( uint id );
    Edge_ptr getEdge( uint id );

    void initialize();


    void addNodeId( uint id, uint y );
    void addEdgeId( const std::vector< uint > & parents, const std::vector< uint > & children, uint id );

    void     preprocess                  ( const set< Node_ptr >& feedbackSet     );
    LogVar   isolatedParentJoint         ( Node_ptr node, parentStates X, uint i  );
    LogVar   jointParentChild  ( Node_ptr node, parentStates X, uint i  );
    void     getStats                    (                                        );
    LogVar   probOfAllNodeObservations   (                                        );
    float    logProbOfAllNodeObservations(                                        );
    void     getCounts                   (                                        );
    void     preprocessWithInt           ( const std::vector< uint >& feedbackSet );
    float    logFullJoint                ( uint nodeId, uint i                    );


    void setPi   ( std::vector< std::vector< float > > _pi );
    void setL    ( std::vector< float > _L  );
    void setTrans( std::vector< std::vector< std::vector < float > > > _trans );

    friend class Node;
    friend class Edge;
};

std::ostream& operator<<( std::ostream &os, const Node& node );
std::ostream& operator<<( std::ostream &os, const Edge& edge );

#endif
