#ifndef __HYPERGRAPH_H__
#define __HYPERGRAPH_H__

#define notInSet( s, x ) s.find( x ) == s.end()
#define inSet( s, x ) s.find( x ) != s.end()

#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>

class Node;
class Edge;

typedef std::unordered_map map;
typedef std::unordered_set set;

typedef Node* Node_ptr;
typedef Edge* Edge_ptr;

class Node {
public:

    set< Node_ptr >                  parents;
    map< Edge_ptr, set< Node_ptr > > childrenForEdge;
    Edge_ptr                         upEdge;
    set< Edge_ptr >                  downEdges;

    int  id     = -1;
    bool isRoot = false;
    bool isLeaf = false;

    Node():
          parents(),
          childrenForEdge(),
          upEdge(nullptr),
          downEdges() {}

    void addUpEdge  ( Edge_ptr edge );
    void addDownEdge( Edge_ptr edge );

    bool operator == (const Node& other) { return id == other.id; }
    bool operator != (const Node& other) { return id != other.id; }
};

class Edge {
public:

    set< Node_ptr > parents;
    set< Node_ptr > children;

    int id = -1;

    Edge():
            parents (),
            children() {}

    void addParent( Node_ptr node );
    void addChild ( Node_ptr node );

    bool operator == (const Node& other) { return id == other.id; }
    bool operator != (const Node& other) { return id != other.id; }
};

class HyperGraph {
public:

    map< int,Node > nodeIds;
    map< int,Edge > edgeIds;
    set< Node_ptr > leaves;
    set< Node_ptr > roots;

    bool initialized = false;

    Node_ptr addNode( int id );
    bool     hasNode( int id );
    Node_ptr getNode( int id );

    Edge_ptr addEdge( const std::vector< Node_ptr > & parents, int id );
    bool     hasEdge( int id );
    Edge_ptr getEdge( int id );

    void initialize();
};

std::ostream& operator<<( std::ostream &os, const Node& node ) {
    return os << node.id;
}

std::ostream& operator<<( std::ostream &os, const Edge& edge ) {
    return os << edge.id;
}



#endif