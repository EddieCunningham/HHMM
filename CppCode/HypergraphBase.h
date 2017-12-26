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

#define map std::unordered_map
#define set std::set

class Node;
class Edge;

typedef Node* Node_ptr;
typedef Edge* Edge_ptr;

class Node {
public:

    set< Node_ptr >                  parents;
    map< Edge_ptr, set< Node_ptr > > childrenForEdge;
    Edge_ptr                         upEdge;
    set< Edge_ptr >                  downEdges;

    uint  id    = -1;
    bool isRoot = false;
    bool isLeaf = false;

    Node():
          parents( set< Node_ptr >() ),
          childrenForEdge( map< Edge_ptr, set< Node_ptr > >() ),
          upEdge( nullptr ),
          downEdges( set< Edge_ptr >() ) {}

    void addUpEdge  ( Edge_ptr edge );
    void addDownEdge( Edge_ptr edge );

    bool operator == (const Node& other) { return id == other.id; }
    bool operator != (const Node& other) { return id != other.id; }
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
};

class HyperGraph {
public:

    map< uint, Node > nodeIds;
    map< uint, Edge > edgeIds;
    set< Node_ptr >  leaves;
    set< Node_ptr >  roots;
    set< Node_ptr >  nodes;
    set< Edge_ptr >  edges;

    bool initialized = false;

    HyperGraph():
        nodeIds( map< uint, Node >() ),
        edgeIds( map< uint, Edge >() ),
        leaves( set< Node_ptr >() ),
        roots( set< Node_ptr >() ),
        nodes( set< Node_ptr >() ),
        edges( set< Edge_ptr >() ) {}

    Node_ptr addNode( uint id );
    bool     hasNode( uint id );
    Node_ptr getNode( uint id );

    Edge_ptr addEdge( const std::vector< Node_ptr > & parents, uint id );
    bool     hasEdge( uint id );
    Edge_ptr getEdge( uint id );

    void initialize();
};

std::ostream& operator<<( std::ostream &os, const Node& node );
std::ostream& operator<<( std::ostream &os, const Edge& edge );

#endif