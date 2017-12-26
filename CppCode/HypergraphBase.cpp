#include "HypergraphBase.h"


void Node::addUpEdge( Edge_ptr edge ) {

    if( this->upEdge != nullptr ) {
        /* Model assumes that node can only come from */
        /* 1 nuclear family                           */
        exit(0);
    }

    this->upEdge = edge;
    edge->children.insert( this );
}

void Node::addDownEdge( Edge_ptr edge ) {

    if( notInSet( this->downEdges, edge ) ) {

        this->downEdges.insert( edge );
        edge->parents  .insert( this );

        this->childrenForEdge.insert(
                                    std::make_pair( edge, set< Node_ptr >() )
                                    );
    }
}

void Edge::addParent( Node_ptr parent ) {

    this->parents.insert( parent );
    parent->addDownEdge ( this   );

    for( Node_ptr child : children ) {
        child->parents.insert( parent );
        parent->childrenForEdge.at( this ).insert( child );
    }
}

void Edge::addChild( Node_ptr child ) {

    this->children.insert( child );
    child->addUpEdge     ( this  );

    for( Node_ptr parent : parents ) {
        child->parents.insert( parent );
    }
}

Node_ptr HyperGraph::addNode( uint id ) {

    if( initialized ) { exit( 0 ); }

    /* Don't re-add node */
    if( this->hasNode( id ) ) {
        return this->getNode( id );
    }

    Node node = Node();
    node.id   = id;

    this->nodeIds.emplace( id, node );
    Node_ptr node_ref = &( this->nodeIds.at( id ) );
    this->nodes.insert( node_ref );
    return node_ref;
}

bool HyperGraph::hasNode( uint id ) {
    return inSet( this->nodeIds, id );
}

Node_ptr HyperGraph::getNode( uint id ) {
    return &this->nodeIds.at( id );
}

Edge_ptr HyperGraph::addEdge( const std::vector< Node_ptr > & parents, uint id ) {

    if( initialized ) { exit( 0 ); }

    if( hasEdge( id ) ) {
        return getEdge( id );
    }

    Edge edge = Edge();
    edge.id = id;

    for( Node_ptr parent : parents ) {
        edge.addParent( parent );
    }
    this->edgeIds.emplace( id, edge );
    Edge_ptr edge_ref = &( this->edgeIds.at( id ) );
    return edge_ref;
}

bool HyperGraph::hasEdge( uint id ) {
    return inSet( this->edgeIds, id );
}

Edge_ptr HyperGraph::getEdge( uint id ) {
    return &this->edgeIds.at( id );
}

void HyperGraph::initialize() {

    for( std::pair< const uint, Node >& idNode : this->nodeIds ) {

        Node_ptr n = &idNode.second;

        /* Mark roots and leaves */
        if( n->parents.size() == 0 ) {
            roots.insert( n );
            n->isRoot = true;
        }
        if( n->downEdges.size() == 0 ) {
            leaves.insert( n );
            n->isLeaf = true;
        }
    }
    this->initialized = true;
}

std::ostream& operator<<( std::ostream &os, const Node& node ) {
    return os << node.id;
}

std::ostream& operator<<( std::ostream &os, const Edge& edge ) {
    return os << edge.id;
}
