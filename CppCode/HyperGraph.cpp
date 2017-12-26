#include "HHMM.h"
#include <algorithm>
#include <cassert>
#include <iostream>


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


/* ------------------------------------------------------------------------------------ */

LogVar HyperGraph::_sortaRootProb( const conditioning& cond ) {

    condKey key = getKey( cond );

    if( inSet( _sortaRootProbs, key ) ) {
        return _sortaRootProbs.at( key );
    }

    LogVar prod( 1 );

    if( _sortaRootDeps.size() > 0 ) {

        for( Node_ptr const &sortaRoot : _sortaRootDeps ) {

            uint j = cond.at( sortaRoot );
            LogVar innerProd = LogVar( ( *_L )( j ) );

            if( sortaRoot->parents.size() == 0 ) {
                innerProd *= LogVar( ( *_pi )( j ) );
            }
            else {
                parentStates X = parentStates();
                for( const std::pair< const Node_ptr, uint > &c : cond ) {
                    X.push_back( c.second );
                }
                innerProd *= LogVar( ( *_trans )( X, j ) );
            }
            prod *= innerProd;
        }
    }
    _sortaRootProbs.at( key ) = prod;
    return prod;
}

void HyperGraph::_computeForPreprocessing() {

    for( const Node_ptr& _node : nodes ) {
        Node_ptr node = ( Node_ptr )_node;

        if( node->inFeedbackSet ) {
            continue;
        }

        ConditioningProduct xIter = ConditioningProduct( std::vector< Node_ptr >( feedbackSet.begin(), feedbackSet.end() ) );
        do {
            conditioning cond = xIter.getConditioning();

            for( uint i = 0; i < node->N; ++ i ) {

                node->_getU( i, cond );

                for( const Edge_ptr& edge : node->downEdges ) {

                    node->_getV( edge, i, cond );
                }
            }
        } while( xIter.next() );
    }
    for( const Node_ptr& _node : nodes ) {
        Node_ptr node = ( Node_ptr )_node;

        if( node->inFeedbackSet ) { continue; }

        node->_accumulateFullJoint( feedbackSet );
    }

    Node_ptr aLeaf = ( Node_ptr )( *( leaves.begin() ) );

    _sortaRootProbs = map< condKey, LogVar >();

    ConditioningProduct xIter = ConditioningProduct( std::vector< Node_ptr >( feedbackSet.begin(), feedbackSet.end() ) );

    do {
        conditioning cond = xIter.getConditioning();
        parentStates X = xIter.getParentStates();

        for( uint i = 0; i < aLeaf->N; ++i ) {
            LogVar val = aLeaf->_getU( i, cond );
            LogVar srProb = _sortaRootProb( cond );

            uint j = 0;
            for( Node_ptr const &node : feedbackSet ) {

                uint x = X.at( j );
                node->_updateFullJoint( x, val * srProb );
                ++j;
            }
        }
    } while( xIter.next() );

    for( Node_ptr const &node : feedbackSet ) {
        LogVar total( 0 );
        for( uint i = 0; i < node->N; ++i ) {
            total += node->_fullJoint.at( i );
        }
    }
}

/* ------------------------------------------------------------------------------------ */

void HyperGraph::preprocess( const set< Node_ptr >& feedbackSet ) {
    this->feedbackSet = feedbackSet;

    for( Node_ptr const &_node : nodes ) {
        Node_ptr node = ( Node_ptr )_node;
        node->_msg = this;
    }

    _sortaRootDeps = set< Node_ptr >();
    for( Node_ptr const &node : feedbackSet ) {
        if( node->parents.size() == 0 ) {
            _sortaRootDeps.insert( node );
            continue;
        }

        bool notSRP = false;
        for( Node_ptr const &_parent : node->parents ) {
            Node_ptr parent = ( Node_ptr )_parent;

            if( notInSet( feedbackSet, parent ) ) {
                notSRP = true;
                break;
            }
        }
        if( notSRP ) {
            continue;
        }
        for( Node_ptr const &_sibling : node->upEdge->children ) {
            Node_ptr sibling = ( Node_ptr )_sibling;

            if( notInSet( feedbackSet, sibling ) ) {
                notSRP = true;
                break;
            }
        }
        if( notSRP ) {
            continue;
        }
        _sortaRootDeps.insert( node );
    }

    for( Node_ptr const &_node : nodes ) {
        Node_ptr node = ( Node_ptr )_node;

        if( inSet( feedbackSet, node ) ) {
            node->inFeedbackSet = true;
        }
    }

    _preprocessing = true;
    _computeForPreprocessing();
    _preprocessing = false;
}

void HyperGraph::getCounts() {

    for( Node_ptr const &_node : nodes ) {
        Node_ptr node = ( Node_ptr )_node;
        node->reset();
    }

    // make this a true message passing
    // algorithm later!
    preprocess( feedbackSet );
}

LogVar HyperGraph::isolatedParentJoint( Node_ptr node, parentStates X, uint i ) {

    Node_ptr aLeaf = ( Node_ptr )( *( leaves.begin() ) );

    LogVar jointProb( 0 );

    std::vector< Node_ptr > toMarginalizeOut = std::vector< Node_ptr >();
    for( Node_ptr const &_node : feedbackSet ) {
        if( node != _node && notInSet( node->parents, _node ) ) {
            toMarginalizeOut.push_back( node );
        }
    }

    conditioning parentCond = conditioning();
    uint j = 0;
    for( Node_ptr const &_parent : node->parents ) {
        Node_ptr parent = ( Node_ptr )_parent;

        parentCond.emplace( parent, X.at( j ) );
        ++j;
    }

    ConditioningProduct xIter = ConditioningProduct( toMarginalizeOut );
    do {
        conditioning cond = xIter.getConditioning();
        cond.insert( parentCond.begin(), parentCond.end() );
        cond.emplace( node, i );

        for( uint j = 0; j < aLeaf->N; ++j ) {
            LogVar val = aLeaf->_getU( j, cond );
            LogVar srProb = _sortaRootProb( cond );
            jointProb += val * srProb;
        }
    } while( xIter.next() );

    return jointProb;
}

LogVar HyperGraph::probOfParentsProducingNode( Node_ptr node, parentStates X, uint i ) {

    if( inSet( _sortaRootDeps, node ) ) {
        return isolatedParentJoint( node, X, i );
    }

    std::vector< uint > fbsNotInParents = std::vector< uint >();
    std::vector< std::pair< uint, uint > > constValues = std::vector< std::pair< uint, uint > >();

    uint k = 0;
    for( Node_ptr const&fbsNode : feedbackSet ) {
        if( inSet( node->parents, fbsNode ) ) {

            uint index = std::distance( node->parents.begin(),
                                        std::find( node->parents.begin(), node->parents.end(), ( Node_ptr )fbsNode )
                                       );
            if( index == node->parents.size() ) {
                assert( 0 );
            }

            constValues.push_back( std::make_pair( k, X.at( index ) ) );
        }
        else {
            fbsNotInParents.push_back( fbsNode->N );
        }
        if( fbsNode == node ) {
            constValues.push_back( std::make_pair( k, i ) );
        }
        ++k;
    }

    LogVar prob( 0 );

    CartesianNodeProduct latentRanges = CartesianNodeProduct( fbsNotInParents, constValues );
    do {
        parentStates _X = latentRanges.getParentStates();
        conditioning cond = conditioning();
        uint k = 0;
        for( Node_ptr const& fbsNode : feedbackSet ) {
            cond.emplace( fbsNode, _X.at( k ) );
            ++k;
        }

        LogVar prod( 1 );
        for( Edge_ptr const &edge : node->downEdges ) {
            prod *= node->_getMarginalizedV( edge, i, cond, feedbackSet );
        }

        for( Node_ptr const &_sibling : node->upEdge->children ) {
            Node_ptr sibling = ( Node_ptr )_sibling;
            if( sibling == node ) {
                continue;
            }
            prod *= sibling->_getMarginalizedB( X, cond, feedbackSet );
        }

        uint j = 0;
        for( Node_ptr const &_parent : node->parents ) {
            Node_ptr parent = ( Node_ptr )_parent;

            prod *= parent->_getMarginalizedA( node->upEdge, X.at( j ), cond, feedbackSet );
            ++j;
        }

        prod *= _sortaRootProb( cond );
        prob += prod;

    } while( latentRanges.next() );

    LogVar transProb    = LogVar( ( *_trans )( X, k ) );
    LogVar emissionProb = LogVar( ( *_L )( k ) );
    prob *= transProb * emissionProb;

    return prob;
}

void HyperGraph::getStats() {

    getCounts();

    for( Node_ptr const &_node : nodes ) {
        Node_ptr node = ( Node_ptr )_node;

        for( uint i = 0; i < node->N; ++i ) {
        }
    }
}

LogVar HyperGraph::probOfAllNodeObservations() {

    Node_ptr aLeaf = ( Node_ptr )( *( leaves.begin() ) );

    LogVar total( 0 );
    for( uint i = 0; i < aLeaf->N; ++i ) {
        total += aLeaf->getFullJoint( i );
    }
    return total;
}

float HyperGraph::logProbOfAllNodeObservations() {

    return probOfAllNodeObservations().logValue();
}

void HyperGraph::preprocessWithInt( const std::vector< uint >& feedbackSet ) {

    set< Node_ptr > fbs = set< Node_ptr >();
    for( uint _id : feedbackSet ) {
        Node_ptr node = ( Node_ptr )getNode( _id );
        fbs.insert( node );
    }
    preprocess( fbs );
}

float HyperGraph::logFullJoint( uint nodeId, uint i) {
    Node_ptr node = ( Node_ptr )getNode( nodeId );
    return node->getFullJoint( i ).logValue();
}

void HyperGraph::addNodeId( uint id, uint y ) {

    Node_ptr node = ( Node_ptr )HyperGraph::addNode( id );
    node->N = N;
    node->y = y;
    node->_msg = this;
}

void HyperGraph::addEdgeId( const std::vector< uint > & parents, const std::vector< uint > & children, uint id ) {

    std::vector< Node_ptr > parentVect = std::vector< Node_ptr >();
    for( uint _id : parents ) {
        Node_ptr node = getNode( _id );
        parentVect.push_back( node );
    }
    Edge_ptr edge = HyperGraph::addEdge( parentVect, id );

    for( uint _id : children ) {
        Node_ptr node = getNode( _id );
        edge->addChild( node );
    }
}
