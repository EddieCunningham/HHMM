#include "HHMM.h"
#include <algorithm>
#include <cassert>

ConditioningProduct::ConditioningProduct( const std::vector< Node_ptr >& nodes ):
_nodes( nodes ) {
    _X = parentStates( 0, nodes.size() );
}

/* Increment the product by 1 */
bool ConditioningProduct::next() {

    for( uint i = 0; i < _nodes.size(); ++i ) {

        uint& val = _X.at( i );
        val += 1;
        if( val == _nodes.at( i )->N ) {
            val = 0;
        }
        else {
            return true;
        }
    }
    return false;
}

/* Return a map with the nodes and their current value */
conditioning ConditioningProduct::getConditioning() {

    conditioning cond = conditioning();
    for( uint i = 0; i < _nodes.size(); ++i ) {
        cond.emplace( _nodes.at( i ), _X.at( i ) );
    }
    return cond;
}

parentStates ConditioningProduct::getParentStates() {
    return _X;
}

/* Like Conditioning Product, but allows
   some values to be constant
*/
CartesianNodeProduct::CartesianNodeProduct( const std::vector< uint >& maxVals, const std::vector< std::pair< uint, uint > >& constValues ):
_maxVals( maxVals ) {

    _X = parentStates( 0, maxVals.size() + constValues.size() );
    for( const std::pair< uint, uint >& val : constValues ) {
        _X.at( val.first ) = val.second;
        _constIndices.push_back( val.first );
    }
}

bool CartesianNodeProduct::next() {

    for( uint i = 0; i < _X.size(); ++i ) {
        if( inVector( _constIndices, i ) ) {
            continue;
        }
        uint& val = _X.at( i );
        val += 1;
        if( val == _maxVals.at( i ) ) {
            val = 0;
        }
        else {
            return true;
        }
    }
    return false;
}

/* Get the states that we're incrementing over */
parentStates CartesianNodeProduct::getNonConstStates( int index ) {

    parentStates X = parentStates();
    for( int i = 0; i < (int)_X.size(); ++i ) {

        if( index == -1 ) {
            if( inVector( _constIndices, i ) ) {
                continue;
            }
        }
        else {
            if( i == index ) {
                continue;
            }
        }
        X.push_back( _X.at( i ) );
    }
    return X;
}

parentStates CartesianNodeProduct::getParentStates() {
    return _X;
}

/* A key for a node that will be conditioned on */
std::string keyString( const Node_ptr node, uint state ) {
    return std::to_string( node->id ) +":"+ std::to_string( state );
}

/* Get the key for condition nodes in cond.  We only need to condition
   on nodes that are in deps
*/
condKey getKey( const conditioning& cond, const set< Node_ptr >& deps ) {

    std::string tmpString;
    std::vector< std::string > keyBuffer = std::vector< std::string >();

    /* Conditioning only applies to the */
    /* intersection of cond and deps    */
    for( const std::pair< const Node_ptr, uint > &c : cond ) {

        const Node_ptr node = c.first;
        uint state          = c.second;

        if( inSet( deps, ( Node_ptr )node ) ) {

            /* Build the key string for this node and */
            /* add it to the key buffer               */
            tmpString = keyString( node, state );
            keyBuffer.push_back( tmpString );
        }
    }

    /* Sort the key buffer for consistency */
    std::sort( keyBuffer.begin(), keyBuffer.end() );

    condKey key = "";
    for( std::string &keySegment : keyBuffer ) {
        key += keySegment+" ";
    }

    return key;
}

/* Get the key without restricting to deps */
condKey getKey( const conditioning& cond ) {

    std::string tmpString;
    std::vector< std::string > keyBuffer = std::vector< std::string >();

    for( const std::pair< const Node_ptr, uint > &c : cond ) {

        const Node_ptr node = c.first;
        uint state    = c.second;

        /* Build the key string for this node and */
        /* add it to the key buffer               */
        tmpString = keyString( node, state );
        keyBuffer.push_back( tmpString );
    }

    /* Sort the key buffer */
    std::sort( keyBuffer.begin(), keyBuffer.end() );

    condKey key = "";
    for( std::string &keySegment : keyBuffer ) {
        key += keySegment+" ";
    }
    return key;
}


std::ostream& operator<<( std::ostream &os, const Node& node ) {
    return os << node.id;
}

std::ostream& operator<<( std::ostream &os, const Edge& edge ) {
    return os << edge.id;
}
