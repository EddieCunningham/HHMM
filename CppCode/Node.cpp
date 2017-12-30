#include "HHMM.h"
#include <algorithm>
#include <cassert>
#include <iostream>

/* Simulates
    a = set(...)
    b = set(...)
    a |= b
*/
template < typename T >
void inplace_set_union( set<T>& a, const set<T>& b ) {
    for( T const &elt : b ) {
        a.insert( elt );
    }
}


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


/* Get the range of latent states node can have when conditioned
   on cond
*/
std::pair< uint, uint > Node::_getN( Node_ptr node , const conditioning& cond ) {
    if( inSet( cond, node ) ) {
        return std::make_pair( cond.at(node), cond.at(node) );
    }
    return std::make_pair( 0, node->N) ;
}

/* ------------------------------------------------------------------------- */

/* Get conditioning key for the a value */
condKey Node::_aKey( const conditioning& cond ) {
    return getKey( cond, _aDependencies );
}


/* Check if we need to compute the a value */
bool Node::_needToComputeA( Edge_ptr edge, uint i, condKey key ) {

    if( inSet( _a                   , edge ) &&
        inSet( _a.at( edge )        , i    ) &&
        inSet( _a.at( edge ).at( i ), key  ) ) {
        return false;
    }
    return true;
}

/* Set the a value */
void Node::_setAVal( Edge_ptr edge, uint i, condKey key, LogVar aVal ) {

    if( notInSet( _a, edge ) ) {
        _a.emplace( edge, map< uint, map< condKey, LogVar > >() );
    }
    if( notInSet( _a.at( edge ), i ) ) {
        _a.at( edge ).emplace( i, map< condKey, LogVar >() );
    }
    if( inSet( _a.at( edge ).at( i ), key ) ) {
        /* We only want to set this value once! */
        assert( 0 );
    }
    _a.at( edge ).at( i ).emplace( key, aVal );
}

LogVar Node::_getAVal( Edge_ptr edge, uint i, condKey key ) {
    return _a.at( edge ).at( i ).at( key );
}

/* Get or compute the a value */
LogVar Node::_getA( Edge_ptr edge, uint i, const conditioning& cond ) {

    condKey key = _aKey( cond );

    LogVar aVal;

    /* a_n_e( i ) = P( Y\!( e, n ), n_x = i ) */
    if( _needToComputeA( edge, i, key ) ) {
        aVal = _computeA( edge, i, cond );
        _setAVal( edge, i, key, aVal );
    }
    else {
        aVal = _getAVal( edge, i, key );
    }

    return aVal;
}

/* Compute the a value */
LogVar Node::_computeA( Edge_ptr edge, uint i, const conditioning& cond ) {

    LogVar aVal( 1 );

    /* All nodes up from this node */
    aVal *= _getU( i, cond );
    inplace_set_union( _aDependencies, _UDependencies );

    /* All nodes down from this node but not down edge */
    if( this->downEdges.size() > 1 ) {

        for( const Edge_ptr& e : this->downEdges ) {
            if( e == edge ) { continue; }

            /* Down each mate branch */
            aVal *= _getV( e, i, cond );
        }
        inplace_set_union( _aDependencies, _VDependencies );
    }

    return aVal;
}

/* Marginalize out the feedback set latent states so that we get the */
/* regular a values                                                  */
LogVar Node::_getMarginalizedA( Edge_ptr edge, uint i, const conditioning& nodesToKeep, const set< Node_ptr >& feedbackSet ) {

    /* Marginalize out all of the nodes in the feedback set */
    std::vector< Node_ptr > toMarginalizeOut = std::vector< Node_ptr >();
    for( Node_ptr const &node : feedbackSet ) {

        if( notInSet( nodesToKeep, node ) && inSet( _aDependencies, node ) ) {
            toMarginalizeOut.push_back( node );
        }
    }

    LogVar aVal( 0 );

    ConditioningProduct xIter = ConditioningProduct( toMarginalizeOut );
    do {
        conditioning cond = xIter.getConditioning();
        cond.insert( nodesToKeep.begin(), nodesToKeep.end() );

        aVal += _getA( edge, i, cond );
    } while( xIter.next() );

    return aVal;
}

/* ------------------------------------------------------------------------- */

/* Get conditioning key for the b value */
condKey Node::_bKey( const conditioning& cond ) {
    return getKey( cond, _bDependencies );
}

/* Check if we need to compute the b value */
bool Node::_needToComputeB( parentStates X, condKey key ) {

    if( inSet( _b        , X   ) &&
        inSet( _b.at( X ), key ) ) {
        return false;
    }
    return true;
}

/* Set the b value */
void Node::_setBVal( parentStates X, condKey key, LogVar bVal ) {

    if( notInSet( _b, X ) ) {
        _b.emplace( X, map< condKey, LogVar >() );
    }
    if( inSet( _b.at( X ), key ) ) {
        /* We only want to set this value once! */
        assert( 0 );
    }
    _b.at( X ).emplace( key, bVal );
}

LogVar Node::_getBVal( parentStates X, condKey key ) {
    return _b.at( X ).at( key );
}

/* Get or compute the b value */
LogVar Node::_getB( parentStates X, const conditioning& cond ) {

    condKey key = _bKey( cond );

    LogVar bVal;

    /* b_n( X ) = P( n_y, Y\^( n )_y|^( n )_x=X ) */
    if( _needToComputeB( X, key ) ) {
        bVal = _computeB( X, cond );
        _setBVal( X, key, bVal );
    }
    else {
        bVal = _getBVal( X, key );
    }
    return bVal;
}

LogVar Node::_computeB( parentStates X, const conditioning& cond ) {

    LogVar bVal( 0 );

    std::pair< uint, uint > startAndEnd = _getN( this, cond );
    for( uint k = startAndEnd.first; k < startAndEnd.second; ++k ) {

        LogVar prod( 1 );

        /* Prob of this sibling */
        LogVar transProb    = transFunc( X, k );
        LogVar emissionProb = emissionFunc( k );

        prod *= transProb * emissionProb;

        /* Branch down siblings */
        for( Edge_ptr const &edge : this->downEdges ) {
            prod *= _getV( edge, k, cond );
        }

        bVal += prod;
    }
    inplace_set_union( _bDependencies, _VDependencies );

    return bVal;
}

LogVar Node::_getMarginalizedB( parentStates X, const conditioning& nodesToKeep, const set<Node_ptr>& feedbackSet ) {

    std::vector< Node_ptr > toMarginalizeOut = std::vector< Node_ptr >();
    for( Node_ptr const &node : feedbackSet ) {

        if( notInSet( nodesToKeep, node ) && inSet( _bDependencies, node ) ) {
            toMarginalizeOut.push_back( node );
        }
    }

    LogVar bVal( 0 );

    ConditioningProduct xIter = ConditioningProduct( toMarginalizeOut );
    do {
        conditioning cond = xIter.getConditioning();
        cond.insert( nodesToKeep.begin(), nodesToKeep.end() );

        bVal += _getB( X, cond );
    } while( xIter.next() );
    return bVal;
}

/* ------------------------------------------------------------------------- */

/* Get conditioning key for the U value */
condKey Node::_UKey( const conditioning& cond ) {
    return getKey( cond, _UDependencies );
}

/* Check if we need to compute the U value */
bool Node::_needToComputeU( uint i, condKey key ) {

    if( inSet( _U        , i   ) &&
        inSet( _U.at( i ), key ) ) {
        return false;
    }
    return true;
}

/* Set the U value */
void Node::_setUVal( uint i, condKey key, LogVar uVal ) {

    if( notInSet( _U, i ) ) {
        _U.emplace( i, map< condKey, LogVar >() );
    }
    if( inSet( _U.at( i ), key ) ) {
        /* We only want to set this value once! */
        assert( 0 );
    }
    _U.at( i ).emplace( key, uVal );
}

LogVar Node::_getUVal( uint i, condKey key ) {
    return _U.at( i ).at( key );
}

/* Get or compute the U value */
LogVar Node::_getU( uint i, const conditioning& cond ) {

    condKey key = _UKey( cond );

    LogVar uVal;

    /* U_n( i ) = P( n_y,^( n )_y, n_x=i ) */
    if( _needToComputeU( i, key ) ) {
        uVal = _computeU( i, cond );
        _setUVal( i, key, uVal );
    }
    else {
        uVal = _getUVal( i, key );
    }
    return uVal;
}

LogVar Node::_computeU( uint i, const conditioning& cond ) {

    LogVar uVal( 0 );

    if( this->parents.size() == 0 || this->inFeedbackSet ) {

        if( this->inFeedbackSet ) {
            /* Prob if we conditioned on this node */
            uVal = LogVar( ( uint )( i == cond.at( this ) ) );
        }
        else {
            /* Root probability */
            LogVar rootProb     = rootFunc( i );
            LogVar emissionProb = emissionFunc( i );
            uVal = rootProb * emissionProb;
        }
    }
    else {

        std::vector< uint > maxVals = std::vector< uint >();
        std::vector< std::pair< uint, uint > > constValues;

        uint j = 0;
        /* Determine what parent latent states to sum over */
        for( Node_ptr const &_parent : this->parents ) {
            Node_ptr parent = ( Node_ptr )_parent;

            std::pair< uint, uint > startAndEnd = _getN( parent, cond );

            /* If we are conditioning on this node, then it will be */
            /* constant in the loop                                 */
            if( startAndEnd.first == startAndEnd.second ) {
                constValues.push_back( std::make_pair( j, startAndEnd.first ) );
            }
            /* Otherwise add the max value to maxVals so that we can */
            /* iterate over range(N)                                 */
            else {
                maxVals.push_back( startAndEnd.first );
            }
            ++j;
        }

        /* Sum over all of the parent latent states */
        CartesianNodeProduct xIter = CartesianNodeProduct( maxVals, constValues );
        do {
            parentStates X = xIter.getParentStates();

            /* Prob of this node */
            LogVar prod = transFunc( X, i );

            /* Branch out from each parent */
            uint j = 0;
            for( Node_ptr const &_parent : this->parents ) {
                Node_ptr parent = ( Node_ptr )_parent;

                prod *= parent->_getA( this->upEdge, X.at( j ), cond );
                inplace_set_union( _UDependencies, _aDependencies );
                ++j;
            }

            /* Branch out from each sibling */
            for( Node_ptr const &_sibling : this->upEdge->children ) {
                if( _sibling == this ) { continue; }
                Node_ptr sibling = ( Node_ptr )_sibling;

                prod *= sibling->_getB( X, cond );
                inplace_set_union( _UDependencies, _bDependencies );
            }
            uVal += prod;
        } while( xIter.next() );

        uVal *= emissionFunc( i );
    }
    return uVal;
}

LogVar Node::_getMarginalizedU( uint i, const conditioning& nodesToKeep, const set<Node_ptr>& feedbackSet ) {

    std::vector< Node_ptr > toMarginalizeOut = std::vector< Node_ptr >();
    for( Node_ptr const &node : feedbackSet ) {

        if( notInSet( nodesToKeep, node ) && inSet( _UDependencies, node ) ) {
            toMarginalizeOut.push_back( node );
        }
    }

    LogVar uVal( 0 );

    ConditioningProduct xIter = ConditioningProduct( toMarginalizeOut );
    do {
        conditioning cond = xIter.getConditioning();
        cond.insert( nodesToKeep.begin(), nodesToKeep.end() );

        uVal += _getU( i, cond );
    } while( xIter.next() );
    return uVal;
}

/* ------------------------------------------------------------------------- */

/* Get conditioning key for the a value */
condKey Node::_VKey( const conditioning& cond ) {
    return getKey( cond, _VDependencies );
}

/* Check if we need to compute the a value */
bool Node::_needToComputeV( Edge_ptr edge, uint i, condKey key ) {

    if( inSet( _V                   , edge ) &&
        inSet( _V.at( edge )        , i    ) &&
        inSet( _V.at( edge ).at( i ), key  ) ) {
        return false;
    }
    return true;
}

/* Set the a value */
void Node::_setVVal( Edge_ptr edge, uint i, condKey key, LogVar vVal ) {

    if( notInSet( _V, edge ) ) {
        _V.emplace( edge, map< uint, map< condKey, LogVar > >() );
    }
    if( notInSet( _V.at( edge ), i ) ) {
        _V.at( edge ).emplace( i, map< condKey, LogVar >() );
    }
    if( inSet( _V.at( edge ).at( i ), key ) ) {
        /* We only want to set this value once! */
        assert( 0 );
    }
    _V.at( edge ).at( i ).emplace( key, vVal );
}

LogVar Node::_getVVal( Edge_ptr edge, uint i, condKey key ) {
    return _V.at( edge ).at( i ).at( key );
}

/* Get or compute the V value */
LogVar Node::_getV( Edge_ptr edge, uint i, const conditioning& cond ) {

    condKey key = _VKey( cond );

    LogVar vVal;

    /* V_n_e( i ) = P( !( e,n ) | n_x=i ) */
    if( _needToComputeV( edge, i, key ) ) {
        vVal = _computeV( edge, i, cond );
        _setVVal( edge, i, key, vVal );
    }
    else {
        vVal = _getVVal( edge, i, key );
    }
    return vVal;
}

LogVar Node::_computeV( Edge_ptr edge, uint i, const conditioning& cond ) {

    if( inFeedbackSet || this->downEdges.size() == 0 ) {
        return LogVar( 1 );
    }

    uint index = -1;
    std::vector< Node_ptr > mates = std::vector< Node_ptr >();

    std::vector< uint > maxVals = std::vector< uint >();
    std::vector< std::pair< uint, uint > > constValues;

    uint j = 0;
    for( Node_ptr const &_mate : edge->parents ) {
        Node_ptr mate = ( Node_ptr )_mate;

        /* We want to condition on this mate so mark that */
        /* the constant index j has value i               */
        if( mate == this ) {
            constValues.push_back( std::make_pair( j, i ) );
            index = j;
            continue;
        }

        /* Update mates */
        mates.push_back( mate );

        std::pair< uint, uint > startAndEnd = _getN( mate, cond );

        /* If we are conditioning on this node, then it will be */
        /* constant in the loop                                 */
        if( startAndEnd.first == startAndEnd.second ) {
            constValues.push_back( std::make_pair( j, startAndEnd.first ) );
        }
        /* Otherwise add the max value to maxVals so that we can */
        /* iterate over range(N)                                 */
        else {
            maxVals.push_back( startAndEnd.first );
        }

        ++j;
    }

    LogVar vVal = LogVar( 0 );

    CartesianNodeProduct xIter = CartesianNodeProduct( maxVals, constValues );
    do {

        /* X is X_ without mate's latent state value */
        parentStates X_ = xIter.getNonConstStates( index );
        parentStates X = xIter.getParentStates();

        LogVar prod = LogVar( 1 );

        /* Branch out from each mate */
        for( uint j = 0; j < X_.size(); ++j ) {

            Node_ptr mate = mates.at( j );
            prod *= mate->_getA( edge, X_.at( j ), cond );
            inplace_set_union( _VDependencies, mate->_aDependencies );
        }

        /* Branch out from each child */
        for( Node_ptr const &_child : edge->children ) {
            Node_ptr child = ( Node_ptr )_child;

            prod *= child->_getB( X, cond );
            inplace_set_union( _VDependencies, child->_bDependencies );
        }
        vVal += prod;
    } while( xIter.next() );
    return vVal;
}

LogVar Node::_getMarginalizedV( Edge_ptr edge, uint i, const conditioning& nodesToKeep, const set<Node_ptr>& feedbackSet ) {

    std::vector< Node_ptr > toMarginalizeOut = std::vector< Node_ptr >();
    for( Node_ptr const &node : feedbackSet ) {

        if( notInSet( nodesToKeep, node ) && inSet( _VDependencies, node ) ) {
            toMarginalizeOut.push_back( node );
        }
    }

    LogVar vVal(0);

    ConditioningProduct xIter = ConditioningProduct( toMarginalizeOut );
    do {
        conditioning cond = xIter.getConditioning();
        cond.insert( nodesToKeep.begin(), nodesToKeep.end() );

        vVal += _getV( edge, i, cond );
    } while( xIter.next() );
    return vVal;
}

/* ------------------------------------------------------------------------------------ */


LogVar Node::_sortaRootProb( const conditioning& cond ) {
    return _msg->_sortaRootProb( cond );
}

void Node::_accumulateFullJoint( const set<Node_ptr>& feedbackSet ) {

    for( uint i = 0; i < this->N; ++i ) {

        LogVar total( 0 );

        ConditioningProduct xIter = ConditioningProduct( std::vector< Node_ptr >( feedbackSet.begin(), feedbackSet.end() ) );
        do {
            conditioning cond = xIter.getConditioning();

            LogVar prod( 1 );

            prod *= _getU( i, cond );

            for( Edge_ptr const &edge : this->downEdges ) {
                prod *= _getV( edge, i, cond );
            }

            prod *= _sortaRootProb( cond );
            total += prod;
        } while( xIter.next() );

        _fullJoint.at( i ) = total;
    }
}

LogVar Node::getFullJoint( uint i ) {
    return _fullJoint.at( i );
}

void Node::_updateFullJoint( uint i, LogVar val ) {
    if( notInVector( _fullJoint, i ) ) {
        _fullJoint.at( i ) = LogVar( 0 );
    }
    _fullJoint.at( i ) += val;
}

/* ------------------------------------------------------------------------------------ */

void Node::reset() {
    _a         = map< Edge_ptr    , map< uint   , map< condKey, LogVar > > >();
    _b         = map< parentStates, map< condKey, LogVar >, StateHash >     ();
    _U         = map< uint        , map< condKey, LogVar > >                ();
    _V         = map< Edge_ptr    , map< uint   , map< condKey, LogVar > > >();
    _fullJoint = std::vector< LogVar >();
}

LogVar Node::transFunc( const parentStates& X, uint k ) {
    /* Assume that there are only 2 parents */
    return _msg->trans.at( X.at( 0 ) ).at( X.at( 1 ) ).at( k );
}

LogVar Node::emissionFunc( uint i ) {
    return _msg->L.at( i );
}

LogVar Node::rootFunc( uint i ) {
    uint rootIndex = _msg->indexOfRoot( id );
    return _msg->_pi.at( rootIndex ).at( i );
}



