#include "HHMMMessagePassing.h"
#include <algorithm>
#include <cassert>
#include <iostream>

/* Replicates
for X in itertools.product(*[range(n.N) for n in nodes]):
    conditioning = {node:i for node,i in zip(nodes,X)}
*/
class ConditioningProduct {
private:
    parentStates _X;
    std::vector< MPNW_ptr > _nodes;
public:
    ConditioningProduct( const std::vector< MPNW_ptr >& nodes ):
    _nodes( nodes ) {
        _X = parentStates( 0, nodes.size() );
    }

    /* Increment the product by 1 */
    bool next() {

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
    conditioning getConditioning() {

        conditioning cond = conditioning();
        for( uint i = 0; i < _nodes.size(); ++i ) {
            cond.emplace( _nodes.at( i ), _X.at( i ) );
        }
        return cond;
    }

    parentStates getParentStates() {
        return _X;
    }
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
    CartesianNodeProduct( const std::vector< uint >& maxVals, const std::vector< std::pair< uint, uint > >& constValues ):
    _maxVals( maxVals ) {

        _X = parentStates( 0, maxVals.size() + constValues.size() );
        for( const std::pair< uint, uint >& val : constValues ) {
            _X.at( val.first ) = val.second;
            _constIndices.push_back( val.first );
        }
    }

    bool next() {

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
    parentStates getNonConstStates( int index = -1 ) {

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

    parentStates getParentStates() {
        return _X;
    }
};

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

/* A key for a node that will be conditioned on */
std::string keyString( const MPNW_ptr node, uint state ) {
    return std::to_string( node->id ) +":"+ std::to_string( state );
}

/* Get the key for condition nodes in cond.  We only need to condition
   on nodes that are in deps
*/
condKey getKey( const conditioning& cond, const set< MPNW_ptr >& deps ) {

    std::string tmpString;
    std::vector< std::string > keyBuffer = std::vector< std::string >();

    /* Conditioning only applies to the */
    /* intersection of cond and deps    */
    for( const std::pair< const MPNW_ptr, uint > &c : cond ) {

        const MPNW_ptr node = c.first;
        uint state          = c.second;

        if( inSet( deps, ( MPNW_ptr )node ) ) {

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

    for( const std::pair< const MPNW_ptr, uint > &c : cond ) {

        const MPNW_ptr node = c.first;
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

/* Get the range of latent states node can have when conditioned
   on cond
*/
std::pair< uint, uint > MPNW::_getN( MPNW_ptr node , const conditioning& cond ) {
    if( inSet( cond, node ) ) {
        return std::make_pair( cond.at(node), cond.at(node) );
    }
    return std::make_pair( 0, node->N) ;
}

/* ------------------------------------------------------------------------- */

/* Get conditioning key for the a value */
condKey MPNW::_aKey( const conditioning& cond ) {
    return getKey( cond, _aDependencies );
}


/* Check if we need to compute the a value */
bool MPNW::_needToComputeA( Edge_ptr edge, uint i, condKey key ) {

    if( inSet( _a                   , edge ) &&
        inSet( _a.at( edge )        , i    ) &&
        inSet( _a.at( edge ).at( i ), key  ) ) {
        return false;
    }
    return true;
}

/* Set the a value */
void MPNW::_setAVal( Edge_ptr edge, uint i, condKey key, LogVar aVal ) {

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

LogVar MPNW::_getAVal( Edge_ptr edge, uint i, condKey key ) {
    return _a.at( edge ).at( i ).at( key );
}

/* Get or compute the a value */
LogVar MPNW::_getA( Edge_ptr edge, uint i, const conditioning& cond ) {

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
LogVar MPNW::_computeA( Edge_ptr edge, uint i, const conditioning& cond ) {

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
LogVar MPNW::_getMarginalizedA( Edge_ptr edge, uint i, const conditioning& nodesToKeep, const set<MPNW_ptr>& feedbackSet ) {

    /* Marginalize out all of the nodes in the feedback set */
    std::vector< MPNW_ptr > toMarginalizeOut = std::vector< MPNW_ptr >();
    for( MPNW_ptr const &node : feedbackSet ) {

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
condKey MPNW::_bKey( const conditioning& cond ) {
    return getKey( cond, _bDependencies );
}

/* Check if we need to compute the b value */
bool MPNW::_needToComputeB( parentStates X, condKey key ) {

    if( inSet( _b        , X   ) &&
        inSet( _b.at( X ), key ) ) {
        return false;
    }
    return true;
}

/* Set the b value */
void MPNW::_setBVal( parentStates X, condKey key, LogVar bVal ) {

    if( notInSet( _b, X ) ) {
        _b.emplace( X, map< condKey, LogVar >() );
    }
    if( inSet( _b.at( X ), key ) ) {
        /* We only want to set this value once! */
        assert( 0 );
    }
    _b.at( X ).emplace( key, bVal );
}

LogVar MPNW::_getBVal( parentStates X, condKey key ) {
    return _b.at( X ).at( key );
}

/* Get or compute the b value */
LogVar MPNW::_getB( parentStates X, const conditioning& cond ) {

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

LogVar MPNW::_computeB( parentStates X, const conditioning& cond ) {

    LogVar bVal( 0 );

    std::pair< uint, uint > startAndEnd = _getN( this, cond );
    for( uint k = startAndEnd.first; k < startAndEnd.second; ++k ) {

        LogVar prod( 1 );

        /* Prob of this sibling */
        LogVar transProb    = LogVar( ( *_trans )( X, k ) );
        LogVar emissionProb = LogVar( ( *_L )( k ) );

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

LogVar MPNW::_getMarginalizedB( parentStates X, const conditioning& nodesToKeep, const set<MPNW_ptr>& feedbackSet ) {

    std::vector< MPNW_ptr > toMarginalizeOut = std::vector< MPNW_ptr >();
    for( MPNW_ptr const &node : feedbackSet ) {

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
condKey MPNW::_UKey( const conditioning& cond ) {
    return getKey( cond, _UDependencies );
}

/* Check if we need to compute the U value */
bool MPNW::_needToComputeU( uint i, condKey key ) {

    if( inSet( _U        , i   ) &&
        inSet( _U.at( i ), key ) ) {
        return false;
    }
    return true;
}

/* Set the U value */
void MPNW::_setUVal( uint i, condKey key, LogVar uVal ) {

    if( notInSet( _U, i ) ) {
        _U.emplace( i, map< condKey, LogVar >() );
    }
    if( inSet( _U.at( i ), key ) ) {
        /* We only want to set this value once! */
        assert( 0 );
    }
    _U.at( i ).emplace( key, uVal );
}

LogVar MPNW::_getUVal( uint i, condKey key ) {
    return _U.at( i ).at( key );
}

/* Get or compute the U value */
LogVar MPNW::_getU( uint i, const conditioning& cond ) {

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

LogVar MPNW::_computeU( uint i, const conditioning& cond ) {

    LogVar uVal( 0 );

    if( this->parents.size() == 0 || this->inFeedbackSet ) {

        if( this->inFeedbackSet ) {
            /* Prob if we conditioned on this node */
            uVal = LogVar( ( uint )( i == cond.at( this ) ) );
        }
        else {
            /* Root probability */
            LogVar rootProb     = LogVar( ( *_pi )( i ) );
            LogVar emissionProb = LogVar( ( *_L )( i ) );
            uVal = rootProb * emissionProb;
        }
    }
    else {

        std::vector< uint > maxVals = std::vector< uint >();
        std::vector< std::pair< uint, uint > > constValues;

        uint j = 0;
        /* Determine what parent latent states to sum over */
        for( Node_ptr const &_parent : this->parents ) {
            MPNW_ptr parent = ( MPNW_ptr )_parent;

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
            LogVar prod = LogVar( ( *_trans )( X, i ) );

            /* Branch out from each parent */
            uint j = 0;
            for( Node_ptr const &_parent : this->parents ) {
                MPNW_ptr parent = ( MPNW_ptr )_parent;

                prod *= parent->_getA( this->upEdge, X.at( j ), cond );
                inplace_set_union( _UDependencies, _aDependencies );
                ++j;
            }

            /* Branch out from each sibling */
            for( Node_ptr const &_sibling : this->upEdge->children ) {
                if( _sibling == this ) { continue; }
                MPNW_ptr sibling = ( MPNW_ptr )_sibling;

                prod *= sibling->_getB( X, cond );
                inplace_set_union( _UDependencies, _bDependencies );
            }
            uVal += prod;
        } while( xIter.next() );

        uVal *= LogVar( ( *_L )( i ) );
    }
    return uVal;
}

LogVar MPNW::_getMarginalizedU( uint i, const conditioning& nodesToKeep, const set<MPNW_ptr>& feedbackSet ) {

    std::vector< MPNW_ptr > toMarginalizeOut = std::vector< MPNW_ptr >();
    for( MPNW_ptr const &node : feedbackSet ) {

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
condKey MPNW::_VKey( const conditioning& cond ) {
    return getKey( cond, _VDependencies );
}

/* Check if we need to compute the a value */
bool MPNW::_needToComputeV( Edge_ptr edge, uint i, condKey key ) {

    if( inSet( _V                   , edge ) &&
        inSet( _V.at( edge )        , i    ) &&
        inSet( _V.at( edge ).at( i ), key  ) ) {
        return false;
    }
    return true;
}

/* Set the a value */
void MPNW::_setVVal( Edge_ptr edge, uint i, condKey key, LogVar vVal ) {

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

LogVar MPNW::_getVVal( Edge_ptr edge, uint i, condKey key ) {
    return _V.at( edge ).at( i ).at( key );
}

/* Get or compute the V value */
LogVar MPNW::_getV( Edge_ptr edge, uint i, const conditioning& cond ) {

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

LogVar MPNW::_computeV( Edge_ptr edge, uint i, const conditioning& cond ) {

    if( inFeedbackSet || this->downEdges.size() == 0 ) {
        return LogVar( 1 );
    }

    uint index = -1;
    std::vector< MPNW_ptr > mates = std::vector< MPNW_ptr >();

    std::vector< uint > maxVals = std::vector< uint >();
    std::vector< std::pair< uint, uint > > constValues;

    uint j = 0;
    for( Node_ptr const &_mate : edge->parents ) {
        MPNW_ptr mate = ( MPNW_ptr )_mate;

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

            MPNW_ptr mate = mates.at( j );
            prod *= mate->_getA( edge, X_.at( j ), cond );
            inplace_set_union( _VDependencies, mate->_aDependencies );
        }

        /* Branch out from each child */
        for( Node_ptr const &_child : edge->children ) {
            MPNW_ptr child = ( MPNW_ptr )_child;

            prod *= child->_getB( X, cond );
            inplace_set_union( _VDependencies, child->_bDependencies );
        }
        vVal += prod;
    } while( xIter.next() );
    return vVal;
}

LogVar MPNW::_getMarginalizedV( Edge_ptr edge, uint i, const conditioning& nodesToKeep, const set<MPNW_ptr>& feedbackSet ) {

    std::vector< MPNW_ptr > toMarginalizeOut = std::vector< MPNW_ptr >();
    for( MPNW_ptr const &node : feedbackSet ) {

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


LogVar MPNW::_sortaRootProb( const conditioning& cond ) {
    return _msg->_sortaRootProb( cond );
}

void MPNW::_accumulateFullJoint( const set<MPNW_ptr>& feedbackSet ) {

    for( uint i = 0; i < this->N; ++i ) {

        LogVar total( 0 );

        ConditioningProduct xIter = ConditioningProduct( std::vector< MPNW_ptr >( feedbackSet.begin(), feedbackSet.end() ) );
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

LogVar MPNW::getFullJoint( uint i ) {
    return _fullJoint.at( i );
}

void MPNW::_updateFullJoint( uint i, LogVar val ) {
    if( notInVector( _fullJoint, i ) ) {
        _fullJoint.at( i ) = LogVar( 0 );
    }
    _fullJoint.at( i ) += val;
}

/* ------------------------------------------------------------------------------------ */

MPNW::MessagePassingNodeWrapper():
    _aDependencies(),
    _bDependencies(),
    _UDependencies(),
    _VDependencies() {
    reset();
}

void MPNW::reset() {
    _a         = map< Edge_ptr    , map< uint   , map< condKey, LogVar > > >();
    _b         = map< parentStates, map< condKey, LogVar >, StateHash >     ();
    _U         = map< uint        , map< condKey, LogVar > >                ();
    _V         = map< Edge_ptr    , map< uint   , map< condKey, LogVar > > >();
    _fullJoint = std::vector< LogVar >();
}
/* ------------------------------------------------------------------------------------ */

LogVar MPHGW::_sortaRootProb( const conditioning& cond ) {

    condKey key = getKey( cond );

    if( inSet( _sortaRootProbs, key ) ) {
        return _sortaRootProbs.at( key );
    }

    LogVar prod( 1 );

    if( _sortaRootDeps.size() > 0 ) {

        for( MPNW_ptr const &sortaRoot : _sortaRootDeps ) {

            uint j = cond.at( sortaRoot );
            LogVar innerProd = LogVar( ( *_L )( j ) );

            if( sortaRoot->parents.size() == 0 ) {
                innerProd *= LogVar( ( *_pi )( j ) );
            }
            else {
                parentStates X = parentStates();
                for( const std::pair< const MPNW_ptr, uint > &c : cond ) {
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

void MPHGW::_computeForPreprocessing() {

    for( const Node_ptr& _node : nodes ) {
        MPNW_ptr node = ( MPNW_ptr )_node;

        if( node->inFeedbackSet ) {
            continue;
        }

        ConditioningProduct xIter = ConditioningProduct( std::vector< MPNW_ptr >( feedbackSet.begin(), feedbackSet.end() ) );
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
        MPNW_ptr node = ( MPNW_ptr )_node;

        if( node->inFeedbackSet ) { continue; }

        node->_accumulateFullJoint( feedbackSet );
    }

    MPNW_ptr aLeaf = ( MPNW_ptr )( *( leaves.begin() ) );

    _sortaRootProbs = map< condKey, LogVar >();

    ConditioningProduct xIter = ConditioningProduct( std::vector< MPNW_ptr >( feedbackSet.begin(), feedbackSet.end() ) );

    do {
        conditioning cond = xIter.getConditioning();
        parentStates X = xIter.getParentStates();

        for( uint i = 0; i < aLeaf->N; ++i ) {
            LogVar val = aLeaf->_getU( i, cond );
            LogVar srProb = _sortaRootProb( cond );

            uint j = 0;
            for( MPNW_ptr const &node : feedbackSet ) {

                uint x = X.at( j );
                node->_updateFullJoint( x, val * srProb );
                ++j;
            }
        }
    } while( xIter.next() );

    for( MPNW_ptr const &node : feedbackSet ) {
        LogVar total( 0 );
        for( uint i = 0; i < node->N; ++i ) {
            total += node->_fullJoint.at( i );
        }
    }
}

/* ------------------------------------------------------------------------------------ */

void MPHGW::preprocess( const set< MPNW_ptr >& feedbackSet ) {
    this->feedbackSet = feedbackSet;

    for( Node_ptr const &_node : nodes ) {
        MPNW_ptr node = ( MPNW_ptr )_node;
        node->_msg = this;
    }

    _sortaRootDeps = set< MPNW_ptr >();
    for( MPNW_ptr const &node : feedbackSet ) {
        if( node->parents.size() == 0 ) {
            _sortaRootDeps.insert( node );
            continue;
        }

        bool notSRP = false;
        for( Node_ptr const &_parent : node->parents ) {
            MPNW_ptr parent = ( MPNW_ptr )_parent;

            if( notInSet( feedbackSet, parent ) ) {
                notSRP = true;
                break;
            }
        }
        if( notSRP ) {
            continue;
        }
        for( Node_ptr const &_sibling : node->upEdge->children ) {
            MPNW_ptr sibling = ( MPNW_ptr )_sibling;

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
        MPNW_ptr node = ( MPNW_ptr )_node;

        if( inSet( feedbackSet, node ) ) {
            node->inFeedbackSet = true;
        }
    }

    _preprocessing = true;
    _computeForPreprocessing();
    _preprocessing = false;
}

void MPHGW::getCounts() {

    for( Node_ptr const &_node : nodes ) {
        MPNW_ptr node = ( MPNW_ptr )_node;
        node->reset();
    }

    // make this a true message passing
    // algorithm later!
    preprocess( feedbackSet );
}

LogVar MPHGW::isolatedParentJoint( MPNW_ptr node, parentStates X, uint i ) {

    MPNW_ptr aLeaf = ( MPNW_ptr )( *( leaves.begin() ) );

    LogVar jointProb( 0 );

    std::vector< MPNW_ptr > toMarginalizeOut = std::vector< MPNW_ptr >();
    for( MPNW_ptr const &_node : feedbackSet ) {
        if( node != _node && notInSet( node->parents, _node ) ) {
            toMarginalizeOut.push_back( node );
        }
    }

    conditioning parentCond = conditioning();
    uint j = 0;
    for( Node_ptr const &_parent : node->parents ) {
        MPNW_ptr parent = ( MPNW_ptr )_parent;

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

LogVar MPHGW::probOfParentsProducingNode( MPNW_ptr node, parentStates X, uint i ) {

    if( inSet( _sortaRootDeps, node ) ) {
        return isolatedParentJoint( node, X, i );
    }

    std::vector< uint > fbsNotInParents = std::vector< uint >();
    std::vector< std::pair< uint, uint > > constValues = std::vector< std::pair< uint, uint > >();

    uint k = 0;
    for( MPNW_ptr const&fbsNode : feedbackSet ) {
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
        for( MPNW_ptr const& fbsNode : feedbackSet ) {
            cond.emplace( fbsNode, _X.at( k ) );
            ++k;
        }

        LogVar prod( 1 );
        for( Edge_ptr const &edge : node->downEdges ) {
            prod *= node->_getMarginalizedV( edge, i, cond, feedbackSet );
        }

        for( Node_ptr const &_sibling : node->upEdge->children ) {
            MPNW_ptr sibling = ( MPNW_ptr )_sibling;
            if( sibling == node ) {
                continue;
            }
            prod *= sibling->_getMarginalizedB( X, cond, feedbackSet );
        }

        uint j = 0;
        for( Node_ptr const &_parent : node->parents ) {
            MPNW_ptr parent = ( MPNW_ptr )_parent;

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

void MPHGW::getStats() {

    getCounts();

    for( Node_ptr const &_node : nodes ) {
        MPNW_ptr node = ( MPNW_ptr )_node;

        for( uint i = 0; i < node->N; ++i ) {
        }
    }
}

LogVar MPHGW::probOfAllNodeObservations() {

    MPNW_ptr aLeaf = ( MPNW_ptr )( *( leaves.begin() ) );

    LogVar total( 0 );
    for( uint i = 0; i < aLeaf->N; ++i ) {
        total += aLeaf->getFullJoint( i );
    }
    return total;
}

float MPHGW::logProbOfAllNodeObservations() {

    return probOfAllNodeObservations().logValue();
}

void MPHGW::preprocessWithInt( const std::vector< uint >& feedbackSet ) {

    set< MPNW_ptr > fbs = set< MPNW_ptr >();
    for( uint _id : feedbackSet ) {
        MPNW_ptr node = ( MPNW_ptr )getNode( _id );
        fbs.insert( node );
    }
    preprocess( fbs );
}

float MPHGW::logFullJoint( uint nodeId, uint i) {
    MPNW_ptr node = ( MPNW_ptr )getNode( nodeId );
    return node->getFullJoint( i ).logValue();
}

void MPHGW::addNodeId( uint id, uint y ) {

    MPNW_ptr node = ( MPNW_ptr )HyperGraph::addNode( id );
    node->N = N;
    node->y = y;
    node->_msg = this;
}

void MPHGW::addEdgeId( const std::vector< uint > & parents, const std::vector< uint > & children, uint id ) {

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
