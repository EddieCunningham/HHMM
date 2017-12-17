#include "HHMMMessagePassing.h"
#include <algorithm>


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

    bool next() {

        for( int i = 0; i < _nodes.size(); ++i ) {
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

    conditioning getConditioning() {

        conditioning cond = conditioning();

        for( int i = 0; i < _nodes.size(); ++i ) {
            cond.insert( {{ _nodes.at( i ), _X.at( i ) }} );
        }
        return cond;
    }

    parentStates getParentStates() {
        return _X;
    }
};

class CartesianNodeProduct {
private:
    parentStates _X;
    std::vector< uint > _maxVals;
public:
    CartesianNodeProduct( const std::vector< uint >& maxVals, uint index = -1, uint constVal = -1 ):
    _maxVals( maxVals ) {
        _X = parentStates( 0, maxVals.size() );
    }

    bool next() {

        for( int i = 0; i < _maxVals.size(); ++i ) {
            uint& val = _X.at( i );

            val += 1;
            if( val == _maxVals.at( i )->N ) {
                val = 0;
            }
            else {
                return true;
            }
        }
        return false;
    }

    parentStates getParentStates() {
        return _X;
    }

    parentStates getExtendedParentStates() {
        if( index == -1 ) {
            assert( 0 );
        }
        parentStates X = parentStates(_X.begin(), _X.end());
        X.insert( X.begin() + index, constVal );
        return X;
    }
};

template < typename T >
void inplace_set_union( set<T>& a, const set<T>& b ) {
    for( T& elt : b ) {
        a.insert( b );
    }
}

std::string keyString( MPNW_ptr node, uint state ) {
    return std::to_string( node->_node->id ) +":"+ std::string( state );
}

condKey getKey( const conditioning& cond, set deps ) {

    std::string tmpString;
    std::vector< std::string > keyBuffer = std::vector< std::string >();

    /* Conditioning only applies to the */
    /* intersection of cond and deps    */
    for( std::pair< const MPNW_ptr, uint > &c : cond ) {

        MPNW_ptr node = c.first;
        uint state    = c.second;

        if( inSet( deps, node ) ) {

            /* Build the key string for this node and */
            /* add it to the key buffer               */
            tmpString = keyString( node, state );
            keyBuffer.push_back( tmpString );
        }
    }

    /* Sort the key buffer */
    std::sort( keyBuffer.begin(), keyBuffer.end() );

    condKey = "";
    for( std::string &keySegment : keyBuffer ) {
        key += keySegment+" ";
    }

    return key;
}

uint MPNW::_getN( MPNW_ptr node , conditioning cond ) {
    if( inSet( cond, node ) ) {
        return cond.at(node);
    }
    return node->N;
}

/* ------------------------------------------------------------------------- */

condKey MPNW::_aKey( const conditioning& cond ) {
    return getKey( cond, _aDependencies );
}

bool MPNW::_needToComputeA( Edge_ptr edge, uint i, condKey key ) {

    if( inSet( _a                   , edge ) &&
        inSet( _a.at( edge )        , i    ) &&
        inSet( _a.at( edge ).at( i ), key  ) ) {
        return false;
    }
    return true;
}

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

LogVar MPNW::_computeA( Edge_ptr edge, uint i, const conditioning& cond ) {

    LogVar aVal( 1 );

    /* All nodes up from this node */
    aVal *= _getU( i, cond );
    inplace_set_union( _aDependencies, _UDependencies );

    /* All nodes down from this node but not down edge */
    if( _node->downEdges.size() > 1 ) {

        for( Edge_ptr& e : _node->downEdges ) {
            if( e == edge ) { continue; }

            /* Down each mate branch */
            aVal *= _getV( e, i, cond )
        }
        inplace_set_union( _aDependencies, _VDependencies );
    }
    return aVal;
}

LogVar MPNW::_getMarginalizedA( Edge_ptr edge, uint i, const conditioning& nodesToKeep, const set<MPNW_ptr>& feedbackSet ) {

    std::vector< MPNW_ptr > toMarginalizeOut = std::vector< MPNW_ptr >();
    for( MPNW_ptr& node : feedbackSet ) {
        if( notInSet( nodesToKeep, node ) && inSet( _aDependencies, node ) ) {
            toMarginalizeOut.push_back( node->N );
        }
    }

    LogVar aVal( 0 );

    ConditioningProduct xIter = ConditioningProduct( toMarginalizeOut );
    while( xIter.next() ) {
        conditioning cond = xIter.getConditioning();
        cond.insert( nodesToKeep.begin(), nodesToKeep.end() );

        aVal += _getA( edge, i, cond );
    }
    return aVal;
}

/* ------------------------------------------------------------------------- */

condKey MPNW::_bKey( const conditioning& cond ) {
    return getKey( cond, _bDependencies );
}

bool MPNW::_needToComputeB( parentStates X, condKey key ) {

    if( inSet( _b        , X   ) &&
        inSet( _b.at( X ), key ) ) {
        return false;
    }
    return true;
}

void MPNW::_setBVal( parentStates X, condKey key, LogVar bVal ) {

    if( notInSet( _b, X ) ) {
        _a.emplace( X, map< condKey, LogVar >() );
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

    for( uint k = 0; k < _getN( this, cond ); ++k ) {

        LogVar prod( 1 );

        /* Prob of this sibling */
        LogVar transProb    = LogVar( ( *_trans )( parentStates, k ) );
        LogVar emissionProb = LogVar( ( *_L ).operator()< EmissionType >( k ) );

        prod *= transProb * emissionProb;

        /* Branch down siblings */
        for( Edge_ptr& e : _node->downEdges ) {
            prod *= _getV( edge, k, cond );
        }

        bVal += prod;
    }
    inplace_set_union( _bDependencies, _VDependencies );

    return bVal;
}

LogVar MPNW::_getMarginalizedB( parentStates X, const conditioning& nodesToKeep, const set<MPNW_ptr>& feedbackSet ) {

    std::vector< MPNW_ptr > toMarginalizeOut = std::vector< MPNW_ptr >();
    for( MPNW_ptr& node : feedbackSet ) {
        if( notInSet( nodesToKeep, node ) && inSet( _bDependencies, node ) ) {
            toMarginalizeOut.push_back( node->N );
        }
    }

    LogVar bVal( 0 );

    ConditioningProduct xIter = ConditioningProduct( toMarginalizeOut );
    while( xIter.next() ) {
        conditioning cond = xIter.getConditioning();
        cond.insert( nodesToKeep.begin(), nodesToKeep.end() );

        bVal += _getB( X, cond );
    }
    return bVal;
}

/* ------------------------------------------------------------------------- */

condKey MPNW::_UKey( const conditioning& cond ) {
    return getKey( cond, _UDependencies );
}

bool MPNW::_needToComputeU( uint i, condKey key ) {

    if( inSet( _U        , i   ) &&
        inSet( _U.at( i ), key ) ) {
        return false;
    }
    return true;
}

void MPNW::_setUVal( uint i, condKey key, LogVar uVal ) {

    if( notInSet( _U, i ) ) {
        _a.emplace( i, map< condKey, LogVar >() );
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

LogVar MPNW::_getU( uint i, const conditioning& cond ) {

    condKey key = _bKey( cond );

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

    if( _parents.size() == 0 || this->inFeedbackSet ) {

        if( this->inFeedbackSet ) {
            /* Prob if we conditioned on this node */
            uVal = LogVar( ( uint )( i == conditioning.at( this ) ) );
        }
        else {
            LogVar rootProb     = LogVar( ( *_pi )( k ) );
            LogVar emissionProb = LogVar( ( *_L ).operator()< EmissionType >( k ) );
            uVal = rootProb * emissionProb;
        }
    }
    else {

        std::vector< uint > maxVals = std::vector< uint >();
        for( MPNW_ptr& parent : _parents ) {
            maxVals.push_back( _getN( parent, cond ) );
        }

        CartesianNodeProduct xIter = CartesianNodeProduct( maxVals );
        while( xIter.next() ) {
            parentStates X = xIter.getParentStates();

            LogVar prod( 1 );

            /* Prob of this node */
            LogVar transProb = LogVar( ( *_trans )( X, i ) );

            /* Branch out from each parent */
            for( MPNW_ptr& parent : _parents ) {
                prod *= parent->_getA( _upEdge, j, cond );
            }

            /* Branch out from each sibling */
            for( MPNW_ptr& sibling : _upEdge->children ) {
                if( sibling == this ) { continue; }

                prod *= sibling->_getB( X, cond );
                inplace_set_union( _UDependencies, _bDependencies );
            }
            uVal += prod;
        }
        uVal *= LogVar( ( *_L ).operator()< EmissionType >( i ) );
    }
    return uVal;
}

LogVar MPNW::_getMarginalizedU( uint i, const conditioning& nodesToKeep, const set<MPNW_ptr>& feedbackSet ) {

    std::vector< MPNW_ptr > toMarginalizeOut = std::vector< MPNW_ptr >();
    for( MPNW_ptr& node : feedbackSet ) {
        if( notInSet( nodesToKeep, node ) && inSet( _UDependencies, node ) ) {
            toMarginalizeOut.push_back( node->N );
        }
    }

    LogVar uVal( 0 );

    ConditioningProduct xIter = ConditioningProduct( toMarginalizeOut );
    while( xIter.next() ) {
        conditioning cond = xIter.getConditioning();
        cond.insert( nodesToKeep.begin(), nodesToKeep.end() );

        uVal += _getU( i, cond );
    }
    return uVal;
}

/* ------------------------------------------------------------------------- */

condKey MPNW::_VKey( const conditioning& cond ) {
    return getKey( cond, _VDependencies );
}

bool MPNW::_needToComputeV( Edge_ptr edge, uint i, condKey key ) {

    if( inSet( _V                   , edge ) &&
        inSet( _V.at( edge )        , i    ) &&
        inSet( _V.at( edge ).at( i ), key  ) ) {
        return false;
    }
    return true;
}

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

    if( inFeedbackSet || _downEdges.size() == 0 ) {
        return LogVar( 1 );
    }

    uint index = -1;
    std::vector< MPNW_ptr > mates = std::vector< MPNW_ptr >();
    std::vector< uint > mateRanges = std::vector< uint >();
    for( int j = 0; j < parents.size(); ++j ) {

        MPNW_ptr mate = parents.at( j );
        if( mate == this ) {
            index = j;
            continue;
        }
        mateRanges.push_back( _getN( mate, cond ) );
    }

    vVal = LogVar( 0 );

    ConditioningProduct xIter = ConditioningProduct( mateRanges, index, i );
    while( xIter.next() ) {
        parentStates X = xIter.getParentStates();
        parentStates X_ = xIter.getExtendedParentStates();

        LogVar prod = LogVar( 1 );

        /* Branch out from each mate */
        for( int j = 0; j < X.size(); ++j ) {

            MPNW_ptr mate = mates.at( j );
            prod *= mate->_getA( edge, j, cond );
            inplace_set_union( _VDependencies, mate->_aDependencies );
        }

        /* Branch out from each child */
        for( MPNW_ptr& child : edge->children ) {

            prod *= child->_getB( X_, cond );
            inplace_set_union( _VDependencies, child->_bDependencies );
        }
        vVal += prod;
    }
    return vVal;
}

LogVar MPNW::_getMarginalizedV( Edge_ptr edge, uint i, const conditioning& nodesToKeep, const set<MPNW_ptr>& feedbackSet ) {

    std::vector< MPNW_ptr > toMarginalizeOut = std::vector< MPNW_ptr >();
    for( MPNW_ptr& node : feedbackSet ) {
        if( notInSet( nodesToKeep, node ) && inSet( _VDependencies, node ) ) {
            toMarginalizeOut.push_back( node->N );
        }
    }

    LogVar vVal(0);

    ConditioningProduct xIter = ConditioningProduct( toMarginalizeOut );
    while( xIter.next() ) {
        conditioning cond = xIter.getConditioning();
        cond.insert( nodesToKeep.begin(), nodesToKeep.end() );

        vVal += _getV( edge, i, cond );
    }
    return vVal;
}

/* ------------------------------------------------------------------------------------ */










