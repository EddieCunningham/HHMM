#ifndef __HHMMMESSAGEPASSING_H__
#define __HHMMMESSAGEPASSING_H__

#include "HypergraphBase.h"
#include "LogVar.h"

#define MPNW     MessagePassingNodeWrapper< EmissionType >
#define MPNW_ptr MessagePassingNodeWrapper< EmissionType >*

#define MPHGW     MessagePassingHyperGraphWrapper< EmissionType >
#define MPHGW_ptr MessagePassingHyperGraphWrapper< EmissionType >*

typedef std::string                          condKey;
typedef std::unordered_map< MPNW_ptr, uint > conditioning;
typedef std::vector< uint >                  parentStates;

class MessagePassingNodeWrapper;
class MessagePassingHyperGraphWrapper;

template < typename EmissionType >
class MessagePassingNodeWrapper : public Node {
private:

    MessagePassingHyperGraphWrapper* _msg;

    map< Edge_ptr    , map< uint   , map< condKey, LogVar > > > _a;
    map< parentStates, map< condKey, LogVar > >                 _b;
    map< uint        , map< condKey, LogVar > >                 _U;
    map< Edge_ptr    , map< uint   , map< condKey, LogVar > > > _V;

    set< MPNW_ptr > _aDependencies;
    set< MPNW_ptr > _bDependencies;
    set< MPNW_ptr > _UDependencies;
    set< MPNW_ptr > _VDependencies;

    std::vector< LogVar > _fullJoint;

    LogVar ( *_pi    )( uint               );
    LogVar ( *_L     )( uint               );
    LogVar ( *_trans )( parentStates, uint );

    condKey _aKey            ( const conditioning& cond                                                                  );
    bool    _needToComputeA  ( Edge_ptr edge, uint i, condKey key                                                        );
    void    _setAVal         ( Edge_ptr edge, uint i, condKey key, LogVar aVal                                           );
    LogVar  _getAVal         ( Edge_ptr edge, uint i, condKey key                                                        );
    LogVar  _getA            ( Edge_ptr edge, uint i, const conditioning& cond                                           );
    LogVar  _computeA        ( Edge_ptr edge, uint i, const conditioning& cond                                           );
    LogVar  _getMarginalizedA( Edge_ptr edge, uint i, const conditioning& nodesToKeep, const set<MPNW_ptr>& feedbackSet );

    condKey _bKey            ( const conditioning& cond                                                           );
    bool    _needToComputeB  ( parentStates X, condKey key                                                        );
    void    _setBVal         ( parentStates X, condKey key, LogVar bVal                                           );
    LogVar  _getBVal         ( parentStates X, condKey key                                                        );
    LogVar  _getB            ( parentStates X, const conditioning& cond                                           );
    LogVar  _computeB        ( parentStates X, const conditioning& cond                                           );
    LogVar  _getMarginalizedB( parentStates X, const conditioning& nodesToKeep, const set<MPNW_ptr>& feedbackSet );

    condKey _uKey            ( const conditioning& cond                                                   );
    bool    _needToComputeU  ( uint i, condKey key                                                        );
    void    _setUVal         ( uint i, condKey key, LogVar uVal                                           );
    LogVar  _getUVal         ( uint i, condKey key                                                        );
    LogVar  _getU            ( uint i, const conditioning& cond                                           );
    LogVar  _computeU        ( uint i, const conditioning& cond                                           );
    LogVar  _getMarginalizedU( uint i, const conditioning& nodesToKeep, const set<MPNW_ptr>& feedbackSet );

    condKey _vKey            ( const conditioning& cond                                                                  );
    bool    _needToComputeV  ( Edge_ptr edge, uint i, condKey key                                                        );
    void    _setVVal         ( Edge_ptr edge, uint i, condKey key, LogVar vVal                                           );
    LogVar  _getVVal         ( Edge_ptr edge, uint i, condKey key                                                        );
    LogVar  _getV            ( Edge_ptr edge, uint i, const conditioning& cond                                           );
    LogVar  _computeV        ( Edge_ptr edge, uint i, const conditioning& cond                                           );
    LogVar  _getMarginalizedV( Edge_ptr edge, uint i, const conditioning& nodesToKeep, const set<MPNW_ptr>& feedbackSet );

    LogVar _sortaRootProb      ( const conditioning& cond         );
    void   _accumulateFullJoint( const set<MPNW_ptr>& feedbackSet );
    void   _updateFullJoint    ( uint i, LogVar val               );

    std::pair< uint, uint > _getN( MPNW_ptr node , const conditioning& cond );

public:


    uint         N;
    EmissionType y;
    bool         inFeedbackSet;

    MessagePassingNodeWrapper( Node_ptr node , uint N );
    void    reset            (                        );
    LogVar  getFullJoint     ( uint i                 );
};


template < typename EmissionType >
class MessagePassingHyperGraphWrapper : public HyperGraph {
private:

    HyperGraph*            _hyperGraph;
    conditioning           _conditioning;
    map< condKey, LogVar > _sortaRootProbs;
    set< MPNW_ptr >        _sortaRootDeps;
    bool                   _preprocessing;

    void   _computeForPreprocessing(                          );
    LogVar _sortaRootProb          ( const conditioning& cond );

public:

    set< MPNW_ptr > feedbackSet;

    void   preprocess                (                                       );
    LogVar isolatedParentJoint       ( Node_ptr node, parentStates X, uint i );
    LogVar probOfParentsProducingNode( Node_ptr node, parentStates X, uint i );
    void   getStats                  (                                       );
    LogVar probOfAllNodeObservations (                                       );
    bool   test                      (                                       );
    void   getCounts                 (                                       );

};

#endif