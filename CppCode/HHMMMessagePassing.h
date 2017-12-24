#ifndef __HHMMMESSAGEPASSING_H__
#define __HHMMMESSAGEPASSING_H__

#include "HypergraphBase.h"
#include "LogVar.h"

#define MPNW     MessagePassingNodeWrapper< EmissionType >
#define MPNW_ptr MessagePassingNodeWrapper< EmissionType >*

#define MPHGW     MessagePassingHyperGraphWrapper< EmissionType >
#define MPHGW_ptr MessagePassingHyperGraphWrapper< EmissionType >*

#define conditioning std::unordered_map< MPNW_ptr, uint >

template < class EmissionType >
class MessagePassingNodeWrapper;

template < class EmissionType >
class MessagePassingHyperGraphWrapper;


typedef std::string                          condKey;
// typedef std::unordered_map< MPNW_ptr, uint > conditioning;
typedef std::vector< uint >                  parentStates;


typedef struct {
    /* Hash function taken from                                                       */
    /* https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector */
    std::size_t operator()( std::vector<uint32_t> const& vec ) const {
        std::size_t seed = vec.size();
        for( auto& i : vec ) {
            seed ^= i + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
        }
        return seed;
    }
} stateHash;

template < class EmissionType >
class MessagePassingNodeWrapper : public Node {
private:

    MessagePassingHyperGraphWrapper< EmissionType >* _msg;

    map< Edge_ptr    , map< uint   , map< condKey, LogVar > > > _a;
    map< parentStates, map< condKey, LogVar >, stateHash >      _b;
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

    condKey _UKey            ( const conditioning& cond                                                   );
    bool    _needToComputeU  ( uint i, condKey key                                                        );
    void    _setUVal         ( uint i, condKey key, LogVar uVal                                           );
    LogVar  _getUVal         ( uint i, condKey key                                                        );
    LogVar  _getU            ( uint i, const conditioning& cond                                           );
    LogVar  _computeU        ( uint i, const conditioning& cond                                           );
    LogVar  _getMarginalizedU( uint i, const conditioning& nodesToKeep, const set<MPNW_ptr>& feedbackSet );

    condKey _VKey            ( const conditioning& cond                                                                  );
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

    MessagePassingNodeWrapper( MPNW_ptr node , uint N );
    void    reset            (                        );
    LogVar  getFullJoint     ( uint i                 );
};


template < class EmissionType >
class MessagePassingHyperGraphWrapper : public HyperGraph {
private:

    HyperGraph*            _hyperGraph;
    conditioning           _conditioning;
    map< condKey, LogVar > _sortaRootProbs;
    set< MPNW_ptr >        _sortaRootDeps;
    bool                   _preprocessing;

    LogVar ( *_pi    )( uint               );
    LogVar ( *_L     )( uint               );
    LogVar ( *_trans )( parentStates, uint );

    void   _computeForPreprocessing(                          );
    LogVar _sortaRootProb          ( const conditioning& cond );

public:

    set< MPNW_ptr > feedbackSet;

    void   preprocess                  ( const set< MPNW_ptr >& feedbackSet    );
    LogVar isolatedParentJoint         ( MPNW_ptr node, parentStates X, uint i );
    LogVar probOfParentsProducingNode  ( MPNW_ptr node, parentStates X, uint i );
    void   getStats                    (                                       );
    LogVar probOfAllNodeObservations   (                                       );
    float  logProbOfAllNodeObservations(                                       );
    bool   test                        (                                       );
    void   getCounts                   (                                       );

};

#endif