#include "HHMM.h"

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