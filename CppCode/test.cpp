#include <iostream>
#include "HHMMMessagePassing.h"

using namespace std;

int main() {

    MessagePassingHyperGraphWrapper mp = MessagePassingHyperGraphWrapper( 2 );

    mp.addNodeId( 1, 0 );
     mp.addNodeId( 3, 0 );
//     mp.addNodeId( 4, 0 );
    // mp.addNodeId( 5, 0 );
    // mp.addNodeId( 7, 0 );

     mp.addEdgeId( vector< uint >( { 1 } ),
                   vector< uint >( { 3 } ),
                   1 );

    // mp.addEdgeId( vector< uint >( { 3 } ),
    //               vector< uint >( { 4 } ),
    //               2 );

    // mp.addEdgeId( vector< uint >( { 3, 4 } ),
    //               vector< uint >( { 5 } ),
    //               3 );

    // mp.addEdgeId( vector< uint >( { 1, 4, 5 } ),
    //               vector< uint >( { 7 } ),
    //               4 );

    return 0;
}
