import time
import itertools
from cycleDetector import identifyCycles
from exampleHG import *
from pedigreeHGs import *

def easyInit( hg ):
    feedbackSet, blockManager = identifyCycles( hg )
    feedbackSetIds = [ node._id for node in feedbackSet ]
    hg.preprocess( feedbackSetIds )
    hg.draw()
    return hg


def HMMMargTest( hg ):

    # Compute the feedback set and preprocess
    hg = easyInit( hg )

    # run the test
    start = time.time()
    hg.marginalizeTest( printStuff=False )
    end = time.time()
    print( 'Test time: '+str( end-start ) )

def marginalizationTests():
    # HMMMargTest( cycleExample1( isHidden=True ) )
    # HMMMargTest( cycleExample3( isHidden=True ) )
    # HMMMargTest( cycleExample5( isHidden=True ) )
    # HMMMargTest( cycleExample5_1( isHidden=True ) )
    # HMMMargTest( cycleExample7( isHidden=True ) )
    # HMMMargTest( cycleExample4( isHidden=True ) )
    # HMMMargTest( cycleExample8( isHidden=True ) )
    # HMMMargTest( cycleExample9( isHidden=True ) )
    # HMMMargTest( cycleExample12( isHidden=True ) )
    HMMMargTest( pedigreeExampleOldJSON( '3818J' ) )
    # HMMMargTest( cycleExample6( isHidden=True ) )

    # These take forever to test!
    # HMMMargTest( cycleExample10( isHidden=True ) )
    # HMMMargTest( cycleExample11( isHidden=True ) )

marginalizationTests()

hg = easyInit( cycleExample5( isHidden=True ) )
hg.resampleGraphStates()