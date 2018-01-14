from Mendel import AutosomalMendelModel
from pedigreeHGs import *
import numpy as np
from exampleHG import *
from mcmcTests import ratioTest, marginalizeTests



from PedigreeHypergraph import PedigreeHG

graphs = []

goodPedigrees = 0
badPedigrees = 0

dataFolder = 'Pedigrees_JSON'
for filename in os.listdir( dataFolder ):

    try:
        filename = os.path.join( dataFolder, filename )

        with open( filename ) as data_file:
            data = json.loads( json.load( data_file ) )

        pedigree = Pedigree( data )

        hg = PedigreeHG( filename )
        hg.initialize( pedigree )

        print( 'Done with %s'%filename )

        goodPedigrees += 1

        graphs.append( hg )

    except Exception as error:
        print( 'FAILED ON %s.  %s'%( filename, str( error ) ) )

        badPedigrees += 1

    # assert 0

print('Good: %d'%goodPedigrees)
print('Bad: %d'%badPedigrees)





A_hyper = np.array([
        [
            [ 1.0, 0.5, 0.5, 0.0 ],
            [ 0.5, .25, .25, 0.0 ],
            [ 0.5, .25, .25, 0.0 ],
            [ 0.0, 0.0, 0.0, 0.0 ]
        ],
        [
            [ 0.0, 0.5, 0.5, 1.0 ],
            [ 0.0, .25, .25, 0.5 ],
            [ 0.0, .25, .25, 0.5 ],
            [ 0.0, 0.0, 0.0, 0.0 ]
        ],
        [
            [ 0.0, 0.0, 0.0, 0.0 ],
            [ 0.5, .25, .25, 0.0 ],
            [ 0.5, .25, .25, 0.0 ],
            [ 1.0, 0.5, 0.5, 0.0 ]
        ],
        [
            [ 0.0, 0.0, 0.0, 0.0 ],
            [ 0.0, .25, .25, 0.5 ],
            [ 0.0, .25, .25, 0.5 ],
            [ 0.0, 0.5, 0.5, 1.0 ]
        ]
    ])

A_hyper = np.array( [ [ [ A_hyper[ k ][ i ][ j ] for k in range( A_hyper.shape[ 2 ] ) ] \
                                                 for j in range( A_hyper.shape[ 1 ] ) ] \
                                                 for i in range( A_hyper.shape[ 0 ] ) ] )

L_hyper = np.array( [ [1.,0.],
                      [1.,0.],
                      [1.,0.],
                      [0.,1.] ] )

A_hyper = 1. + ( 5 * A_hyper )
L_hyper = 1. + ( 5 * L_hyper )
pi_hyper = np.ones( 4 )

model = AutosomalMendelModel( graphs, A_hyper, L_hyper, pi_hyper )
model.resample()
model.resample()
model.resample()

# # marginalizeTests( model )
# ratioTest( model )
# ratioTest( model )
# ratioTest( model )
# ratioTest( model )