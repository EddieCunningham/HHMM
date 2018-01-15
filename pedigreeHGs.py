from HHMMMessagePasser import HiddenMarkovModelMessagePasser
from HHMMHG import MessagePassingHG
from model import Pedigree
from PedigreeHypergraph import PedigreeHG
import json
import os
import numpy as np

def pedigreeExampleOldJSON( name ):
    pedigreeFolderName = '/Users/Eddie/kec-bot/app/pedigreeDataOLDBUTWORKS/'

    pedigreeNames = [ name ]

    allPedigrees = {}
    for name in pedigreeNames:

        filename = os.path.join( pedigreeFolderName, name+'.json' )
        with open( filename ) as data_file:
            data = json.loads( json.load( data_file ) )
        pedigree = Pedigree( data )
        """ This does nothing,  just haven't gotten around
            to changing how model.py works """
        pedigree.setAffectedFunctions( lambda x: None )

        hg = PedigreeHG( name )
        hg.initialize( pedigree )


        allPedigrees[ name ] = hg

    hg = allPedigrees[ name ]

    return hg

def allPedigrees( nPedigrees=-1 ):

    graphs = []

    goodPedigrees = 0
    badPedigrees = 0


    IPMap = { 'AD' : 0, 'AR' : 1, 'XL' : 2, 'M' : 3 }
    nIPGood = np.zeros( 4 )
    nIP = np.zeros( 4 )

    dataFolder = 'Pedigrees_JSON_Fixed_Label'
    for filename in os.listdir( dataFolder ):

        try:
            filename = os.path.join( dataFolder, filename )

            with open( filename ) as data_file:
                data = json.loads( json.load( data_file ) )

            pedigree = Pedigree( data )
            nIP[ IPMap[ pedigree.inheritancePattern ] ] += 1

            hg = PedigreeHG( filename )
            hg.initialize( pedigree )

            print( 'Done with %s (%s)'%( filename, pedigree.inheritancePattern ) )

            goodPedigrees += 1

            nIPGood[ IPMap[ pedigree.inheritancePattern ] ] += 1

            graphs.append( hg )

            if( nPedigrees != -1 and len( graphs ) >= nPedigrees ):
                break

        except Exception as error:
            print( 'FAILED ON %s.  %s'%( filename, str( error ) ) )

            badPedigrees += 1

        # assert 0

    # print('Good: %d'%goodPedigrees)
    # print('Bad: %d'%badPedigrees)

    # print('AD: %d, AR: %d, XL: %d, M: %d'%( nIP[ 0 ], nIP[ 1 ], nIP[ 2 ], nIP[ 3 ] ) )
    # print('AD: %d, AR: %d, XL: %d, M: %d'%( nIPGood[ 0 ], nIPGood[ 1 ], nIPGood[ 2 ], nIPGood[ 3 ] ) )
    return graphs