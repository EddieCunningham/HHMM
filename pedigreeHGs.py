from AutosomalDistribution import generic2DParameters, sampleHyperGraphParameters
from HHMMMessagePasser import HiddenMarkovModelMessagePasser
from HHMMHG import MessagePassingHG
from MMUpDown import MarkovModelMessagePasser
from model import Pedigree
from PedigreeHypergraph import PedigreeHG
import json
import os

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

def pedigreeExample( name ):

    pedigreeFolderName = '.'