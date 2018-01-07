from HHMMMessagePasser import HiddenMarkovModelMessagePasser
from HHMMHG import MessagePassingHG
from HHMMNode import NodeForHMM
from AutosomalDistribution import hyperGraphBaseHyperParameters

def getY( person ):
    if( 'diagnoses' in dir( person ) ):
        return int( len( person.diagnoses ) > 0 )
    return 0

def getSex( person ):
    return person.sex

def calcN( person ):
    return 2

class PersonNode( NodeForHMM ):
    def __init__( self, y, N, sex ):
        super( PersonNode, self ).__init__( y, N )
        self.sex = sex

class PedigreeHG( MessagePassingHG ):

    def __init__( self, filename ):
        # Set N to None so that we can have different Ns for different nodes
        super( PedigreeHG, self ).__init__( N=None, NodeType=PersonNode )
        self._filename = filename

        self._parentsToEdge = {}
        self._addedNodes = set()

        self.setParentSortKey( lambda x:[ 'female', 'male', 'unknown' ].index( x.sex ) )

    def addNode( self, ID, y, N, sex ):
        node = super( MessagePassingHG, self ).addNode( ID, y, N, sex )
        return node

    def checkNode( self, person ):
        if( len( person.parents + person.mateKids ) == 0 ):
            return False

        personId = person.Id

        if( self.hasNode( personId ) ):
            currentNode = super( PedigreeHG, self ).getNode( personId )
        else:
            y = getY( person )
            N = calcN( person )
            sex = getSex( person )
            currentNode = self.addNode( personId, y, N, sex )

        return currentNode

    def checkEdge( self, parents ):
        sortedParents = tuple( self.parentSort( parents ) )
        if( sortedParents not in self._parentsToEdge ):
            edgeNumber = len( self._parentsToEdge )
            self._parentsToEdge[ sortedParents ] = edgeNumber
            return super( PedigreeHG, self ).addEdge( parents, edgeNumber )
        else:
            edgeNumber = self._parentsToEdge[ sortedParents ]
            return super( PedigreeHG, self ).getEdge( edgeNumber )

    def initialize( self, pedigree ):

        for person in pedigree.family:

            if( len( person.parents + person.mateKids ) == 0 ):
                continue

            personNode = self.checkNode( person )
            if( not personNode ): continue

            for mate, children in person.mateKids:

                if( len( mate.parents + mate.mateKids ) == 0 ):
                    mate.parents = person.parents
                    mate.mateKids = [person, children]

                mateNode = self.checkNode( mate )
                if( not mateNode ): continue
                familyEdge = self.checkEdge( [ personNode, mateNode ] )

                for child in children:

                    if( len( child.parents + child.mateKids ) == 0 ):
                        continue

                    childNode = self.checkNode( child )
                    if( not childNode ): continue
                    familyEdge.addChild( childNode )

        super( PedigreeHG, self ).initialize()

    # get rid of this later.  this is specific to bayesian stuff
    def initHyperParams( self, transType, emissionType, priorStrength=1 ):
        self._initialHyperParams = hyperGraphBaseHyperParameters( self, transType, emissionType, priorStrength )
        self._hyperParams = hyperGraphBaseHyperParameters( self, transType, emissionType, priorStrength )
