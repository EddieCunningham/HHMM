from HHMMMessagePasser import HiddenMarkovModelMessagePasser
from HHMMHG import MessagePassingHG
from HHMMNode import NodeForHMM

class PersonNode( NodeForHMM ):
    def __init__( self, y, N, NObs, sex ):
        super( PersonNode, self ).__init__( y, N, NObs )
        self.sex = sex

class PedigreeHGBase( MessagePassingHG ):

    def __init__( self, filename ):
        # Set N to None so that we can have different Ns for different nodes
        super( PedigreeHGBase, self ).__init__( N=None, NodeType=PersonNode )
        self._filename = filename

        self._parentsToEdge = {}
        self._addedNodes = set()

        self.setParentSortKey( lambda x:[ 'female', 'male', 'unknown' ].index( x.sex ) )

    def addNode( self, ID, y, N, NObs, sex ):
        node = super( MessagePassingHG, self ).addNode( ID, y, N, NObs, sex )
        return node

    def calcN( self, person ):
        assert 0, 'Implement this'
        return 4

    def calcNObs( self, person ):
        assert 0, 'Implement this'
        return 2

    def getY( self, person ):
        if( 'diagnoses' in dir( person ) ):
            return int( len( person.diagnoses ) > 0 )
        return 0

    def getSex( self, person ):
        return person.sex

    def checkNode( self, person ):
        if( len( person.parents + person.mateKids ) == 0 ):
            return False

        personId = person.Id

        if( self.hasNode( personId ) ):
            currentNode = super( PedigreeHGBase, self ).getNode( personId )
        else:
            y = self.getY( person )
            N = self.calcN( person )
            NObs = self.calcNObs( person )
            sex = self.getSex( person )
            currentNode = self.addNode( personId, y, N, NObs, sex )

        return currentNode

    def checkEdge( self, parents ):
        sortedParents = tuple( self.parentSort( parents ) )
        if( sortedParents not in self._parentsToEdge ):
            edgeNumber = len( self._parentsToEdge )
            self._parentsToEdge[ sortedParents ] = edgeNumber
            return super( PedigreeHGBase, self ).addEdge( parents, edgeNumber )
        else:
            edgeNumber = self._parentsToEdge[ sortedParents ]
            return super( PedigreeHGBase, self ).getEdge( edgeNumber )

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

        super( PedigreeHGBase, self ).initialize()

class AutosomalPedigree( PedigreeHGBase ):

    def calcN( self, person ):
        return 4

    def calcNObs( self, person ):
        return 2

class XLinkedPedigree( PedigreeHGBase ):

    def calcN( self, person ):
        if( person.sex == 'female' ):
            return 2
        elif( person.sex == 'male' ):
            return 4
        elif( person.sex == 'unknown' ):
            return 6
        else:
            assert 0, 'invalid sex'

    def calcNObs( self, person ):
        return 2

