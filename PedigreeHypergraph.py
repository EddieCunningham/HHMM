# from CythonCode.HHMMUpDownFast import HiddenMarkovModelMessagePasser,MessagePassingHG
from HHMMUpDown import HiddenMarkovModelMessagePasser,MessagePassingHG
from AutosomalDistribution import *

class PedigreeHG(MessagePassingHG):

    def __init__(self,filename):
        super(PedigreeHG,self).__init__()
        self._filename = filename

        self._parentsToEdge = {}
        self._addedNodes = set()

    def checkNode(self,person):
        if(len(person.parents + person.mateKids) == 0):
            return False
        if(self.hasNode(person.Id)):
            currentNode = super(PedigreeHG,self).getNode(person.Id)
        else:
            currentNode = super(PedigreeHG,self).addNode(person.Id)
        return currentNode

    def checkEdge(self,parents):
        sortedParents = tuple(sorted(parents))
        if(sortedParents not in self._parentsToEdge):
            edgeNumber = len(self._parentsToEdge)
            self._parentsToEdge[sortedParents] = edgeNumber
            return super(PedigreeHG,self).addEdge(parents,edgeNumber)
        else:
            edgeNumber = self._parentsToEdge[sortedParents]
            return super(PedigreeHG,self).getEdge(edgeNumber)

    def initialize(self,pedigree):

        for person in pedigree.family:

            if(len(person.parents + person.mateKids) == 0):
                continue

            personNode = self.checkNode(person)
            if(not personNode): continue

            for mate,children in person.mateKids:

                if(len(mate.parents + mate.mateKids) == 0):
                    mate.parents = person.parents
                    mate.mateKids = [person,children]

                mateNode = self.checkNode(mate)
                if(not mateNode): continue
                familyEdge = self.checkEdge([personNode,mateNode])

                for child in children:

                    if(len(child.parents + child.mateKids) == 0):
                        continue

                    childNode = self.checkNode(child)
                    if(not childNode): continue
                    familyEdge.addChild(childNode)

        # for node in sorted(self._nodes):
        #     print('\n')
        #     print(node)
        #     print(node._upEdge)
        #     print(node._downEdges)
        #     assert len(node._downEdges) > 0 or node._upEdge

        # assert 0

        super(PedigreeHG,self).initialize()

    # get rid of this later.  this is specific to bayesian stuff
    def initHyperParams(self,transType,emissionType,priorStrength=1):
        self._initialHyperParams = hyperGraphBaseHyperParameters(self,transType,emissionType,priorStrength)
        self._hyperParams = hyperGraphBaseHyperParameters(self,transType,emissionType,priorStrength)
