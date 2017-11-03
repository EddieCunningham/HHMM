# from HHMMUpDown4 import *
from HHMMUpDown9 import *
# from HHMMUpDown8 import *
# from HHMMUpDown7 import *
# from HHMMUpDown6 import *
# from HHMMUpDown5 import *
# from HHMMLoopy import *
# from HHMMUpDown3 import *
# from HHMMUpDown2 import *
# from HHMMUpDown import *
from AutosomalDistribution import *

class PedigreeHG(MessagePassingHG):

    def __init__(self,filename):
        super(PedigreeHG,self).__init__()
        self._filename = filename

        self._parentsToEdge = {}
        self._addedNodes = set()

    def checkNode(self,person):
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

            personNode = self.checkNode(person)

            for mate,children in person.mateKids:

                mateNode = self.checkNode(mate)
                familyEdge = self.checkEdge([personNode,mateNode])

                for child in children:

                    childNode = self.checkNode(child)
                    familyEdge.addChild(childNode)

        super(PedigreeHG,self).initialize()

    # get rid of this later.  this is specific to bayesian stuff
    def initHyperParams(self,transType,emissionType,priorStrength=1):
        self._initialHyperParams = hyperGraphBaseHyperParameters(self,transType,emissionType,priorStrength)
        self._hyperParams = hyperGraphBaseHyperParameters(self,transType,emissionType,priorStrength)
