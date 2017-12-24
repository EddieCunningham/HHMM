import itertools
import random
import operator


class Block(object):

    def __init__(self,obj,objType,bm):
        self._obj = obj
        self._chains = set()
        self._count = -1
        self._objType = objType
        if(objType == 'node'):
            self._id = str('n('+str(obj._id)+')')
        elif(objType == 'edge'):
            self._id = str('e('+str(obj._id)+')')
        else:
            assert 0
        self._bm = bm

    def __repr__(self):
        if(self._objType == 'node'):
            return str(self._id)
        return str(self._id)

    def setCount(self,count):
        self._count = count

    def incCount(self):
        self._count += 1

    def addChain(self,chain):
        self._chains.add(chain)

    def accumulateChains(self):
        # turns all of the chains into a set
        # containing all of the blocks
        self._allAncestors = set()
        for chain in self._chains:
            self._allAncestors |= set(chain._blocks)


    def getAncenstors(self):
        if(self._objType == 'node'):
            if(len([n for n in self._obj._parents if n not in self._bm._ex]) > 0):
                return [self._bm.getBlock(self._obj._upEdge)]
            else:
                return []
        elif(self._objType == 'edge'):
            return [self._bm.getBlock(n) for n in self._obj._parents if n not in self._bm._ex]
        else:
            assert 0
    def getDescendants(self,endNodes=[]):
        if(self._objType == 'node'):
            if(self._obj in endNodes):
                return []
            return [self._bm.getBlock(e) for e in self._obj._downEdges \
                       if len([n for n in e._children if n not in self._bm._ex]) > 0]
        elif(self._objType == 'edge'):
            return [self._bm.getBlock(n) for n in self._obj._children if n not in self._bm._ex]
        else:
            assert 0

    def printAllChains(self):
        print(str(self._id)+' {')
        for chain in self._chains:
            print('\t'),
            chain.printChain()
        print('}')

class BlockManager(object):

    def __init__(self,hg,visitOnce=False,exclude=[]):
        self._hg = hg
        self._ex = exclude
        self._allBlocks = []
        self._mapper = {}

        for n in hg._nodes:
            if(n in exclude):
                continue
            block = Block(n,'node',self)
            self._allBlocks.append(block)
            block.setCount(0)
            self._mapper[n] = block

        for e in hg._edges:
            if(len([n for n in e._parents if n not in exclude]) == 0 or\
               len([n for n in e._children if n not in exclude]) == 0):
                continue
            block = Block(e,'edge',self)
            self._allBlocks.append(block)
            block.setCount(0)
            self._mapper[e] = block

        self._startNodes = self.setCounts(visitOnce)

    def getChainsForBlocks(self,nodes):
        allChains = []
        for node in nodes:
            block = self._mapper[node]
            for chain in block._chains:
                found = True
                for _node in nodes:
                    if(chain.nodeInChain(_node) == False):
                        found = False
                        break
                if(found):
                    allChains.append(chain)
        return allChains

    def getLatentNodes(self,nodes):

        latentNodes = []

        for node in nodes:
            for _node in nodes:
                chains = self.getChainsForBlocks([node,_node])
                for chain in chains:

                    nodeBuffer = []
                    recording = False
                    for block in chain._blocks:
                        if('e' in block._id):
                            continue
                        currentNode = block._obj

                        if(recording):
                            if(currentNode in [node,_node] or currentNode in nodes):
                                for n in nodeBuffer:
                                    if(n not in latentNodes):
                                        latentNodes.append(n)
                                nodeBuffer = []
                            else:
                                assert currentNode not in nodeBuffer
                                nodeBuffer.append(currentNode)
                        else:
                            if(currentNode in [node,_node]):
                                recording = True

        return latentNodes


    def accumulateAllChains(self):
        self._blocksToAncestors = {}

        for block in self._allBlocks:
            block.accumulateChains()
            self._blocksToAncestors[block] = block._allAncestors

    def ancestorCmp(self,node1,node2):
        # returns true if node2 is an ancestor of node1
        if(node1 == node2):
            return False
        block1 = self._mapper[node1]
        block2 = self._mapper[node2]
        return block2 in block1._allAncestors

    def clearCache(self):
        # will erase all of the chains at each block
        for block in self._allBlocks:
            block._chains = set()


    def setCounts(self,visitOnce=False):
        if(visitOnce == False):

            # if a root in to be excluded, need to add its children
            current = self._hg._roots
            notGood = True
            while(notGood):

                notGood = False
                newCurrent = []

                relevantEdges = []
                [relevantEdges.extend(n._downEdges) for n in current]
                relevantEdges = list(set(relevantEdges))
                for e in relevantEdges:
                    # THIS MIGHT INDUCE A BUG. BUT FUCK IT, NOT TRYNA THINK
                    # OF A BETTER WAY AT THE MOMENT.............
                    if(len([p for p in e._parents if p not in self._ex]) == 0):
                        newCurrent.extend(e._children)
                        notGood = True
                    else:
                        for p in e._parents:
                            if(p not in self._ex and p in current):
                                newCurrent.append(p)
                newCurrent = set(list(newCurrent))
                current = newCurrent

            startNodes = list(current)

            if(DEBUG):
                print('Excluding: '+str(self._ex))
                print('Starting chain gathering with nodes: '+str(current))

            for n in current:
                assert n not in self._ex
            # current = [n for n in self._hg._roots if n not in self._ex]
            while(len(current) > 0):
                [self._mapper[c].incCount() for c in current]
                for c in current:
                    for e in c._downEdges:
                        if(len([n for n in e._parents if n not in self._ex]) == 0 or\
                            len([n for n in e._children if n not in self._ex]) == 0):
                            continue
                        self._mapper[e].incCount()
                # [[self._mapper[e].incCount() for e in c._downEdges] for c in current]
                newNodes = []
                for c in current:
                    for e in c._downEdges:
                        for _c in e._children:
                            if(_c not in self._ex):
                                newNodes.append(_c)
                current = newNodes

            return startNodes
        else:
            for node in self._hg._nodes:
                if(node not in self._ex):
                    self._mapper[node].incCount()
            for edge in self._hg._edges:
                [self._mapper[edge].incCount() for n in edge._parents if n not in self._ex]

            return None


    def getBlock(self,obj,chain=None):
        block = self._mapper[obj]
        if(chain):
            block.addChain(chain)
        return block

    def printBlockChains(self):
        for block in self._allBlocks:
            block.printAllChains()
            print('\n')

class Chain(object):

    def __init__(self,blockManager):
        self._bm = blockManager
        self._blocks = []
        self._dec = 0

    def nodeInChain(self,node):
        for b in self._blocks:
            if(node == b._obj):
                return True
        return False

    def copy(self):
        newChain = Chain(self._bm)
        newChain._blocks = list(self._blocks)
        newChain._dec = 0
        return newChain

    def addObj(self,obj):
        block = self._bm.getBlock(obj,self)
        self._blocks.append(block)
        self._dec = 1
        return self

    def addBlock(self,block):
        self._blocks.append(block)
        block.addChain(self)
        self.resetDec()
        return self

    def resetDec(self):
        self._dec = 1

    def decrementHeadCount(self):
        self._blocks[-1]._count -= self._dec
        self._dec = 0

    def getHead(self):
        return self._blocks[-1]

    def printChain(self):
        print(str(self))

    def __str__(self):
        ans = '[ '
        for b in self._blocks:
            ans += str(b._id)+' '
        ans += ']'
        return ans

    def __repr__(self):
        return str(self)

class Cycle(object):
    def __init__(self):
        self._members = set()

    def addHead(self,block):
        self.addMember(block)
        self._head = block

    def addBase(self,block):
        self.addMember(block)
        self._base = block

    def addMember(self,block):
        self._members.add(block)

    def printCycle(self):
        print('h: '+str(self._head._id)),
        print(' - b: '+str(self._base._id)),
        print(' - members: [ '),
        for m in sorted(list(self._members)):
            print(str(m._id)+' '),
        print(']')

    def __hash__(self):
        theHash = hash(str(sorted([str(b)+',' for b in self._members])))
        return theHash


def continueOnBlock(block,chains,bm,endNodes=[],copy=True):
    # chain = mergeChains(block,chains,bm)
    # return splitChain(chain,block.getDescendants(),bm)

    ans = []
    for chain in chains:
        if(copy):
            ans += [chain.copy().addBlock(_block) for _block in block.getDescendants(endNodes)]
        else:
            ans += [Chain(bm).addBlock(_block) for _block in block.getDescendants(endNodes)]
    return ans

def getAllChains(hg,endNodes=[],exclude=[]):

    # returns a new block manager with all of the chains
    bm = BlockManager(hg,exclude=exclude)
    startNodes = bm._startNodes

    chains = [Chain(bm).addObj(obj) for obj in startNodes if obj not in exclude]


    while(len(chains) > 0):

        if(DEBUG):
            print('\n-------- CHAINS SO FAR --------\n')
            for chain in chains:
                chain.printChain()

        [chain.decrementHeadCount() for chain in chains]

        readyBlocks = {}

        newChains = []
        for chain in chains:
            head = chain.getHead()
            if(head._count == 0):
                if(head not in readyBlocks):
                    readyBlocks[head] = set()
                readyBlocks[head].add(chain)
            else:
                newChains.append(chain)

        for block,_chains in readyBlocks.items():
            __chains = continueOnBlock(block,_chains,bm,endNodes)
            newChains.extend(__chains)

        chains = newChains

    return bm

def findCycles(block):

    allCycles = []
    allHashes = set()

    for c1,c2 in itertools.combinations(block._chains,2):

        chain1 = c1._blocks
        chain2 = c2._blocks

        _intersection = set(chain1) & set(chain2)
        # now need to order the intersection nodes.  also if
        # there isn't a node in between intersection nodes,
        # then we don't have a cycle
        intersection = []
        lastIndex = -1
        lastBlock = None
        for i,block in enumerate(chain1):

            if(block in _intersection):

                # check chain2 to see if there is actually a cycle
                chain2Dist = 0
                if(lastBlock != None):
                    for j,block2 in enumerate(chain2[chain2.index(lastBlock):]):
                        if(block2 == block):
                            chain2Dist = j
                            break

                if(lastIndex != -1 and (i - lastIndex + chain2Dist) > 2):
                    if(len(intersection) == 0):
                        intersection.append(lastBlock)
                    intersection.append(block)

                lastBlock = block
                lastIndex = i


        # there is no intersection
        if(len(intersection) <= 1):
            continue

        lastBlock = intersection[0]
        for _block in intersection[1:]:

            head = lastBlock
            base = _block
            currentCycle = Cycle()
            currentCycle.addHead(head)
            currentCycle.addBase(base)

            # now go down each of the chains to add members
            for block in chain1[chain1.index(head)+1:]:
                if(block == base):
                    break
                currentCycle.addMember(block)

            for block in chain2[chain2.index(head)+1:]:
                if(block == base):
                    break
                currentCycle.addMember(block)

            if(currentCycle.__hash__() in allHashes):
                continue
            allCycles.append(currentCycle)
            allHashes.add(currentCycle.__hash__())

            lastBlock = _block

    return allCycles

def identifyCycles(hg,exclude=[],tryToKeepOutOfFeedbackSet=[]):

    bm = getAllChains(hg,exclude=exclude)

    if(DEBUG):
        print('Found all chains')
        bm.printBlockChains()

    # assert 0

    # look at leaves to identify cycles
    allCycles = []


    initialLeaves = hg._leaves
    notGood = True
    while(notGood):

        notGood = False
        newInitialLeaves = []

        relevantEdges = []
        [relevantEdges.append(n._upEdge) for n in initialLeaves if n._upEdge]
        relevantEdges = list(set(relevantEdges))
        for e in relevantEdges:
            if(len([c for c in e._children if c not in exclude]) == 0):
                newInitialLeaves.extend(e._parents)
                notGood = True
            else:
                for c in e._children:
                    if(c not in exclude and c in initialLeaves):
                        newInitialLeaves.append(c)
        newInitialLeaves = set(list(newInitialLeaves))
        initialLeaves = newInitialLeaves

    startLeaves = list(initialLeaves)

    if(DEBUG):
        print('Starting cycle detection with nodes: '+str(startLeaves))

    for leaf in startLeaves:
        if(leaf in exclude):
            continue
        block = bm._mapper[leaf]
        cycles = findCycles(block)
        allCycles.extend(cycles)


    # now identify the feedback vertex set by finding the
    # node that is in the most cycles incrementally
    feedbackSet = set()

    useMostFreqNode = len(tryToKeepOutOfFeedbackSet) == 0


    def getNodesToCutFromSortedNodeCounts(sortedNodeCounts,cutIndex,ttkofs,umfn):

        if(not umfn):

            newNodeCounts = []
            for block,count in sortedNodeCounts:
                if('n' not in block._id):
                    ancestors = block.getAncenstors()
                    descendants = block.getDescendants()
                    if(len([n for n in ancestors if n._obj in ttkofs]) > 0 and \
                       len([n for n in descendants if n._obj in ttkofs]) > 0):
                        continue
                else:
                    if(block._obj in ttkofs):
                        continue

                newNodeCounts.append((block,count))
            sortedNodeCounts = newNodeCounts



        blockToCut = sortedNodeCounts[cutIndex][0]
        if('n' not in blockToCut._id):
            ancestors = blockToCut.getAncenstors()
            descendants = blockToCut.getDescendants()

            if(len([n for n in ancestors if n._obj in ttkofs]) > 0):

                if(len([n for n in descendants if n._obj in ttkofs]) > 0):
                    assert 0

                nodesToCut = descendants
                print('In here 1')

            elif(len([n for n in descendants if n._obj in ttkofs]) > 0):

                if(len([n for n in ancestors if n._obj in ttkofs]) > 0):
                    assert 0
                nodesToCut = ancestors
                print('In here 2')
            else:
                nodesToCut = ancestors if len(ancestors) < len(descendants) else descendants
                print('In here 3')
        else:
            nodesToCut = [blockToCut]


        for n in nodesToCut:
            if(n._obj in ttkofs):
                assert 0

        if(umfn):
            return nodesToCut,sortedNodeCounts
        return nodesToCut,None



    lastNCycles = -1
    while(len(allCycles) > 0):

        # if(lastNCycles == len(allCycles)):
        #     useMostFreqNode = True

        nodeCounts = {}
        for cycle in allCycles:
            if(DEBUG):
                print('---- CYCLE MEMBERS -----')
                print(cycle._members)
                print('---------')
            for block in cycle._members:
                if('n' not in block._id):
                    for child in block.getDescendants():
                        if('n' in child._id):
                            if(child not in nodeCounts):
                                nodeCounts[child] = 0
                            nodeCounts[child] += 1
                    for parent in block.getAncenstors():
                        if('n' in parent._id):
                            if(parent not in nodeCounts):
                                nodeCounts[parent] = 0
                            nodeCounts[parent] += 1

                # its ok to include edges too
                if(block not in nodeCounts):
                    nodeCounts[block] = 0
                nodeCounts[block] += 1

        cutIndex = 0

        sortedNodeCounts = sorted(nodeCounts.items(), key=operator.itemgetter(1), reverse=True)
        nodesToCut,snc = getNodesToCutFromSortedNodeCounts(sortedNodeCounts,cutIndex,tryToKeepOutOfFeedbackSet,useMostFreqNode)
        if(snc): sortedNodeCounts = snc



        if(DEBUG):
            print('\nSORTED NODE COUNTS')
            print(sortedNodeCounts)


        # remove all cycles that contain this node
        newAllCycles = []
        notFound = True
        while(notFound):

            for cycle in allCycles:

                # if none of the nodes to cut are in the cycle, then have to
                # add this cycle to the next loop
                if(len([n for n in nodesToCut if n in cycle._members]) == 0): newAllCycles.append(cycle)
                else: notFound = False

            if(notFound):
                cutIndex += 1

                if(cutIndex >= len(sortedNodeCounts)):
                    # then we have to use some nodes that we didn't want to
                    assert 0, 'FAILED IN THIS FUNCTION'
                    assert useMostFreqNode == False
                    tryToKeepOutOfFeedbackSet.pop()

                    if(len(tryToKeepOutOfFeedbackSet) == 0):
                        useMostFreqNode = True
                        sortedNodeCounts = sorted(nodeCounts.items(), key=operator.itemgetter(1), reverse=True)

                    cutIndex = 0
                    newAllCycles = []
                    continue

                nodesToCut,nothing = getNodesToCutFromSortedNodeCounts(sortedNodeCounts,cutIndex,[],False)

                newAllCycles = []

        for node in nodesToCut:
            feedbackSet.add(node._obj)


        allCycles = newAllCycles

        lastNCycles = len(allCycles)

        if(DEBUG):
            print('Feedback vertex set: [ '),
            for n in feedbackSet:
                print(str(n._id)+' '),
            print(']')


    feedbackSet = list(feedbackSet)

    return feedbackSet,bm

DEBUG=False

# hg = constructGraph2()
# hg.draw()
# identifyCycles(hg)














