import scipy.stats
import numpy as np
import os
import time
import matplotlib.pyplot as plt
import json
import itertools
from HGMP import MessagePasser

hg = MessagePasser(2)
hg.addNode(1)
hg.addNode(3)
hg.addNode(4)
hg.addNode(5)
hg.addNode(7)

hg.addEdge( parentIds=[1],\
            childIds=[3],\
            edgeId=1 )

hg.addEdge( parentIds=[3],\
            childIds=[4],\
            edgeId=2 )

# hg.addEdge( parentIds=[3,4],\
#             childIds=[5],\
#             edgeId=3 )

# hg.addEdge( parentIds=[1,4,5],\
#             childIds=[7],\
#             edgeId=4 )

# hg.initialize()
