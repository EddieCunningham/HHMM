import scipy.stats
import numpy as np

import time
from scipy.special import digamma
import matplotlib.pyplot as plt
from pyLogVar import LogVar
import json
import itertools
from model import *
from PedigreeHypergraph import *
from AutosomalDistribution import *

from HHMMUpDown import *



def allPedigrees():
    pedigreeFolderName = '/Users/Eddie/kec-bot/app/pedigreeDataOLDBUTWORKS/'


    # 3239PB, 3818J, 5092AD, 5546EL, 5596IN, 5712CS, 5713BS, 5865MH, 5992VM, 5992VM, 6050MM
    weird = ['5777AH','235TL','3239PB','5092AD','5546EL',\
    '5596IN','5712CS','5713BS','5865MH','5992VM','5992VM','6050MM']

    pedigreeNames = []
    for filename in os.listdir(pedigreeFolderName):
        if('json' in filename and 'test' not in filename):
            if(filename.strip('.json') in weird):
                continue
            pedigreeNames.append(filename.strip('.json'))

    allPedigrees = {}
    for name in pedigreeNames:

        filename = os.path.join(pedigreeFolderName,name+'.json')
        with open(filename) as data_file:
            data = json.loads(json.load(data_file))
        pedigree = Pedigree(data)
        """ This does nothing, just haven't gotten around
            to changing how model.py works """
        pedigree.setAffectedFunctions(lambda x: None)

        hg = PedigreeHG(name)
        hg.initialize(pedigree)
        allPedigrees[name] = hg

    ans = []
    for name,hg in allPedigrees.items():
        hg.initHyperParams('autosome','dominant',1000)
        msg = HiddenMarkovModelMessagePasser(hg,sampleHyperGraphParameters)
        ans.append({'name':name,'hg':hg,'msg':msg})

    return ans

allObjs = allPedigrees()
print('Done loading pedigrees')
start = time.time()

for i in range(len(allObjs)):
    allObjs[i]['msg'].preprocess()

end = time.time()
print('\nPreprocess time: '+str(end-start))
print(str(len(allObjs))+' pedigrees -> '+str((end-start)/float(len(allObjs)))+' avg')

for it in range(10):
    start = time.time()

    for i in range(len(allObjs)):
        allObjs[i]['msg'].getStats()

    end = time.time()
    print('\nBatch test time: '+str(end-start))
    print(str(len(allObjs))+' pedigrees -> '+str((end-start)/float(len(allObjs)))+' avg')

