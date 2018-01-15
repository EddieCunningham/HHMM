from Mendel import AutosomalMendelModel
from pedigreeHGs import allPedigrees
import numpy as np
from mcmcTests import ratioTest, marginalizeTests
from Parameters import autosomalHyperParameters
from PedigreeHypergraph import PedigreeHG

graphs = allPedigrees( nPedigrees=1 )

A_hyper, L_hyper, pi_hyper = autosomalHyperParameters()
model = AutosomalMendelModel( graphs, A_hyper, L_hyper, pi_hyper )
marginalizeTests( model )
log_mean, log_var = model.modelEvidence()
print('Mean: %f, var: %f'%(log_mean, log_var))

# # marginalizeTests( model )
# ratioTest( model )
# ratioTest( model )
# ratioTest( model )
# ratioTest( model )