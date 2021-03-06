from Mendel import AutosomalMendelModel
from pedigreeHGs import allPedigrees, allXLPedigrees, allARPedigrees, allADPedigrees
import numpy as np
from mcmcTests import ratioTest, marginalizeTests
from Parameters import autosomalHyperParameters


# from Distributions import Categorical

# k = Categorical( alpha=np.array( [ 1,2,2,1 ] ) )
# print( k.log_likelihood() )
# print( k.logpdf( 1 ) )
# print( k.logpdf( 2 ) )
# print( k.logpdf( 3 ) )
# assert 0



# graphs = allPedigrees( nPedigrees=1 )
# graphs = allXLPedigrees( nPedigrees=1 )
# graphs = allARPedigrees( nPedigrees=1 )
graphs = allADPedigrees( nPedigrees=1 )

A_hyper, L_hyper, pi_hyper = autosomalHyperParameters()
model = AutosomalMendelModel( graphs, A_hyper, L_hyper, pi_hyper )

transPost, emissionPost = model.posteriorEstimate()

print( transPost.mean() )
print( np.sqrt( transPost.variance() ) )

print('\n')

print( emissionPost.mean() )
print( np.sqrt( emissionPost.variance() ) )

assert 0
marginalizeTests( model )
log_mean, log_var = model.modelEvidence( samples=5000 )
print('Mean: %f, var: %f'%(log_mean, log_var))

# # marginalizeTests( model )
# ratioTest( model )
# ratioTest( model )
# ratioTest( model )
# ratioTest( model )
