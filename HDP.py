class CRP():
    # chinese restaurant process
    def __init__( self, alpha ):
        self.alpha = alpha
        self.nTables = 0
        self.nCustomers = 0
        self.tableCounts = []

    def sample( self, size=1 ):
        # Sample from P( θ_i | θ_1,...,θ_i-1, G_0, alpha_0 )

        samples = []
        for i in range( size ):

            probs = np.array( tableCounts + [ self.alpha ] ) / ( self.nCustomers - 1 + self.alpha )
            tableChoice = np.random.choice( len( probs ), 1, p=probs )
            if( tableChoice == len( probs ) - 1 ):
                self.nTables += 1
                self.tableCounts.append( 1 )
            else:
                self.tableCounts[ tableChoice ] += 1
            self.nCustomers += 1
            samples.append( tableChoice )
        return np.array( samples )


class GEM():
    # stick breaking distribution
    def __init__( self ):
        pass

class DP():
    # dirichlet process
    def __init__( self ):
        pass

class HDP():
    # hierarchical dirichlet process
    def __init__( self ):
        pass