from HHMMModel import HHMMModelBase

class AutosomalMendelModel( HHMMModelBase ):

    def _transPdfClosure( self, graph ):

        def transPdf( parents, child, X, i ):
            assert parents[ 0 ].sex == 'female'
            return self.transDists[ X[ 0 ] ][ X[ 1 ] ].ipdf( i )

        return transPdf

    def _emissionPdfClosure( self, graph ):
        def emissionPdf( person, i ):
            return self.emissionDists[ i ].ipdf( person.y )

        return emissionPdf

    def _rootPdfClosure( self, graph ):
        roots = sorted( graph.roots )
        graphIndex = self.graphs.index( graph )

        def rootPdf( person, i ):
            personIndex = roots.index( person )
            return self.rootDists[ graphIndex ][ personIndex ].ipdf( i )

        return rootPdf
