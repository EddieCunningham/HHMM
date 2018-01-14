
import csv

pedigreeToIP = {}

with open( 'Complete_pedigree_data.csv', 'r' ) as csvfile:
    reader = csv.DictReader( csvfile )
    for row in reader:
        name = row[ 'Patient ID' ]
        IP   = row[ 'Inheritance Pattern' ]
        pedigreeToIP[ name ] = IP


inFolder = 'Pedigrees_JSON_Original'
outFolder = 'Pedigrees_JSON_Fixed_Label'
for filename in os.listdir( inFolder ):

    actualIP = pedigreeToIP[ filename.replace( '.json', '' ) ]

    inFile = os.path.join( inFolder, filename )

    with open(inFile, 'r') as f:
        data = f.read()

    findStr = '\\"inheritancePattern\\":'
    start = data.index( findStr ) + len( findStr )
    end = data.find( ',', start )
    data = data[ :start ] + '\\"' + actualIP + '\\"' + data[ end: ]

    outFile = os.path.join( outFolder, filename )

    with open(outFile, 'w') as f:
        f.write( data )

assert 0

