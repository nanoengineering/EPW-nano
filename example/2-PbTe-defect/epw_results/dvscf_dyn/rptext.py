import os
import sys
import fileinput

print ("Text to search for:")
textToSearch = input( "> " ) 

print ("Text to replace it with:")
textToReplace = input( "> " )

print ("File to perform Search-Replace on:")

for i in range(16):
	fileToSearch  = "PbTe.dyn_q"+str(i+1)+".xml"
	#fileToSearch = 'D:\dummy1.txt'

	tempFile = open( fileToSearch, 'r+' )

	for line in fileinput.input( fileToSearch ):
	    if textToSearch in line :
	        print('Match Found')
	    else:
	        print('Match Not Found!!')
	    tempFile.write( line.replace( textToSearch, textToReplace ) )
	tempFile.close()


input( '\n\n Press Enter to exit...' )
