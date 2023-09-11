#Inserts whatever data comes in
#indexesToInsert and insertStuff has to be the same size (insertStuff can also be a scalar)
#output size is data.size + indexesToInsert.size
#this only does 1D stuff don't 2D me thx

import numpy as np

def subfun_insert( data, indexesToInsert, insertStuff):
    dataOut = np.empty(data.size+indexesToInsert.size,dtype=data.dtype); #preallocate
    
    #----- Insert stuff -----
    indexesForData = np.isin(np.arange(0,dataOut.size,step=1), indexesToInsert, assume_unique=True,invert=True); #get data indexes
    dataOut[indexesForData] = data; #insert that stuff
    dataOut[indexesToInsert] = insertStuff; #insert that stuff
    
    return dataOut
#END DEF