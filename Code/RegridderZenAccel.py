"""
This function regrids an array to a new array size through matrix math
inputs:
    inArray - input array of any size
    yMatrix - y matrix that allows for regridding, produced by RegridderZenMats
    xMatrix - x matrix that allows for regridding, produced by RegridderZenMats
@author - RD & TS 2019
"""

import numpy as np

def RegridderZenAccel(inArray, yMatrix, xMatrix):

    # inSize = inArray.shape # Get the size of the in array
    # xFact = inSize[0]/outSize[0]; # Get the y resize factor
    # yFact = inSize[1]/outSize[1]; # Get the x resize factor
    # outArray = np.zeros(outSize); # Create output array
    
    # # ========== Y Matrix Creation (Uses all x variables, naming flipped) ==========
    # xRange = np.arange(0,inSize[0]+xFact,xFact); # Get a range of pixel slices to make
    # xRangeAdj = np.copy(xRange); # Copy to avoid python variable issues
    # for i in range(0,xRangeAdj.size): # Run through all of the x ranges and adjust as needed
    #     if( xRangeAdj[i] == int(xRangeAdj[i]) ): # Catch true integers
    #         xRangeAdj[i] = xRangeAdj[i] - .1; # Force any int values to floor to the previous integer for counting reasons
    #     # END IF
    # # END FOR i
    
    # if( inSize[0] != outSize[0] ):   
    #     yMatrix = np.zeros( (outSize[0],inSize[0]) ); # Preallocate
    #     for i in range(0,outSize[0]):
    #         xPixelsInvolved = np.arange(int(xRangeAdj[i]), int(xRangeAdj[i+1])+1,1); # Get the pixels involved
    #         xPixelPercent = np.ones(xPixelsInvolved.size); # Get the percentages involved
    #         if( xPixelsInvolved.size == 1 ): # If size is one, then need the difference between the xRange[1] and xRange[0] to cover the divide
    #             xPixelPercent[0] = xRange[i+1] - xRange[i]; # Sets the only value to be the difference between the percentages
    #         else:
    #             if( i == 0 ): # If i == 0 it means that xRange[0] = 0, so special +1 added
    #                 xPixelPercent[0] = np.ceil(xRange[i]) + 1 - xRange[i]; # Sets the first value to be the first percentage
    #             else:
    #                 xPixelPercent[0] = np.ceil(xRange[i]) - xRange[i]; # Sets the first value to be the first percentage
    #             # END IF
    #             xPixelPercent[-1] = xRange[i+1] - int(xRangeAdj[i+1]); # Sets the last value to be the last percentage
    #             # Middle numbers are 1 since full pixels are taken, and if there's only 1 pixel involved the last percentage (one that matters )
    #         # END IF
    #         for j in range(0,xPixelsInvolved.size):
    #             yMatrix[i,xPixelsInvolved[j]] = xPixelPercent[j]; # Put in the pixel percent (the extra magic is using xPixelsInvolved as an index)
    #         # END FOR j
    #     # END FOR i
        
    # else:
    #     yMatrix = np.diag( np.tile(1 , inSize[0]), k = 0  ); # I matrix so nothing changes, since shiftX == 0
    # # END IF
    
    
    # # ========== X Matrix Creation (Uses all y variables, naming flipped) ==========
    # yRange = np.arange(0,inSize[1]+yFact,yFact); # Get a range of pixel slices to make
    # yRangeAdj = np.copy(yRange); # Copy to avoid python variable issues
    # for i in range(0,yRangeAdj.size): # Run through all of the y ranges and adjust as needed
    #     if( yRangeAdj[i] == int(yRangeAdj[i]) ): # Catch true integers
    #         yRangeAdj[i] = yRangeAdj[i] - .1; # Force any int values to floor to the previous integer for counting reasons
    #     # END IF
    # # END FOR i
    
    # if( inSize[1] != outSize[1] ):   
    #     xMatrix = np.zeros( (inSize[1],outSize[1]) ) # Preallocate
    #     for i in range(0,outSize[1]):
    #         yPixelsInvolved = np.arange(int(yRangeAdj[i]), int(yRangeAdj[i+1])+1,1); # Get the pixels involved
    #         yPixelPercent = np.ones(yPixelsInvolved.size); # Get the percentages involved
    #         if( yPixelsInvolved.size == 1 ): # If size is one, then need the difference between the yRange[1] and yRange[0] to cover the divide
    #             yPixelPercent[0] = yRange[i+1] - yRange[i]; # Sets the only value to be the difference between the percentages
    #         else:
    #             if( i == 0 ): # If j == 0 it means that yRange[0] = 0, so special +1 added
    #                 yPixelPercent[0] = np.ceil(yRange[i]) + 1 - yRange[i]; # Sets the first value to be the first percentage
    #             else:
    #                 yPixelPercent[0] = np.ceil(yRange[i]) - yRange[i]; # Sets the first value to be the first percentage
    #             # END IF
    #             yPixelPercent[-1] = yRange[i+1] - int(yRangeAdj[i+1]); # Sets the last value to be the last percentage
    #             # Middle numbers are 1 since full pixels are taken, and if there's only 1 pixel involved the last percentage (one that matters )
    #         # END IF
    #         for j in range(0,yPixelsInvolved.size):
    #             xMatrix[yPixelsInvolved[j],i] = yPixelPercent[j]; # Put in the pixel percent (the extra magic is using yPixelsInvolved as an index)
    #         # END FOR j
    #     # END FOR i
        
    # else:
    #     xMatrix = np.diag( np.tile(1 , inSize[1]), k = 0  ); # I matrix so nothing changes, since shiftX == 0
    # # END IF 
    
    # ========== Resize ========== 
    # Check for the need to convert to sparse matricies for speed (if array is mostly 0's)
    # if( (inSize[0]>1000)|(inSize[1]>1000)|(outSize[0]>1000)|(outSize[1]>1000) ):
    #     # If array large, use sparse arrays - big speed up (but gotta be big to be worth it)
    #     from scipy.sparse import csr_matrix
    #     yMatrix = csr_matrix(yMatrix); # The resize will use sparse matrix math to use these right, 3.5x speed up on large array
    #     xMatrix = csr_matrix(xMatrix);
    #     # Resize
    #     outArray = csr_matrix.dot(csr_matrix.dot(yMatrix, inArray),xMatrix); # yMatrix and xMatrix are sparse matricies - do sparse math to go fast
    # else:
        # Smaoll array and regular, don't use sparse
        # Resize
        # if (np.all(xMatrix == yMatrix) == False):  # Prevents the math if a (size1,size1) -> (size1,size1) was requested (so nothing changes and xMatrix/yMatrix = Identity Matrix)
    outArray = np.matmul(np.matmul(yMatrix, inArray), xMatrix);  # One shot calculate the resize
        # else:
        #     outArray = inArray; #would be same, just set it to be same
        # # END IF (commented out assuming time saving didn't matter much if rarely will have samesize->samesize)
        # outArray = yMatrix@inArray@xMatrix; # One shot calculate the resize (matmul may be faster)
    # END IF
    
    # Returns: Resized array
    return outArray
#END DEF