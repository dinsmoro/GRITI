"""
This function recovers an original array by the removal of a decimal shift through matrix math.
inputs:
    inArray - input array of any size
    shiftXY - (shiftX, shiftY) decimal shift to undo (e.g., shifted 0.2 ->right, put in 0.2 and it will shift it <-left) - can be any size, any decimal, no integers required
    padVal *UNSUPPORTED* - default 0, can set the value that it was padded with when shifted (or the background mean in a real-world sense)
    FLG_forceSum - default 0, flag that forces the sum to be kept the same if 1
        Commentary: Why isn't the sum always preserved?
        This unshift method involves inverting poorly conditioned matricies.
        In ideal conditions, this will give the original matrix back.
        But noise will heavily distort the returned matrix, and with noise the in-sum and out-sum will no longer be the same.
@author - RD & TS 2019
"""

import numpy as np

def RegridderUnShift(inArray, shiftXY, padVal=0, FLG_forceSum=0):
    
    inSum = np.sum(inArray) # Get the sum of the in array (only used if FLG_forceSum == 1)
    inSize = inArray.shape # Get the size of the in array
    shiftX = shiftXY[0] # Get the X shift
    shiftXpix = np.int64(np.ceil(np.abs(shiftX))) # Number of pixels to shift in X by
    if( shiftXpix == 0 ):
        shiftXdir = np.int64(0) # 0 for zero
    else:
        shiftXdir = np.int64(np.round(shiftX/np.abs(shiftX))) # 1 for positive, -1 for negative
        if( shiftXdir == -1 ): # Fix so negative shift code isn't directly needed
            shiftX = 1 + shiftX # Fix it up so it's the positive mirror
        # END IF
    # END IF
        
    shiftY = shiftXY[1] # Get the Y shift
    # Fix for wrong direction Y shift (was shifting down for +values)
    shiftYpix = np.int64(np.ceil(np.abs(shiftY))) # Number of pixels to shift in X by
    if( shiftYpix == 0 ):
        shiftYdir = np.int64(0) # 0 for zero
    else:
        shiftYdir = np.int64(np.round(shiftY/np.abs(shiftY))) # 1 for positive, -1 for negative
        if( shiftYdir == -1 ): # Fix so negative shift code isn't directly needed
            shiftY = 1 + shiftY # Fix it up so it's the positive mirror
        # END IF
    # END IF
    
    outSize = (int(inSize[0] - shiftYpix) , int(inSize[1] - shiftXpix)) # Makes Numba's nopython happy
        
    # Do the X and Y separately to keep it simple
    # Create the matrix to shift the X
    if( shiftX != 0 ):
        if(shiftX < 0.5 ): # This is where the matrix condition (np.linalg.cond() ) starts to go to infinity (65 at 0.5, inf at 0.8)
            shiftMatX = np.diag( np.tile((1-shiftX) , inSize[1]), k = 0 ) # Diagonal of (1-shiftX)
            shiftMatX = shiftMatX + np.diag( np.tile(shiftX , outSize[1]), k = 1  )  # Offset diagonal of shiftX
    #        shiftMatX = np.linalg.inv(shiftMatX) # Inverse the matrix (square matricies are inversable only)
            # Better inverse by improving the condition of the matrix using matrix maths
            U, S, V = np.linalg.svd(shiftMatX) # Singular value decomposition
            # Note that V is a transpose of MATLAB's
            # Singular value decomposition (SVD) will improve the condition
            R = np.linalg.cholesky(np.diag(S)) # Cholesky decomposition
            # Cholesky decomposition on the S (singular values) of the SVD! Note S = R'*R
            shiftMatX = np.transpose(V)@(np.linalg.inv(R))@np.transpose(np.linalg.inv(R))@np.transpose(U) # Inverse with better condition (no pseudoinverse here)
            shiftMatX = shiftMatX[:,0:outSize[1]]  # Cut off the right column
            inArray = inArray@shiftMatX # Do the un-shift calculation in one shot
        else:
            # Go about it another way - rotate the inArray matrix so it's a  (1-shiftX) shift that is stable (so 0.8 -> 0.2 on a rotated matrix)
            inArray = np.rot90(inArray,k=2)
            shiftMatX = np.diag( np.tile(shiftX , inSize[1]), k = 0 ) # Diagonal of (1-shiftX)
            shiftMatX = shiftMatX + np.diag( np.tile((1-shiftX) , outSize[1]), k = 1  )  # Offset diagonal of shiftX
    #        shiftMatX = np.linalg.inv(shiftMatX) # Inverse the matrix (square matricies are inversable only)
            # Better inverse by improving the condition of the matrix using matrix maths
            U, S, V = np.linalg.svd(shiftMatX) # Singular value decomposition
            # Note that V is a transpose of MATLAB's
            # Singular value decomposition (SVD) will improve the condition
            R = np.linalg.cholesky(np.diag(S)) # Cholesky decomposition
            # Cholesky decomposition on the S (singular values) of the SVD! Note S = R'*R
            shiftMatX = np.transpose(V)@(np.linalg.inv(R))@np.transpose(np.linalg.inv(R))@np.transpose(U) # Inverse with better condition (no pseudoinverse here)
            shiftMatX = shiftMatX[:,0:outSize[1]]  # Cut off the bottom row
                    
            inArray = np.rot90(inArray@shiftMatX,k=2) # Do the un-shift calculation in one shot, rotate back
        # This creates a matrix that undoes the shift previously done
    # else:
    #    pass
    #    shiftMatX = np.diag( np.tile(1 , inSize[1]), k = 0  ) # I matrix so nothing changes, since shiftX == 0
    # END IF
    
    # Create the matrix to shift the Y
    if( shiftY != 0 ):
        
        if( shiftY > 0.5 ): # This is where the matrix condition (np.linalg.cond() ) starts to go to infinity (65 at 0.5, inf at 0.8)
            shiftMatY = np.diag( np.tile(shiftY , inSize[0]), k = 0  ) # Diagonal of shiftY
            shiftMatY = shiftMatY + np.diag( np.tile((1-shiftY) , outSize[0]), k = -1  ) # Offset diagonal of (1-shiftY)
    #        shiftMatY = np.linalg.inv(shiftMatY) # Inverse the matrix (square matricies are inversable only)
            # Better inverse by improving the condition of the matrix using matrix maths
            U, S, V = np.linalg.svd(shiftMatY) # Singular value decomposition
            # Note that V is a transpose of MATLAB's
            # Singular value decomposition (SVD) will improve the condition
            R = np.linalg.cholesky(np.diag(S)) # Cholesky decomposition
            # Cholesky decomposition on the S (singular values) of the SVD! Note S = R'*R
            shiftMatY = np.transpose(V)@(np.linalg.inv(R))@np.transpose(np.linalg.inv(R))@np.transpose(U) # Inverse with better condition (no pseudoinverse here)
            shiftMatY = shiftMatY[0:outSize[0],:]  # Cut off the bottom row
            outArray = shiftMatY@inArray # Do the un-shift calculation in one shot
        else:
            # Go about it another way - rotate the inArray matrix so it's an X shift
            inArray = np.rot90(inArray,k=2)
    #            inArrayRot = np.copy(inArray)
            shiftMatY = np.diag( np.tile((1-shiftY) , inSize[0]), k = 0 ) # Diagonal of (1-shiftX)
            shiftMatY = shiftMatY + np.diag( np.tile(shiftY , outSize[0]), k = -1  )  # Offset diagonal of shiftX
    #        shiftMatY = np.linalg.inv(shiftMatY) # Inverse the matrix (square matricies are inversable only)
            # Better inverse by improving the condition of the matrix using matrix maths
            U, S, V = np.linalg.svd(shiftMatY) # Singular value decomposition
            # Note that V is a transpose of MATLAB's
            # Singular value decomposition (SVD) will improve the condition
            R = np.linalg.cholesky(np.diag(S)) # Cholesky decomposition
            # Cholesky decomposition on the S (singular values) of the SVD! Note S = R'*R
            shiftMatY = np.transpose(V)@(np.linalg.inv(R))@np.transpose(np.linalg.inv(R))@np.transpose(U) # Inverse with better condition (no pseudoinverse here)
            shiftMatY = shiftMatY[0:outSize[0],:]  # Cut off the bottom row
            
            outArray = np.rot90(shiftMatY@inArray,k=2) # Do the un-shift calculation in one shot, rotate back
        # END IF
        
        # This creates a matrix that undoes the shift previously done
    else:
        outArray = inArray # Nothing going on in Y
    #    shiftMatY = np.diag( np.tile(1 , inSize[0]), k = 0  ) # I matrix so nothing changes, since shiftY == 0
    # END IF
    
    # NOTE THAT PADVAL IS UNUSED - ASSUMED TO BE 0. Will need to figure out how to implement if not 0.
    
    if( FLG_forceSum == 1 ):
        outArray = outArray*inSum/np.sum(outArray) # Normalize to the original inArray sum
    # END IF

    # Returns: unshifted array
    return outArray