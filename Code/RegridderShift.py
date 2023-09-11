"""
This function regrids an array by a decimal shift through matrix math and shifting for integers
inputs:
    inArray - input array of any size
    shiftXY - (shiftX, shiftY) decimal shift - can be any size, any decimal, no integers required
    padVal *UNSUPPORTED* - default 0, can set the value it pads with on the edges to a constant number
    FLG_origSize - default 0, 0 means that the outArray will NOT have the original inArray size (as the array spreads out when shifting)
                    1 means that the outArray will have the original inArray size using RegridderZen to readjust back (RegridderZen is flux conserving as well)
@author - RD & TS 2019
"""

import numpy as np
from RegridderZen import RegridderZen

def RegridderShift(inArray, shiftXY, padVal=0, FLG_origSize=0):
    inSize = inArray.shape # Get the size of the in array
    shiftX = shiftXY[0] # Get the X shift
    shiftXpix = np.int64(np.ceil(np.abs(shiftX))) # Number of pixels to shift in X by
    if( shiftXpix == 0 ):
        shiftXdir = np.int64(0) # 0 for zero
    else:
        shiftXdir = np.int64(np.round(shiftX/np.abs(shiftX))) # 1 for positive, -1 for negative
        if( shiftXdir == -1 ): # Fix so negative shift code isn't directly needed
            shiftX = np.remainder(shiftX,1) # Get just the decimal portion (if a shift of like 2.2 occured, it'll do 2.2-2=0.2)
            if( shiftX < 0 ):
                shiftX = 1 + shiftX # Fix it up so it's the positive mirror
            # END IF
        else:
            shiftX = np.remainder(shiftX,1) # Get just the decimal portion (if a shift of like 2.2 occured, it'll do 2.2-2=0.2)
        # END IF
    # END IF
    shiftY = shiftXY[1] # Get the Y shift
    # Fix for wrong direction Y shift (was shifting down for +values)
    # shiftY = -shiftY # Fixes it
    shiftYpix = np.int64(np.ceil(np.abs(shiftY))) # Number of pixels to shift in X by
    if( shiftYpix == 0 ):
        shiftYdir = np.int64(0) # 0 for zero
    else:
        shiftYdir = np.int64(np.round(shiftY/np.abs(shiftY))) # 1 for positive, -1 for negative
        if( shiftYdir == -1 ): # Fix so negative shift code isn't directly needed
            shiftY = np.remainder(shiftY,1) # Get just the decimal portion (if a shift of like 2.2 occured, it'll do 2.2-2=0.2)
            if( shiftY < 0 ):
                shiftY = 1 + shiftY # Fix it up so it's the positive mirror
            # END IF
        else:
            shiftY = np.remainder(shiftY,1) # Get just the decimal portion (if a shift of like 2.2 occured, it'll do 2.2-2=0.2)
        # END IF
    # END IF
    
    outSize = ( int(inSize[0] + np.abs(shiftYdir)) , int(inSize[1] + np.abs(shiftXdir)) ) # Makes Numba's nopython happy
    
    # Create the matrix to shift the X
    if( shiftX != 0 ):
        # This creates a matrix that combines two values in X
        shiftMatX = np.diag( np.tile((1-shiftX) , outSize[1]), k = 0 ) # Diagonal of (1-shiftX)
        shiftMatX = shiftMatX + np.diag( np.tile(shiftX , inSize[1]), k = 1  )  # Offset diagonal of shiftX
        shiftMatX = shiftMatX[0:inSize[1],:]  # Cut off the bottom row
    else:
        shiftMatX = np.diag( np.ones( (inSize[1],) ), k = 0  ) # I matrix so nothing changes, since shiftX == 0
    # END IF
    
    # Create the matrix to shift the Y
    if( shiftY != 0 ):
        # This creates a matrix that combines two values in Y - just how the shift does
        shiftMatY = np.diag( np.tile(shiftY , outSize[0]), k = 0  ) # Diagonal of shiftY
        shiftMatY = shiftMatY + np.diag( np.tile((1-shiftY) , inSize[0]), k = -1  ) # Offset diagonal of (1-shiftY)
        shiftMatY = shiftMatY[:,0:inSize[0]]  # Cut off the right column
    else:
        shiftMatY = np.diag( np.ones( (inSize[0],) ), k = 0  ) #I matrix so nothing changes, since shiftY == 0
    # END IF
    
    # NOTE THAT PADVAL IS UNUSED - ASSUMED TO BE 0. Will need to figure out how to implement if not 0.
    outArray = shiftMatY@inArray@shiftMatX # Do the shift calculation in one shot
    
    # Shift by whole pixels if whole pixel shifts are needed as well
    if( ((shiftX != 0) & ((shiftXpix-1) > 0)) & (shiftXdir == 1) ): # There's a whole pixel shift in the positive direction (->)
        outArray = np.pad(outArray, ( (0,0),((shiftXpix-1),0) ) , 'constant', constant_values=padVal) # Pad with 0's on the edge that's needed (left gets padded with 0's since the whole thing moved right (+) by X pixels)
    elif( ((shiftX != 0) & ((shiftXpix-1) > 0)) & (shiftXdir == -1) ): # There's a whole pixel shift in the negative direction (<-)
        outArray = np.pad(outArray, ( (0,0),(0,(shiftXpix-1)) ) , 'constant', constant_values=padVal) # Pad with 0's on the edge that's needed (right gets padded with 0's since the whole thing moved left (-) by X pixels)
    # END IF
    if( ((shiftY != 0) & ((shiftYpix-1) > 0)) & (shiftYdir == 1) ): # There's a whole pixel shift in the positive direction (->)
        outArray = np.pad(outArray, ( (0,(shiftYpix-1)),(0,0) ) , 'constant', constant_values=padVal) # Pad with 0's on the edge that's needed (down gets padded with 0's since the whole thing moved up (+) by X pixels)
    elif( ((shiftY != 0) & ((shiftYpix-1) > 0)) & (shiftYdir == -1) ): # There's a whole pixel shift in the negative direction (<-)
        outArray = np.pad(outArray, ( ((shiftYpix-1),0),(0,0) ) , 'constant', constant_values=padVal) # Pad with 0's on the edge that's needed (up gets padded with 0's since the whole thing moved down (-) by X pixels)
    # END IF
    #special case for 1 pixel shift, maybe I could do it better but this is what I got :3
    if( ((shiftX == 0) & (shiftXpix > 0)) & (shiftXdir == 1) ): # There's a whole pixel shift in the positive direction (->)
        outArray = np.pad(outArray, ( (0,0),((shiftXpix),0) ) , 'constant', constant_values=padVal) # Pad with 0's on the edge that's needed (left gets padded with 0's since the whole thing moved right (+) by X pixels)
    elif( ((shiftX == 0) & (shiftXpix > 0)) & (shiftXdir == -1) ): # There's a whole pixel shift in the negative direction (<-)
        outArray = np.pad(outArray, ( (0,0),(0,(shiftXpix)) ) , 'constant', constant_values=padVal) # Pad with 0's on the edge that's needed (right gets padded with 0's since the whole thing moved left (-) by X pixels)
    # END IF
    if( ((shiftY == 0) & (shiftYpix > 0)) & (shiftYdir == 1) ): # There's a whole pixel shift in the positive direction (->)
        outArray = np.pad(outArray, ( (0,(shiftYpix)),(0,0) ) , 'constant', constant_values=padVal) # Pad with 0's on the edge that's needed (down gets padded with 0's since the whole thing moved up (+) by X pixels)
    elif( ((shiftY == 0) & (shiftYpix > 0)) & (shiftYdir == -1) ): # There's a whole pixel shift in the negative direction (<-)
        outArray = np.pad(outArray, ( ((shiftYpix),0),(0,0) ) , 'constant', constant_values=padVal) # Pad with 0's on the edge that's needed (up gets padded with 0's since the whole thing moved down (-) by X pixels)
    # END IF
    
    if( FLG_origSize == 1 ): # If original size flag is on, keep the original size intact
        outArray = RegridderZen(outArray, tuple(inSize) ) # Regrid for the original input size
    # END IF

    # Returns: shifted array
    return outArray