# Fits a figure to the window size (be sure to set window size before calling this so it fits to the window size you want)
# bufferPix is the amount of pixels of white space you want around the figure, autocalculated at 0.21% of the widest figure size if not defined

# interactive = True is for interactive plots
# interactive = False is for ioff() plots that never show up and are saved to be viewed (or explictly called)
# interactive is automatically determined if not explicitly defined

# fast = False uses the slower version that uses fine adjustments (slower b/c matplotlib has to re-render after every adjustment)
# fast = True uses the faster version that does a few big adjustments (it drops the fine loop adjustments)
# default is fast as the fine adjustment loops can add several seconds to plot creation

import numpy as np
import matplotlib.pyplot as plt

def figFitter(fig, bufferPix = None, interactive = None, fast = True, returnOffsets = False):
    
    if( interactive == None ): #if interactive is none, find out if it is interactive or not
        interactive = plt.isinteractive(); #record if interactive is true or false
    #END IF
    
    fig.canvas.draw(); #key for all instances
    if( interactive == True ): #only needs fig.canvas.draw() for interactive, flush and show and pause aren't req
        fig.canvas.flush_events(); #only needed for interactive plotting
        plt.show(); #req to make plot show up
        plt.pause(0.01); #this is oddly highly required for plotting in an IDE interactively (seems 0.01 is lowest I can go)
    #END IF
    figSize = np.flip(np.int64(fig.get_size_inches()*fig.dpi)); #pixels, fig size
    if( bufferPix == None ): #if bufferPix is none, estimate a good bufferPix
        if( np.max(figSize) < 3000 ):
            bufferPix = np.int64(np.round(np.max(figSize)*0.0021)); #pixels, automagically estimate the buffer pixel amount
        else:
            #if it's a lotta pixels, increase the buffer b/c it gets weird otherwise I'm not sure
            bufferPix = np.int64(np.round(np.max(figSize)*0.0041)); #pixels, automagically estimate the buffer pixel amount
        #END IF
    #END IF
    buffGood = False; #final logical check
    buffGoodArray = np.zeros(4, dtype=bool); #used to track all of the logicals
    # offsets = np.float64(np.array( (0,1,0,1) )); #the offsets to apply with subplots_adjust
    offsets = np.float64(np.array( (fig.subplotpars.left, fig.subplotpars.right, fig.subplotpars.bottom, fig.subplotpars.top) )); #get the current offsets
    # fig.subplots_adjust(left = offsets[0], right = offsets[1], bottom = offsets[2], top = offsets[3]); #start by setting the padding to nothing
    offsetLatch = np.array( ((0,1,0,1),(0.25,0.75,0.25,0.75)) ); #the minimum/maximum allowed offset before giving up
    offsetAdj = offsetAdj = np.array( ((np.max(figSize)*0.0005/figSize[1],np.max(figSize)*0.0005/figSize[1],np.max(figSize)*0.0005/figSize[0],np.max(figSize)*0.0005/figSize[0])) ); #increment to adjust by
    reqVals = np.array( (figSize[0]*bufferPix*255*3, figSize[0]*bufferPix*255*3, figSize[1]*bufferPix*255*3, figSize[1]*bufferPix*255*3) ); #the reference pixel sums to reach
    reqValsBuff = np.array( (figSize[0]*(bufferPix+1)*255*3, figSize[0]*(bufferPix+1)*255*3, figSize[1]*(bufferPix+1)*255*3, figSize[1]*(bufferPix+1)*255*3) ); #the reference pixel sums to reach
    tryLim = 4; #number of tries for big expand
    #----- Initial Big Expand -----
    #mostly it should start from a smol picture in the center of a big white zone, this makes big moves for that
    rgb_data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep=''); #get the rgb data in a numpy array, based on https://stackoverflow.com/questions/40433211/how-can-i-get-the-pixel-colors-in-matplotlib @ Jean-Sébastien
    rgb_data = rgb_data.reshape((int(len(rgb_data)/3), 3)).reshape(np.hstack((figSize,3))); #get the rgb data in a numpy array, based on https://stackoverflow.com/questions/40433211/how-can-i-get-the-pixel-colors-in-matplotlib @ Jean-Sébastien
    tryOut = 1;
    booster = np.array( (1,1,1,1) ,dtype=np.float64); #prep
    boosterVal = 1.5; #amount to increase buffer size by each try
    while( (buffGood == False) & (tryOut <= tryLim) ):
        #----- One side at a time -----
        #--- Left check ---
        if( np.sum(rgb_data[:,0:(bufferPix+1),:]) == reqValsBuff[0] ):
            pixLim = np.int64(np.round(offsets[0]*figSize[1])); #get a pixel limit to check to (the pixel limit is the offset converted to pixels)
            pixSums = np.flip(np.sum(np.sum(rgb_data[:,0:pixLim,:],axis=0),axis=1)); #get the pixel row sums
            pixFind = np.where(pixSums != figSize[0]*255*3)[0]; #find the pix
            if( pixFind.size != 0 ): #only adjust if it's not all white space, if it is just remove that I guess?
                pixLim = pixFind[-1] + 1 + bufferPix; #index where text starts, adjust by buffer
            else:
                pixLim = 1 + bufferPix; #set to 1+ buffer if it's all whitespace then?
            #END IF
            offsets[0] = pixLim/figSize[1]; #set the offset
            buffGoodArray[0] = False;
        elif( (np.sum(rgb_data[-(bufferPix+1):,:,:]) != reqValsBuff[2]) & (tryOut == 1) ):
            #if on 1st loop through it is already clipping off, move it in
            offsets[0] = offsetLatch[1,0]; #set the offset to the offset latch
            buffGoodArray[0] = False;
        elif( np.all(np.sum(rgb_data[:,0:(bufferPix+1),:],axis=1) != figSize[0]*255*3) == True ):
            booster[0] = booster[0]*boosterVal; #boost 
            pixLim = np.int64(np.round(offsets[0]*figSize[1])); #get a pixel limit to check to (the pixel limit is the offset converted to pixels)
            pixSums = np.flip(np.sum(np.sum(rgb_data[:,0:pixLim,:],axis=0),axis=1)); #get the pixel row sums
            pixFind = np.where(pixSums != figSize[0]*255*3)[0]; #find the pix
            if( pixFind.size != 0 ): #only adjust if it's not all white space, if it is just remove that I guess?
                pixLim = pixFind[-1] + 1 + bufferPix*booster[0]; #index where text starts, adjust by buffer
            else:
                pixLim = 1 + bufferPix*booster[0]; #set to 1+ buffer if it's all whitespace then?
            #END IF
            offsets[0] = pixLim/figSize[1]; #set the offset
            buffGoodArray[0] = False;
        else:
            buffGoodArray[0] = True;
        #END IF
        #--- Bottom check ---
        if( np.sum(rgb_data[-(bufferPix+1):,:,:]) == reqValsBuff[2] ):
            pixLim = np.int64(np.round(offsets[2]*figSize[0])); #get a pixel limit to check to (the pixel limit is the offset converted to pixels)
            pixSums = np.sum(np.sum(rgb_data[-pixLim:,:,:],axis=1),axis=1); #get the pixel row sums
            pixFind = np.where(pixSums != figSize[1]*255*3)[0]; #find the pix
            if( pixFind.size != 0 ): #only adjust if it's not all white space, if it is just remove that I guess?
                pixLim = pixFind[-1] + 1 + bufferPix; #index where text starts, adjust by buffer
            else:
                pixLim = 1 + bufferPix; #set to 1+ buffer if it's all whitespace then?
            #END IF
            offsets[2] = pixLim/figSize[0]; #set the offset
            buffGoodArray[2] = False;
        elif( (np.sum(rgb_data[-(bufferPix+1):,:,:]) != reqValsBuff[2]) & (tryOut == 1) ):
            #if on 1st loop through it is already clipping off, move it in
            offsets[2] = offsetLatch[1,2]; #set the offset to the offset latch
            buffGoodArray[2] = False;
        elif( np.all(np.sum(rgb_data[-(bufferPix+1):,:,:],axis=0) != figSize[1]*255*3) == True ):
            booster[2] = booster[2]*boosterVal; #boost 
            pixLim = np.int64(np.round(offsets[2]*figSize[0])); #get a pixel limit to check to (the pixel limit is the offset converted to pixels)
            pixSums = np.sum(np.sum(rgb_data[-pixLim:,:,:],axis=1),axis=1); #get the pixel row sums
            pixFind = np.where(pixSums != figSize[1]*255*3)[0]; #find the pix
            if( pixFind.size != 0 ): #only adjust if it's not all white space, if it is just remove that I guess?
                pixLim = pixFind[-1] + 1 + bufferPix*booster[2]; #index where text starts, adjust by buffer
            else:
                pixLim = 1 + bufferPix*booster[2]; #set to 1+ buffer if it's all whitespace then?
            #END IF
            offsets[2] = pixLim/figSize[0]; #set the offset
            buffGoodArray[2] = False;
        else:
            buffGoodArray[2] = True;
        #END IF
        # Refresh for opposite sides
        fig.subplots_adjust(left = offsets[0], right = offsets[1], bottom = offsets[2], top = offsets[3]); #sets padding to small numbers for minimal white space
        fig.canvas.draw(); #key for all instances
        if( interactive == True ): #only needs fig.canvas.draw() for interactive, flush isn't req
            fig.canvas.flush_events(); #only needed for interactive plotting
        #END IF
        rgb_data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep=''); #get the rgb data in a numpy array, based on https://stackoverflow.com/questions/40433211/how-can-i-get-the-pixel-colors-in-matplotlib @ Jean-Sébastien
        rgb_data = rgb_data.reshape((int(len(rgb_data)/3), 3)).reshape(np.hstack((figSize,3))); #get the rgb data in a numpy array, based on https://stackoverflow.com/questions/40433211/how-can-i-get-the-pixel-colors-in-matplotlib @ Jean-Sébastien
        
        #----- One side at a time -----
        #--- Right check ---
        if( np.sum(rgb_data[:,-(bufferPix+1):,:]) == reqValsBuff[1] ):
            pixLim = np.int64(np.round((1-offsets[1])*figSize[1])); #get a pixel limit to check to
            pixSums = np.sum(np.sum(rgb_data[:,-pixLim:,:],axis=0),axis=1); #get the pixel row sums
            pixFind = np.where(pixSums != figSize[0]*255*3)[0]; #find the pix
            if( pixFind.size != 0 ): #only adjust if it's not all white space, if it is just remove that I guess?
                pixLim = pixFind[-1] + 1 + bufferPix; #index where text starts, adjust by buffer
            else:
                pixLim = 1 + bufferPix; #set to 1+ buffer if it's all whitespace then?
            #END IF
            offsets[1] = 1 - pixLim/figSize[1]; #set the offset
            buffGoodArray[1] = False;
        elif( (np.sum(rgb_data[:,-(bufferPix+1):,:]) != reqValsBuff[1]) & (tryOut == 1) ):
            #if on 1st loop through it is already clipping off, move it in
            offsets[1] = offsetLatch[1,1]; #set the offset to the offset latch
            buffGoodArray[1] = False;
        elif( np.all(np.sum(rgb_data[:,-(bufferPix+1):,:],axis=1) != figSize[0]*255*3) == True ):
            booster[1] = booster[1]*boosterVal; #boost 
            pixLim = np.int64(np.round((1-offsets[1])*figSize[1])); #get a pixel limit to check to
            pixSums = np.sum(np.sum(rgb_data[:,-pixLim:,:],axis=0),axis=1); #get the pixel row sums
            pixFind = np.where(pixSums != figSize[0]*255*3)[0]; #find the pix
            if( pixFind.size != 0 ): #only adjust if it's not all white space, if it is just remove that I guess?
                pixLim = pixFind[-1] + 1 + bufferPix*booster[1]; #index where text starts, adjust by buffer
            else:
                pixLim = 1 + bufferPix*booster[1]; #set to 1+ buffer if it's all whitespace then?
            #END IF
            offsets[1] = 1 - pixLim/figSize[1]; #set the offset
            buffGoodArray[1] = False;
        else:
            buffGoodArray[1] = True;
        #END IF
        #--- Top check ---
        if( np.sum(rgb_data[0:(bufferPix+1),:,:]) == reqValsBuff[3] ):
            pixLim = np.int64(np.round((1-offsets[3])*figSize[0])); #get a pixel limit to check to
            pixSums = np.flip(np.sum(np.sum(rgb_data[0:pixLim,:,:],axis=1),axis=1)); #get the pixel row sums
            pixFind = np.where(pixSums != figSize[1]*255*3)[0]; #find the pix
            if( pixFind.size != 0 ): #only adjust if it's not all white space, if it is just remove that I guess?
                pixLim = pixFind[-1] + 1 + bufferPix; #index where text starts, adjust by buffer
            else:
                pixLim = 1 + bufferPix; #set to 1+ buffer if it's all whitespace then?
            #END IF
            offsets[3] = 1 - pixLim/figSize[0]; #set the offset
            buffGoodArray[3] = False;
        elif( (np.sum(rgb_data[0:(bufferPix+1),:,:]) != reqValsBuff[3]) & (tryOut == 1) ):
            #if on 1st loop through it is already clipping off, move it in
            offsets[3] = offsetLatch[1,3]; #set the offset to the offset latch
            buffGoodArray[3] = False;
        elif( np.all(np.sum(rgb_data[0:(bufferPix+1),:,:],axis=0) == figSize[1]*255*3) == True ):
            booster[3] = booster[3]*boosterVal; #boost 
            pixLim = np.int64(np.round((1-offsets[3])*figSize[0])); #get a pixel limit to check to
            pixSums = np.flip(np.sum(np.sum(rgb_data[0:pixLim,:,:],axis=1),axis=1)); #get the pixel row sums
            pixFind = np.where(pixSums != figSize[1]*255*3)[0]; #find the pix
            if( pixFind.size != 0 ): #only adjust if it's not all white space, if it is just remove that I guess?
                pixLim = pixFind[-1] + 1 + bufferPix*booster[3]; #index where text starts, adjust by buffer
            else:
                pixLim = 1 + bufferPix*booster[3]; #set to 1+ buffer if it's all whitespace then?
            #END IF
            offsets[3] = 1 - pixLim/figSize[0]; #set the offset
            buffGoodArray[3] = False;
        else:
            buffGoodArray[3] = True;
        #END IF
        # Refresh for opposite sides
        fig.subplots_adjust(left = offsets[0], right = offsets[1], bottom = offsets[2], top = offsets[3]); #sets padding to small numbers for minimal white space
        fig.canvas.draw(); #key for all instances
        if( interactive == True ): #only needs fig.canvas.draw() for interactive, flush isn't req
            fig.canvas.flush_events(); #only needed for interactive plotting
        #END IF
        buffGood = np.all(buffGoodArray); #they all have to approve
        rgb_data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep=''); #get the rgb data in a numpy array, based on https://stackoverflow.com/questions/40433211/how-can-i-get-the-pixel-colors-in-matplotlib @ Jean-Sébastien
        rgb_data = rgb_data.reshape((int(len(rgb_data)/3), 3)).reshape(np.hstack((figSize,3))); #get the rgb data in a numpy array, based on https://stackoverflow.com/questions/40433211/how-can-i-get-the-pixel-colors-in-matplotlib @ Jean-Sébastien
        tryOut += 1; #increment
    #END WHILE
    if( fast == False ): #only do the fine adjustments if fast is off
        if( (np.sum(rgb_data[:,0:bufferPix,:]) != reqVals[0]) | (np.sum(rgb_data[:,-bufferPix:,:]) != reqVals[1]) | (np.sum(rgb_data[-bufferPix:,:,:]) != reqVals[2]) | (np.sum(rgb_data[0:bufferPix,:,:]) != reqVals[3]) ):
            buffGood = False; #set to false to run the condenser
            buffGoodArray = np.zeros(4, dtype=bool); #set to false to run the condenser
        #END IF
        #----- Expander -----
        rgb_data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep=''); #get the rgb data in a numpy array, based on https://stackoverflow.com/questions/40433211/how-can-i-get-the-pixel-colors-in-matplotlib @ Jean-Sébastien
        rgb_data = rgb_data.reshape((int(len(rgb_data)/3), 3)).reshape(np.hstack((figSize,3))); #get the rgb data in a numpy array, based on https://stackoverflow.com/questions/40433211/how-can-i-get-the-pixel-colors-in-matplotlib @ Jean-Sébastien
        while( buffGood == False ):
            #--- Left check ---
            if(  (np.sum(rgb_data[:,0:bufferPix,:]) == reqVals[0]) & (offsets[0] > offsetLatch[0,0]) & (offsets[0] < offsetLatch[1,0]) ):
                offsets[0] -= offsetAdj[0]; #adjust the offset to try to get only whitespace there
            else:
                buffGoodArray[0] = True;
            #END IF
            #--- Right check ---
            if(  (np.sum(rgb_data[:,-bufferPix:,:]) == reqVals[1]) & (offsets[1] < offsetLatch[0,1]) & (offsets[1] > offsetLatch[1,1]) ):
                offsets[1] += offsetAdj[1]; #adjust the offset to try to get only whitespace there
            else:
                buffGoodArray[1] = True;
            #END IF
            #--- Bottom check ---
            if(  (np.sum(rgb_data[-bufferPix:,:,:]) == reqVals[2]) & (offsets[2] > offsetLatch[0,2]) & (offsets[2] < offsetLatch[1,2]) ):
                offsets[2] -= offsetAdj[2]; #adjust the offset to try to get only whitespace there 
            else:
                buffGoodArray[2] = True;
            #END IF
            #--- Top check ---
            if( (np.sum(rgb_data[0:bufferPix,:,:]) == reqVals[3]) & (offsets[3] < offsetLatch[0,3]) & (offsets[3] > offsetLatch[1,3]) ):
                offsets[3] += offsetAdj[3]; #adjust the offset to try to get only whitespace there
            else:
                buffGoodArray[3] = True;
            #END IF
    
            fig.subplots_adjust(left = offsets[0], right = offsets[1], bottom = offsets[2], top = offsets[3]); #sets padding to small numbers for minimal white space
            fig.canvas.draw(); #key for all instances
            if( interactive == True ): #only needs fig.canvas.draw() for interactive, flush isn't req
                fig.canvas.flush_events(); #only needed for interactive plotting
            #END IF
            buffGood = np.all(buffGoodArray); #they all have to approve
            rgb_data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep=''); #get the rgb data in a numpy array, based on https://stackoverflow.com/questions/40433211/how-can-i-get-the-pixel-colors-in-matplotlib @ Jean-Sébastien
            rgb_data = rgb_data.reshape((int(len(rgb_data)/3), 3)).reshape(np.hstack((figSize,3))); #get the rgb data in a numpy array, based on https://stackoverflow.com/questions/40433211/how-can-i-get-the-pixel-colors-in-matplotlib @ Jean-Sébastien
        #END WHILE
        rgb_data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep=''); #get the rgb data in a numpy array, based on https://stackoverflow.com/questions/40433211/how-can-i-get-the-pixel-colors-in-matplotlib @ Jean-Sébastien
        rgb_data = rgb_data.reshape((int(len(rgb_data)/3), 3)).reshape(np.hstack((figSize,3))); #get the rgb data in a numpy array, based on https://stackoverflow.com/questions/40433211/how-can-i-get-the-pixel-colors-in-matplotlib @ Jean-Sébastien
        if( (np.sum(rgb_data[:,0:bufferPix,:]) != reqVals[0]) | (np.sum(rgb_data[:,-bufferPix:,:]) != reqVals[1]) | (np.sum(rgb_data[-bufferPix:,:,:]) != reqVals[2]) | (np.sum(rgb_data[0:bufferPix,:,:]) != reqVals[3]) ):
            buffGood = False; #set to false to run the condenser
            buffGoodArray = np.zeros(4, dtype=bool); #set to false to run the condenser
        #END IF
        #----- Condenser -----
        rgb_data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep=''); #get the rgb data in a numpy array, based on https://stackoverflow.com/questions/40433211/how-can-i-get-the-pixel-colors-in-matplotlib @ Jean-Sébastien
        rgb_data = rgb_data.reshape((int(len(rgb_data)/3), 3)).reshape(np.hstack((figSize,3))); #get the rgb data in a numpy array, based on https://stackoverflow.com/questions/40433211/how-can-i-get-the-pixel-colors-in-matplotlib @ Jean-Sébastien 
        while( buffGood == False ):
            rgb_data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep=''); #get the rgb data in a numpy array, based on https://stackoverflow.com/questions/40433211/how-can-i-get-the-pixel-colors-in-matplotlib @ Jean-Sébastien
            rgb_data = rgb_data.reshape((int(len(rgb_data)/3), 3)).reshape(np.hstack((figSize,3))); #get the rgb data in a numpy array, based on https://stackoverflow.com/questions/40433211/how-can-i-get-the-pixel-colors-in-matplotlib @ Jean-Sébastien
            
            #--- Left check ---
            if(  (np.sum(rgb_data[:,0:bufferPix,:]) != reqVals[0]) & (offsets[0] > offsetLatch[0,0]) & (offsets[0] < offsetLatch[1,0]) ):
                offsets[0] += offsetAdj[0]; #adjust the offset to try to get only whitespace there      
            else:
                buffGoodArray[0] = True;
            #END IF
            #--- Right check ---
            if(  (np.sum(rgb_data[:,-bufferPix:,:]) != reqVals[1]) & (offsets[1] < offsetLatch[0,1]) & (offsets[1] > offsetLatch[1,1]) ):
                offsets[1] -= offsetAdj[1]; #adjust the offset to try to get only whitespace there
            else:
                buffGoodArray[1] = True;
            #END IF
            #--- Bottom check ---
            if(  (np.sum(rgb_data[-bufferPix:,:,:]) != reqVals[2]) & (offsets[2] > offsetLatch[0,2]) & (offsets[2] < offsetLatch[1,2]) ):
                offsets[2] += offsetAdj[2]; #adjust the offset to try to get only whitespace there
            else:
                buffGoodArray[2] = True;
            #END IF
            #--- Top check ---
            if( (np.sum(rgb_data[0:bufferPix,:,:]) != reqVals[3]) & (offsets[3] < offsetLatch[0,3]) & (offsets[3] > offsetLatch[1,3]) ):
                offsets[3] -= offsetAdj[3]; #adjust the offset to try to get only whitespace there
            else:
                buffGoodArray[3] = True;
            #END IF
    
            fig.subplots_adjust(left = offsets[0], right = offsets[1], bottom = offsets[2], top = offsets[3]); #sets padding to small numbers for minimal white space
            fig.canvas.draw(); #key for all instances
            if( interactive == True ): #only needs fig.canvas.draw() for interactive, flush isn't req
                fig.canvas.flush_events(); #only needed for interactive plotting
            #END IF
            buffGood = np.all(buffGoodArray); #they all have to approve
            rgb_data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep=''); #get the rgb data in a numpy array, based on https://stackoverflow.com/questions/40433211/how-can-i-get-the-pixel-colors-in-matplotlib @ Jean-Sébastien
            rgb_data = rgb_data.reshape((int(len(rgb_data)/3), 3)).reshape(np.hstack((figSize,3))); #get the rgb data in a numpy array, based on https://stackoverflow.com/questions/40433211/how-can-i-get-the-pixel-colors-in-matplotlib @ Jean-Sébastien
        #END WHILE
    #END IF
    if( returnOffsets == True ):
        return offsets
    #END IF    
#END DEF