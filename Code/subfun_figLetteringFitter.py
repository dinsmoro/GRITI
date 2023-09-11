"""
Adjusts a., b., c.,... subplot lettering to fit automagically
Coded while listening to https://www.youtube.com/watch?v=9lfzBrvhFoI on loop
"""
import numpy as np
import matplotlib.pyplot as plt

def subfun_figLetteringFitter(fig, ax, bufferPix = None, bufferPixEdge = None, tryLim = 4):
    #right now it works with the summed colors, switching to the 3 colorbands may be req'd if more similar colors are needed to be differentiated (but not my problem now, rip future me)
    letteringPix_letteringColor = 255; #sum of the lettering color
    letteringPix_backgroundColor = 765; #sum of the background color
    letteringPix_figTextColor = 0; #sum of the fig text color
    interactive = plt.isinteractive(); #record if interactive is true or false
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
    if( bufferPixEdge == None ): #if bufferPixEdge is none, estimate a good bufferPixEdge
        if( np.max(figSize) < 3000 ):
            bufferPixEdge = np.int64(np.round(np.max(figSize)*0.0016)); #pixels, automagically estimate the buffer pixel amount
        else:
            #if it's a lotta pixels, increase the buffer b/c it gets weird otherwise I'm not sure
            bufferPixEdge = np.int64(np.round(np.max(figSize)*0.0026)); #pixels, automagically estimate the buffer pixel amount
        #END IF
    #END IF
    bufferPixEdgeSq = bufferPixEdge**2; #pre-calc this
    jRec = 0; #prep
    for i in range(0,len(ax)): #check all axes
        for j in range(0,len(ax[i].texts)):
            if( ax[i].texts[j].get_text().find('.') > 0 ): #mayt need to improve this sometime
                jRec += 1; #record if it's real
        #END FOR j
    #END FOR i
    letteringPos = np.zeros((len(ax)*jRec,2)); #get the list to hold the positions ready
    letteringGood  = np.zeros((len(ax)*jRec), dtype=bool); #used to track all of the logicals
    tryOut = 1;
    while( np.any(letteringGood == False) & (tryOut <= tryLim) ):
        rgb_data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep=''); #get the rgb data in a numpy array, based on https://stackoverflow.com/questions/40433211/how-can-i-get-the-pixel-colors-in-matplotlib @ Jean-Sébastien
        rgb_data = rgb_data.reshape((int(len(rgb_data)/3), 3)).reshape(np.hstack((figSize,3))); #get the rgb data in a numpy array, based on https://stackoverflow.com/questions/40433211/how-can-i-get-the-pixel-colors-in-matplotlib @ Jean-Sébastien
        cntr = 0; #prep
        for i in range(0,len(ax)): #check all axes
            for j in range(0,len(ax[i].texts)):
                if( ax[i].texts[j].get_text().find('.') > 0 ): #mayt need to improve this sometime
                    letteringPositionX, letteringPositionY = np.array(ax[i].texts[j].get_position()); #get the position
                    bboxz = ax[i].texts[j].get_window_extent().get_points(); #get the size of the text based on Q https://stackoverflow.com/questions/44700065/matplotlib-direct-way-to-get-axes-coordinates-of-annotation-boxes
                    letteringSize = ax[i].transAxes.inverted().transform(bboxz); #get the box size
                    letteringPix = rgb_data[figSize[0]-np.int64(np.floor(bboxz[1,1])):figSize[0]-np.int64(np.ceil(bboxz[0,1])),np.int64(np.floor(bboxz[0,0])):np.int64(np.ceil(bboxz[1,0])),:]; #get the box where the text is
                    letteringPixSums = np.sum(letteringPix,axis=2);
                    #Remove anti-aliasing "muddled" pixels that aren't the full sums
                    letteringPixSums_antiAliased = (letteringPixSums != letteringPix_letteringColor) & (letteringPixSums != letteringPix_backgroundColor) & (letteringPixSums != letteringPix_figTextColor); #get any colors not the 3 primaries
                    #this all assumes lettering sum < background sum - may need extra logic for other situations
                    #letter & background mix
                    letteringPixSums_antiAliased_letterBackgroundMix_letter = (letteringPixSums > letteringPix_letteringColor) & (letteringPixSums <= np.int64((letteringPix_backgroundColor-letteringPix_letteringColor)/2)+letteringPix_letteringColor); #get where lettering and background mixing would be pushed to lettering
                    letteringPixSums_antiAliased_letterBackgroundMix_background = (letteringPixSums < letteringPix_backgroundColor) & (letteringPixSums > np.int64((letteringPix_backgroundColor-letteringPix_letteringColor)/2)+letteringPix_letteringColor); #get where lettering and background mixing would be pushed to background
                    letteringPixSums[letteringPixSums_antiAliased & letteringPixSums_antiAliased_letterBackgroundMix_letter] = letteringPix_letteringColor; #try to remove antialiasing
                    letteringPixSums[letteringPixSums_antiAliased & letteringPixSums_antiAliased_letterBackgroundMix_background] = letteringPix_backgroundColor; #try to remove antialiasing
                    #letter & fig text mix
                    letteringPixSums_antiAliased_letterFigTextMix_letter = (letteringPixSums < letteringPix_letteringColor) & (letteringPixSums > np.int64((letteringPix_letteringColor-letteringPix_figTextColor)/2)+letteringPix_figTextColor); #get where lettering and fig text mixing would be pushed to lettering
                    letteringPixSums_antiAliased_letterfigTextMix_figText = (letteringPixSums > letteringPix_figTextColor) & (letteringPixSums <= np.int64((letteringPix_letteringColor-letteringPix_figTextColor)/2)+letteringPix_figTextColor); #get where lettering and fig text mixing would be pushed to fig text
                    letteringPixSums[letteringPixSums_antiAliased & letteringPixSums_antiAliased_letterFigTextMix_letter] = letteringPix_letteringColor; #try to remove antialiasing
                    letteringPixSums[letteringPixSums_antiAliased & letteringPixSums_antiAliased_letterfigTextMix_figText] = letteringPix_figTextColor; #try to remove antialiasing
                    #recalc b/c background and fig text color straddles the letter color so need to recalc with that removed
                    letteringPixSums_antiAliased = (letteringPixSums != letteringPix_letteringColor) & (letteringPixSums != letteringPix_backgroundColor) & (letteringPixSums != letteringPix_figTextColor); #get any colors not the 3 primaries
                    #fig text & background mix (last assuming that letter color is between background and fig text colors)
                    letteringPixSums_antiAliased_figTextBackgroundMix_figText = (letteringPixSums > letteringPix_figTextColor) & (letteringPixSums <= np.int64((letteringPix_backgroundColor-letteringPix_figTextColor)/2)+letteringPix_figTextColor); #get where fig text and background mixing would be pushed to fig text
                    letteringPixSums_antiAliased_figTextBackgroundMix_background = (letteringPixSums < letteringPix_backgroundColor) & (letteringPixSums > np.int64((letteringPix_backgroundColor-letteringPix_figTextColor)/2)+letteringPix_figTextColor); #get where fig text and background mixing would be pushed to background                
                    letteringPixSums[letteringPixSums_antiAliased & letteringPixSums_antiAliased_figTextBackgroundMix_figText] = letteringPix_figTextColor; #try to remove antialiasing
                    letteringPixSums[letteringPixSums_antiAliased & letteringPixSums_antiAliased_figTextBackgroundMix_background] = letteringPix_backgroundColor; #try to remove antialiasing
                    #locate fig text
                    letteringPixSums_figText = letteringPixSums == letteringPix_figTextColor; #get where the fig text is at
                    if( letteringPixSums_figText.sum() != 0 ):
                        letteringPixSums_letter = np.copy(letteringPixSums); #copy it over
                        letteringPixSums_letter[letteringPixSums_figText] = letteringPix_backgroundColor; #fill the fig text in with the background color
                        #note that letteringPixSums is uint32 so np.diff on it will make big #s instead of negatives, not a problem for this implementation
                        letteringPixSums_letter_edges = np.vstack( (np.zeros((1,letteringPixSums.shape[1]),dtype=np.bool_),(np.diff(letteringPixSums_letter,axis=0) > 0)) ) | \
                            np.hstack( (np.zeros((letteringPixSums.shape[0],1),dtype=np.bool_),(np.diff(letteringPixSums_letter,axis=1) > 0)) ); #calc the locations of the edges of the letter
                        letteringPixSums_figText_edges = np.vstack( (np.zeros((1,letteringPixSums.shape[1]),dtype=np.bool_),(np.diff(letteringPixSums_figText,axis=0) > 0)) ) | \
                            np.hstack( (np.zeros((letteringPixSums.shape[0],1),dtype=np.bool_),(np.diff(letteringPixSums_figText,axis=1) > 0)) ); #calc the locations of the edges of the figure text
                        #calc distances between letter edges and fig text edges
                        letteringPixSums_letter_edges_where = np.where(letteringPixSums_letter_edges);
                        letteringPixSums_figText_edges_where = np.where(letteringPixSums_figText_edges);
                        letteringPixSums_letterFigText_edges_distance_X = np.repeat(letteringPixSums_letter_edges_where[0],letteringPixSums_figText_edges_where[0].size)-np.tile(letteringPixSums_figText_edges_where[0],letteringPixSums_letter_edges_where[0].size); #get the X component
                        letteringPixSums_letterFigText_edges_distance_Y = np.repeat(letteringPixSums_letter_edges_where[1],letteringPixSums_figText_edges_where[1].size)-np.tile(letteringPixSums_figText_edges_where[1],letteringPixSums_letter_edges_where[1].size); #get the Y component
                        letteringPixSums_letterFigText_edges_distance = (letteringPixSums_letterFigText_edges_distance_X)**2 + (letteringPixSums_letterFigText_edges_distance_Y)**2; #calc the distance w/o sqrt
                        if( letteringPixSums_letterFigText_edges_distance.min() < bufferPixEdgeSq ):
                            indexz = np.where(letteringPixSums_letterFigText_edges_distance.min() == letteringPixSums_letterFigText_edges_distance)[0]; #get where mins were
                            # anglez = np.arctan2(letteringPixSums_letterFigText_edges_distance_X[indexz],letteringPixSums_letterFigText_edges_distance_Y[indexz])*180/np.pi; #calc the angle
                            #adjust by the dist (even if it's not horizontal - it's just easier rn sorry future me if this really matters)
                            bboxzNew = np.copy(bboxz); #copy it over
                            bboxzNew[:,0] = bboxzNew[:,0] - bufferPixEdge-np.abs(letteringPixSums_letterFigText_edges_distance_Y[indexz]).max(); #get the delta to adjust horizontally by
                            if( bboxzNew[0,0] > bufferPixEdge ):
                                letteringSizeNew = ax[i].transAxes.inverted().transform(bboxzNew); #get the box size
                                letteringPos[cntr,:] = (letteringSizeNew[0,0]-letteringSize[0,0]+letteringPositionX,letteringPositionY); #record the positions
                            else:
                                #try again on Y axis
                                bboxzNew = np.copy(bboxz); #copy it over
                                maxz = np.abs(letteringPixSums_letterFigText_edges_distance_X[indexz]) == np.abs(letteringPixSums_letterFigText_edges_distance_X[indexz]).max(); #get the max value index, sprinkle maxes everywhere to deal with the issue where it could be multiple points or not - just want a scalar but don't want to do a lot of code fo rit
                                bboxzNew[:,1] = bboxzNew[:,1] + bufferPixEdge-letteringPixSums_letterFigText_edges_distance_X[indexz][maxz].max(); #get the delta to adjust vertically by
                                letteringSizeNew = ax[i].transAxes.inverted().transform(bboxzNew); #get the box size
                                letteringPos[cntr,:] = (letteringPositionX,letteringPositionY+letteringSizeNew[1,1]-letteringSize[1,1]); #record the positions
                            #END IF
                        else:
                            letteringPos[cntr,:] = (letteringPositionX,letteringPositionY); #record, no change needed if there's no fig text in the cutout
                            letteringGood[cntr] = True; #upgrade to true
                        #END IF
                    else:
                        letteringPos[cntr,:] = (letteringPositionX,letteringPositionY); #record, no change needed if there's no fig text in the cutout
                        letteringGood[cntr] = True; #upgrade to true
                    #END IF
                    cntr += 1; #increment
                #END IF
            #END FOR j
        #END FOR i
        if( np.any(letteringGood == False) == True ):
            letteringPositionXmax = letteringPos[:,0].min(); #get the smallest one (it's negative)
            cntr = 0; #prep cntr
            for i in range(0,len(ax)): #check all axes
                for j in range(0,len(ax[i].texts)):
                    if( ax[i].texts[j].get_text().find('.') > 0 ): #mayt need to improve this sometime
                        #this is designed so y can vary but the x position is the left most it could ever be so it looks consistent
                        ax[i].texts[j].set_position(np.array((letteringPositionXmax,letteringPos[cntr,1])) ); #set the position
                        cntr += 1; #increment cntr at end
                    #END IF
                #END FOR j
            #END FOR i
            #refresh figure to check again
            fig.canvas.draw(); #key for all instances
            if( interactive == True ): #only needs fig.canvas.draw() for interactive, flush and show and pause aren't req
                fig.canvas.flush_events(); #only needed for interactive plotting
            #END IF
        #END IF
        tryOut += 1; #increment try limit
    #END WHILE
    #--old implementation--
    # tryLim = 5; #number of tries to allow
    # bufferLim = 0.005; #axis coordinate buffer between the lettering and an axis label
    # for i in range(0,1): #only check a., b. doesn't need it and the alg'll have to be improved otherwise
    #     yTickCoords = np.vstack((np.zeros_like(ax[0].get_yticks()), ax[0].get_yticks())).T; #based on A https://stackoverflow.com/questions/27913674/matplotlib-get-axis-relative-tick-positions/27919217
    #     yTicks = np.flip(ax[0].transAxes.inverted().transform(ax[0].transData.transform(yTickCoords))[:,1]); #gets the yticks positions in relative axis values
        
    #     for j in range(0,len(ax[i].texts)):
    #        bboxz = ax[i].texts[j].get_window_extent(); #get the size of the text based on Q https://stackoverflow.com/questions/44700065/matplotlib-direct-way-to-get-axes-coordinates-of-annotation-boxes
    #        letteringSize = ax[i].transAxes.inverted().transform(bboxz); #get the box size
           
    #        tryCnt = 1; #prep
    #        while( np.any( (np.min(letteringSize[:,1]) <= yTicks) & (np.max(letteringSize[:,1]) >= yTicks)) & (tryCnt <= tryLim) ):
    #            letteringPos = np.array(ax[i].texts[j].get_position()); #get the position

    #            letteringPos[1] = yTicks[(np.min(letteringSize[:,1]) <= yTicks) & (np.max(letteringSize[:,1]) >= yTicks)][0] + bufferLim*tryCnt; #try moving it up
               
    #            ax[i].texts[j].set_position(letteringPos); #set the position
               
    #            tryCnt += 1; #increment
    #            bboxz = ax[i].texts[j].get_window_extent(); #get the size of the text based on Q https://stackoverflow.com/questions/44700065/matplotlib-direct-way-to-get-axes-coordinates-of-annotation-boxes
    #            letteringSize = ax[i].transAxes.inverted().transform(bboxz); #get the box size
    #        #END WHILE
    #    #END FOR j
    # #END FOR i
#END DEF