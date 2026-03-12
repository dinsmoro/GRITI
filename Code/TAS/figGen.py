#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from os import path as ospath, listdir as oslistdir, getcwd as osgetcwd # For if fileName_save is not provided but FLG_interactivePlots is False
from re import compile as recompile # For if fileName_save is not provided but FLG_interactivePlots is False
import matplotlib.pyplot as plt
from TAS.subfun_figFitter import figFitter

def figGen_starter(nrows=1, ncols=1, figsize=(14,8.5), dpi=150, FLG_interactivePlots=True):
    # --- Plot it up ---
    if( FLG_interactivePlots == True ):
        fig, ax = plt.subplots(nrows=nrows, ncols=ncols); #use instead of fig because it inits an axis too
        figManager = plt.get_current_fig_manager(); #req to maximize
        try:
            figManager.window.showMaximized(); #force maximized
        except:
            plt.close(); # Close figure to restart
            fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize); # Open figure with a set size
        # END TRYING
    else:
        if( FLG_interactivePlots == False ):
            plt.ioff() #disable showing the plot as its size will be larger than the screen, which cannot happen if the plot is shown
        # END IF
        fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize, dpi=dpi); #use instead of fig because it inits an axis too (I think I dunno)
    #END IF
    if( (nrows == 1) and (ncols == 1) ):
        ax = [ax]; # Standardize into a list, makes it easier to add more in the future
    else:
        ax = ax.ravel(); # Flatten so its easier to use lists to hold it
    # END IF
    
    return fig, ax
# END DEF

def figGen_ender(fig, offsets=None, fileName_save=None, FLG_tricky=True, FLG_interactivePlots=True, FLG_verbose=2):
    # --- Finalize ---
    if( offsets is None ):
        offsets = figFitter(fig, FLG_tricky=FLG_tricky, returnOffsets=True); #fit the fig fast
    else:
        fig.subplots_adjust(left = offsets[0], right = offsets[1], bottom = offsets[2], top = offsets[3]); # Use offsets from a previous run for SPEED
    # END IF
    
    if( FLG_interactivePlots != True ):
        if( fileName_save is None ):
            if( FLG_verbose >= 1 ):
                print('WARNING in figGen_ender: fileName_save is None but FLG_interactivePlots is NOT True. No file name to save as, defaulting to `fig#` saved in CWD as a ONG with # incrementing from last saved fig in folder.');
            # END IF
            filezDir2use = osgetcwd(); # Get it
            regexr = recompile(r'fig[0-9]+\.png'); # Compile it
            filezFound = oslistdir(filezDir2use); # Find some files
            filez = [filez for filez in filezFound if ospath.isfile(ospath.join(filezDir2use, filez)) and regexr.search(filez) is not None]; # Only get the right ones
            if( len(filez) > 0 ):
                filez.sort(key=lambda f: int(''.join(filter(str.isdigit, f)))); # Cursed filter stuff from https://stackoverflow.com/a/33159707
                fileName_save = ospath.join(filezDir2use, 'fig'+str(int(filez[-1][3:filez[-1].find('.')])+1)+'.png'); # Prep the name to save as (one more than previous max)
            else:
                fileName_save = ospath.join(filezDir2use, 'fig1.png'); # Prep the name to save as
            # END IF
        # END IF
        if( FLG_verbose >= 2 ):
            print('Plotting: '+fileName_save);
        # END IF
        #fig.savefig(os.path.join(smb.local, dir_sub_plots,os.path.splitext(filez[i])[0]+settings['plot']['save file type'])); #save the figure
        fig.savefig(fileName_save); #save the figure
        if( FLG_interactivePlots == False ):
            plt.close(); #close figure b/c it lurks apparently
            plt.ion(); #re-enable it for later stuff
        # END IF
    else:
        plt.draw(); # Seems this is needed
        plt.show(block=False); # Wait 
    #END IF
    
    return offsets
# END DEF