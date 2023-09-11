import os
import matplotlib.pyplot as plt
import hashlib

# Creates subfun_pltPause dynamically (keeps it up to date with Matplotlib so I never have to care again)
#basically an autopatcher, more lines than the pause function itself
def subfun_pltPause_create(cwd):
    code_path = os.path.join(cwd,'Code','subfun_pltPause.py');
    
    pyplot_path = plt.__file__; #path to the pyplot file which holds the pause function
    
    #--- Get hash of pyplot.py ---
    with open(pyplot_path, 'rb') as file:
        hashslingingslasher = hashlib.blake2b(); #blake2b is ez pz for quick looksie (equiv to sha3_256)
        while chunk := file.read(1024*768):
            hashslingingslasher.update(chunk); #fill 'er up
        #END WHILE
    #END WITH
    pyplot_hash = hashslingingslasher.hexdigest(); #hash it aight
    
    #--- Figure out if pltPause needs to be remade ---
    if( os.path.isfile(code_path) ):
        with open(code_path, 'r') as pause_file:
            if( pause_file.readline() == '#'+pyplot_hash+'\n' ):
                FLG_pltPause_make = False; #don't make it, pause function based off of same file
            else:
                FLG_pltPause_make = True; #make it, hash different so maybe pause function different
            #END IF
        #END WITH
    else:
        FLG_pltPause_make = True; #make it
    #END IF
    
    #--- Make pltPause if needed ---
    if( FLG_pltPause_make == True ):
        #--- Extract pause code from latest matplotlib ---
        FLG_saver = False; #prep saver flag
        saver = ['import time\n','from matplotlib import _pylab_helpers\n']; #prep saver list with time import and _pylab_helpers imports
        with open(pyplot_path, 'r') as pyplot_file:
            for line in pyplot_file:
                if( 'def pause(' in line ):
                    FLG_saver = True; #turn on the saver
                    saver.append(line); #start recording
                elif( FLG_saver == True ):
                    if( (FLG_saver == True) and ((line[0] != '@') and ('def ' not in line)) ):
                        saver.append(line); #keep recording
                    elif( (FLG_saver == True) and ((line[0] == '@') or ('def ' in line)) ):
                        FLG_saver = False; #quit recording
                        break;
                    #END IF
                #END IF
            #END FOR line
        #END WITH
        
        #--- Remove show call ---   
        try:
            saver.pop([idx for idx, s in enumerate(saver) if 'show(' in s][0]);
        except:
            pass; #nothing to do
        #END TRYING
        
        #--- Insert at 1st value the pyplot hash ---
        saver.insert(0, '#'+pyplot_hash+'\n');
        
        #--- Write to file ---
        with open(code_path, 'w') as pause_file:
            for line in saver:
                pause_file.write(line); #write it in
            #END FOR line
        #END WITH
    #END IF
#END DEF