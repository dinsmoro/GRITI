# Renames files AND TEXT INSIDE FILES from the original string to the new string
# It's got no safety checks and will destroy everything if used wronk
# --- var definitions ---
# folder is the folder to operate on
# strOrig is the original string to find
# strRep is the replacement string to replace
# fileEnding restricts the file list to just the provided file ending including the period (e.g. '.py'), default is off (an untennable file ending '@nope?')
# recursive allows the renamer go to into subfolders if true, default is false
# renameTextInFiles enables renaming the text in the files, default true
# renameFileNames enables renaming the file names, default true
# verbose reports stats on renaming, default true

import os
from subfun_strstr import strstr

# from renamer import renamer
def renamer(folder,strOrig,strRep,fileEnding='.@nope?',recursive=False,renameTextInFiles=True,renameFileNames=True,verbose=True):
    #--- Get the file names and their full paths ---
    files = []; #prep a list
    for (dirpath, dirnames, filenames) in os.walk(folder): # adapted from https://stackoverflow.com/a/3207973/2403531
        ditcher = []; #prep a ditch list (lists are silly sometimes)
        for i in range(0,len(filenames)):
            if( (fileEnding == '.@nope?') | (filenames[i].endswith(fileEnding) == True) ):
                filenames[i] = os.path.join(dirpath, filenames[i]); #tack that on
            else:
                ditcher.append(i); #record to ditch it
            #END IF
        #END FOR i
        filenames = [v for i, v in enumerate(filenames) if i not in ditcher]; #ditch em, adapted from https://stackoverflow.com/a/21032057/2403531 b/c lists am I right
        
        files += filenames; #append the full list
        if( recursive == False ): #prevents going into subdirectories
            break;
        #END IF
    #END FOR (dirpath, dirnames, filenames)
    
    #--- Rename text inside each file ---
    if( renameTextInFiles == True ):
        cntr = 0; #prep cntr
        cntrOccur = 0; #prep cntr
        for i in range(0,len(files)):
            try:
                cntrSub = 0; #prep cntr
                textFileEdited = []; #prep a list to hold the text file
                with open(files[i],'r') as textFile:
                    for line in textFile:
                        cntrSub += strstr(line,strOrig).size; #record how many changed
                        textFileEdited.append(line.replace(strOrig,strRep)); #append the line and replace any instances of strOrig with strRep in each line
                    #END FOR line
                #END WITH
            except(UnicodeDecodeError):
                cntrSub = 0; #prep cntr
                textFileEdited = []; #prep a list to hold the text file
                with open(files[i],'r', encoding='utf8') as textFile:
                    for line in textFile:
                        cntrSub += strstr(line,strOrig).size; #record how many changed
                        textFileEdited.append(line.replace(strOrig,strRep)); #append the line and replace any instances of strOrig with strRep in each line
                    #END FOR line
                #END WITH
            #END TRYING
            if( cntrSub > 0 ): #only rewrite if needed
                with open(files[i],'w') as textFile: #replace the file with a new, edited file
                    for j in range(0,len(textFileEdited)):
                        textFile.write('%s' % textFileEdited[j]); #write in each edited line
                    #END FOR j
                #END WITH
                if( verbose == True ):
                    print(files[i]+' replaced "'+strOrig+'" '+str(cntrSub)+' times.'); #report
                #END IF
                cntr += cntrSub; #increment main cntr
                cntrOccur += 1; #increment occurance cntr
            #END IF
        #END FOR i
        if( verbose == True ):
            print('In '+str(cntrOccur)+' of '+str(len(files))+' files replaced "'+strOrig+'" with "'+strRep+'" a total of '+str(cntr)+' times.'); #report
        #END IF
    #END IF
    
    #--- Rename file names ---
    if( renameFileNames == True ):
        cntr = 0; #prep cntr
        for i in range(0,len(files)):
            if( files[i].replace(strOrig,strRep) != files[i] ):
                os.rename(files[i], files[i].replace(strOrig,strRep)); #rename the files to end
                cntr += 1; #increment
            #END IF
        #END FOR i
        if( verbose == True ):
            print('Renamed '+str(cntr)+' files out of '+str(len(files))+' from "'+strOrig+'" to "'+strRep+'".'); #report
        #END IF
    #END IF
#END DEF