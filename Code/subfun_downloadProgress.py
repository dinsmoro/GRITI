#GOAL: Display progress percentage of downloading file - corresponds to urlretrieve in urllib.request library
#RD on 8/27/2018
#Greatly inspired by https://stackoverflow.com/questions/22676/how-do-i-download-a-file-over-http-using-python
#
#INPUT: urlretrieve takes care of it, call from urlretrieve
#OUTPUT: Says percentage downloaded

import sys

def downloadProgress(count, blockSize, totalSize): #
    if( totalSize != -1 ):
        percent = int(count*blockSize*100/totalSize);
        sys.stdout.write("\rDownload %%: %d%%" % percent);
    else:
        dl_mb = count*blockSize/1E6;
        sys.stdout.write("\rDownloaded MB: %.2f" % dl_mb);
    #END IF
    sys.stdout.flush();
 #END DEF