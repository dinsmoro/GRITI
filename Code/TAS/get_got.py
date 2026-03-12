#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from math import factorial
from os import path as ospath, remove as osremove
import sys
from time import time, sleep as timesleep
import urllib.request as urlreq
from urllib.error import URLError
import html2text
from warnings import warn

def get_got_webzine(webobj, web_retryMax=3, web_retryWait=[1,3,1], urlHandlers=None, FLG_rendered=True, FLG_showDownloadedAmnt=True, FLG_verbose=2):
    if( isinstance(webobj, str) ):
        FLG_webWhat = 'page'; # webobj is a string that's just the link to the webpage
    else:
        FLG_webWhat = 'file'; # webobj needs to be a list/tuple/numpy obj of size 2 (weblink2file2dl, webpath4file2dl2locally)
    # END IF
    
    web_retryNum = 1; #prep the retry num
    web_retryWait_val = web_retryWait[0]; #prep the retry wait period
    FLG_webWork = False; #lets while exit
    if( urlHandlers is not None ):
        # debugS_handler = urlreq.HTTPSHandler(debuglevel=10); # This does debug output on HTTPS
        url_opener = urlreq.build_opener(*urlHandlers); # Apply the urlHandlers (they must be in a list or tuple)
        urlreq.install_opener(url_opener);
    # END IF
    while( ((web_retryMax >= web_retryNum) | (web_retryMax == -1)) & (FLG_webWork == False) ):
        try:
            if( FLG_webWhat == 'page' ):
                page = urlreq.urlopen(webobj); # Get raw HTML
                FLG_webWork = True; # Flag for web working
            elif( FLG_webWhat == 'file' ):
                tic = time(); # For DL speed records
                if( (FLG_verbose >= 2) and (FLG_showDownloadedAmnt == True) ):
                    urlreq.urlretrieve(webobj[0], webobj[1], reporthook=downloadProgress); #download the file in question to the data directory
                else:
                    urlreq.urlretrieve(webobj[0], webobj[1]); #download the file in question to the data directory
                # END IF
                toc2tic = time() - tic; # For DL speed records
                FLG_webWork = True; # Exit the loop
            # END IF
        except URLError as err:
            FLG_reportError = True; # Set the reporter to off
            if( not callable(getattr(err, 'status', None)) ):
                reportError = 'likely no internet (`'+str(err)+'`)'; # Dunno what, sorry
                if( web_retryNum == 1 ):
                    if( FLG_webWhat == 'page' ):
                        print('WARNING in get_got.py: '+webobj+' connection failed.');
                    elif( FLG_webWhat == 'page' ):
                        print('WARNING in get_got.py: '+webobj[0]+' connection failed.');
                    # END IF
                # END IF
            elif( (err.status == 404) or (err.status == 429) ): # Treat 429 as a 404 effectively (seems to be used as 404 not as a real 429 for NOAA CORS, i.e., waiting will never yield it. fix handling for future 429, sorry future me)
                FLG_webWork = 2; # Flag for website not found error
                web_retryNum = web_retryMax; # If it is not there, there is no reason to look for it
                FLG_reportError = False; # Set to not report
            else:
                reportError = 'some reason (`'+str(err)+'`)'; # Dunno what, sorry
                
            # END IF
            if( FLG_reportError == True ):
                if( FLG_verbose >= 1 ):
                    sys.stdout.write('\rWARNING in get_got.py: Download failed for '+reportError+' ¯\_(ツ)_/¯. Try #{}/{}, will wait {} sec\t'.format(web_retryNum, web_retryMax, web_retryWait_val)); #report
                    sys.stdout.flush();
                # END IF
                FLG_webWork = False; # Flag for web no work
                timesleep(web_retryWait_val);
                web_retryNum += 1; #increment try
                if( web_retryWait_val < web_retryWait[1] ):
                    web_retryWait_val += web_retryWait[2]; # Increment by amount
                #END IF
            # END IF
        #END TRYING
    #END WHILE
    
    if( FLG_webWhat == 'page' ):
        return page, FLG_webWork, web_retryNum
    elif( FLG_webWhat == 'file' ):
        return FLG_webWork, toc2tic, web_retryNum
    # END IF
# END DEF

def get_got_webpg(weblink, web_retryMax=3, web_retryWait=[1,3,1], urlHandlers=None, FLG_rendered=True, FLG_verbose=2):
    page, FLG_webWork, web_retryNum = get_got_webzine(weblink, web_retryMax=web_retryMax, web_retryWait=web_retryWait, urlHandlers=urlHandlers, FLG_rendered=FLG_rendered, FLG_verbose=FLG_verbose); # Call a function that handles both pages and files - deals with the various errors better that way
    
    # web_retryNum = 1; # Prep the retry num
    # web_retryWait_val = web_retryWait[0]; # Prep the retry wait period
    # FLG_webWork = False; # Lets while exit
    # if( urlHandlers is not None ):
    #     # debugS_handler = urlreq.HTTPSHandler(debuglevel=10); # This does debug output on HTTPS
    #     url_opener = urlreq.build_opener(*urlHandlers); # Apply the urlHandlers (they must be in a list or tuple)
    #     urlreq.install_opener(url_opener);
    # # END IF
    # while( ((web_retryMax >= web_retryNum) | (web_retryMax == -1)) & (FLG_webWork == False) ):
    #     try:
    #         page = urlreq.urlopen(weblink); # Get raw HTML
    #         FLG_webWork = True; # Flag for web working
    #     except URLError as err:
    #         if( (err.status == 404) or (err.status == 429) ): # Treat 429 as a 404 effectively (seems to be used as 404 not as a real 429 for NOAA CORS, i.e., waiting will never yield it. fix handling for future 429, sorry future me)
    #             FLG_webWork = 2; # Flag for website not found error
    #             web_retryNum = web_retryMax; # If it is not there, there is no reason to look for it
    #         # elif( err.status == 429 ): # Something to deal with 429, but for now just treat it as 404
    #         #     if( FLG_verbose >= 2 ):
    #         #         if( web_retryMax < 5 ):
    #         #             web_retryMax = 5; # Increase the max # of waits for this
    #         #         # END IF
    #         #         sys.stdout.write('\rINFOin get_got.py: Web link `'+weblink+'` says "too many requests", chilling out #'+str(cntr_chill)+'.'); #report
    #         #         sys.stdout.flush();
    #         #         timesleep(60*factorial(web_retryNum))
    #         #     # END IF
    #         else:
    #             if( FLG_verbose >= 1 ):
    #                 sys.stdout.write('\rWARNING in get_got.py: No web access for data avail search. Try #{}/{}, will wait {} sec until next try\t'.format(web_retryNum,web_retryMax,web_retryWait_val)); #report
    #                 sys.stdout.flush();
    #             # END IF
    #             timesleep(web_retryWait_val);
    #             web_retryNum += 1; #increment try
    #             if( web_retryWait_val < web_retryWait[1] ):
    #                 web_retryWait_val += web_retryWait[2]; # Increment by amount
    #             #END IF
    #         # END IF
    #     #END TRYING
    # #END WHILE
    if( FLG_webWork == True ):
        if( (web_retryNum > 1) and (FLG_verbose >= 1) ):
            print('\n'); # For the no web access info
        #END IF
        html_content = page.read(); # Read off the HTML from whatever the page holder is
        charset = page.headers.get_content_charset(); # Get the charset from the page holder, w/o it doesn't work
        if( charset is None ):
            charset = 'utf-8'; # Assume utf-8
        #END IF
        html_content = html_content.decode(charset); # "decode" the HTML content so it's easily legible
        if( FLG_rendered == True ):
            return html2text.html2text(html_content) # Render the HTML like webpage would and get the real stuff
        else:
            return html_content # Return the direct HTML
        # END IF
    elif( FLG_webWork == 2 ):
        if( FLG_verbose >= 1 ):
            print('\nWARNING in get_got.py: Web page `'+weblink+'` returned "404: page not found", check the link validity.');
            return None
        # END IF
    else:
        raise Exception('\nERROR in get_got.py: No internet access, can\'t get data from site `'+weblink+'`.');
    #END IF
# END DEF

def downloadProgress(count, blockSize, totalSize):
    if( totalSize != -1 ):
        percent = int(count*blockSize*100/totalSize);
        sys.stdout.write("\rDownload %%: %d%%" % percent);
    else:
        dl_mb = count*blockSize/1E6;
        sys.stdout.write("\rDownloaded MB: %.2f" % dl_mb);
    #END IF
    sys.stdout.flush();
 #END DEF

def get_got_webfile(weblink2file, web_filePath, web_dlSpeed=None, web_retryMax=3, web_retryWait=[1,3,1], urlHandlers=None, FLG_getFileSize=True, FLG_showDownloadedAmnt=True, FLG_unzip=None, FLG_test4hdf5=True, FLG_ignoreFailure=False, FLG_errorOn404=False, FLG_verbose=2):
    if( FLG_getFileSize ):
        try:
            web_fileSize = np.asarray(urlreq.urlopen(weblink2file).info()['Content-Length'],dtype=np.int64);
        except TypeError:
            web_fileSize = None; # Website won't report file size
        # END TRYING
    else:
        web_fileSize = None; # Skip getting a file size, most relevant use is if you know the website won't give a file size
    # END IF
    if( web_fileSize is not None ):
        web_fileSize_report = web_fileSize/(1024**3);
        web_fileSize_reportText = 'GB';
        if( web_fileSize_report < 1 ): # If file size is actually in MB, report that
            web_fileSize_report = web_fileSize/(1024**2);
            web_fileSize_reportText = 'MB';
        # END IF
        if( FLG_verbose >= 2 ):
            if( web_dlSpeed is not None ):
                print("INFO in get_got.py: Downloading {} to \"{}\".\nFile size is {:.2f} {}. At {:.2f} MB/s ({:.2f} Mbps) expect it to take {:.2f} min.\nFrom site {}".format( ospath.basename(web_filePath), ospath.dirname(web_filePath), web_fileSize_report ,web_fileSize_reportText, web_dlSpeed, web_dlSpeed*8., (web_fileSize/1024**2*(1/web_dlSpeed))/60, weblink2file ));
            else:
                web_dlSpeed = 12.5; # MB/s assumption (100 Mbps)
                print("INFO in get_got.py: Downloading {} to \"{}\".\nFile size is {:.2f} {}. Download speed unknown; assuming {:.2f} MB/s ({:.2f} Mbps) expect it to take {:.2f} min.\nFrom site {}".format( ospath.basename(web_filePath), ospath.dirname(web_filePath), web_fileSize_report, web_fileSize_reportText, web_dlSpeed,web_dlSpeed*8.0,(web_fileSize/1024**2*(1/web_dlSpeed))/60, weblink2file ));
            # END IF
        # END IF
    # END IF
    FLG_gotPage, toc2tic, web_retryNum = get_got_webzine([weblink2file, web_filePath], web_retryMax=web_retryMax, web_retryWait=web_retryWait, urlHandlers=urlHandlers, FLG_showDownloadedAmnt=FLG_showDownloadedAmnt, FLG_verbose=FLG_verbose); # Call a function that handles both pages and files - deals with the various errors better that way
    
    # web_retryNum = 1; #prep the retry num
    # web_retryWait_val = web_retryWait[0]; #prep the retry wait period
    # FLG_gotPage = False; #lets while exit
    # if( urlHandlers is not None ):
    #     # debugS_handler = urlreq.HTTPSHandler(debuglevel=10); # This does debug output on HTTPS
    #     url_opener = urlreq.build_opener(*urlHandlers); # Apply the urlHandlers (they must be in a list or tuple)
    #     urlreq.install_opener(url_opener);
    # # END IF
    # while( ((web_retryMax >= web_retryNum) | (web_retryMax == -1)) & (FLG_gotPage == False) ):
    #     try:
    #         tic = time(); #for time testing
    #         if( FLG_verbose >= 2 ):
    #             urlreq.urlretrieve(weblink2file, web_filePath, reporthook=downloadProgress); #download the file in question to the data directory
    #         else:
    #             urlreq.urlretrieve(weblink2file, web_filePath); #download the file in question to the data directory
    #         # END IF
    #         toc2tic = time() - tic; #for time testing
    #         FLG_gotPage = True; #exit the loop
    #     except URLError as err:
    #         if( (err.status == 404) or (err.status == 429) ): # Treat 429 as a 404 effectively (seems to be used as 404 not as a real 429 for NOAA CORS, i.e., waiting will never yield it. fix handling for future 429, sorry future me)
    #             FLG_gotPage = 2; # Flag for website not found error
    #             web_retryNum = web_retryMax; # If it is not there, there is no reason to look for it
    #         else:
    #             if( FLG_verbose >= 1 ):
    #                 sys.stdout.write('\rWARNING in get_got.py: Download failed for some reason ¯\_(ツ)_/¯. Try #{}/{}, will wait {} sec\t'.format(web_retryNum, web_retryMax, web_retryWait_val)); #report
    #                 sys.stdout.flush();
    #             # END IF
    #             FLG_gotPage = False; # Flag for web no work
    #             timesleep(web_retryWait_val);
    #             web_retryNum += 1; #increment try
    #             if( web_retryWait_val < web_retryWait[1] ):
    #                 web_retryWait_val += web_retryWait[2]; # Increment by amount
    #             #END IF
    #         # END IF
    #     #END TRYING
    # #END WHILE
    if( FLG_unzip is not None ):
        if( FLG_unzip != False ): # False could also disable I guess
            if( FLG_unzip == True ):
                FLG_unzip = ospath.splitext(web_filePath)[1]; # Try to automagically determine
            # END IF
            FLG_unzip = FLG_unzip.lower().replace('.',''); # Lower it, remove periods
            if( FLG_unzip in ['gz', 'gzip'] ):
                # Soft req imports
                import gzip
                import shutil
                with gzip.open(web_filePath, 'rb') as finput:
                    with open(ospath.splitext(web_filePath)[0], 'wb') as foutput:
                        shutil.copyfileobj(finput, foutput); # Woozle
                    # END WITH
                # END WITH
                osremove(web_filePath); # Delete the zipped
                web_filePath = ospath.splitext(web_filePath)[0]; # Remove that extension
            else:
                raise Exception('ERROR in get_got.py: Unzip format requested `'+str(FLG_unzip)+'` that is NOT supported.');
            # END IF
        # END IF
    # END IF
    if( (FLG_test4hdf5 == True) and (FLG_gotPage == True) ):
        import h5py, netCDF4 # Soft req
        
        # Test the file to make sure it is useable
        def test_file(file2test, typez):
            if( typez in ['h5', 'hdf5'] ):
                if( 'bytes' in str(type(file2test)) ):
                    import io # Soft import req
                    file2test = io.BytesIO(file2test); # Convert to bytesIO for h5py
                # END IF
                with h5py.File(file2test, 'r', rdcc_nbytes=500*1024*1024) as testFile:
                    testFile.keys(); # Tries to check out some stuff in the file
                # END WITH
            elif( typez in ['nc', 'nc4', 'netcdf4'] ):
                if( "'str'" in str(type(file2test)) ):
                    with netCDF4.Dataset(file2test, 'r') as testFile:
                        testFile.variables.keys(); # Tries to check out some stuff in the file
                    # END WITH
                elif( 'bytes' in str(type(file2test)) ):
                    with netCDF4.Dataset('inmem.nc', 'r', memory=file2test) as testFile:
                        testFile.variables.keys(); # Tries to check out some stuff in the file
                    # END WITH
                else:
                    raise Exception('ERROR in get_got.py: File downloaded link\n'+weblink2file+'\nto path\n'+web_filePath+'\nis NOT known file format (`'+str(type(file2test))+'` provided).');
                # END IF
            else:
                raise Exception('ERROR in get_got.py: File downloaded link\n'+weblink2file+'\nto path\n'+web_filePath+'\nis NOT known file type (`'+typez+'` provided).');
            # END IF
        # END DEF
        
        # Check if zipped
        FLG_unzip = ospath.splitext(web_filePath)[1].lower().replace('.',''); # Lower it, remove periods
        try:
            if( FLG_unzip in ['gz', 'gzip'] ):
                # Soft req imports
                import gzip
                with gzip.open(web_filePath, 'rb') as finput:
                    finputObj = finput.read(); # Read in-place
                    test_file(finputObj, ospath.splitext(ospath.splitext(web_filePath)[0])[1].lower().replace('.','')); # Test that file!
                # END WITH
            else:
                test_file(web_filePath, ospath.splitext(web_filePath)[1].lower().replace('.','')); # Test that file!
            # END IF
        except: # (OSError , KeyError):
            warn('WARNING in get_got.py: File downloaded link\n'+weblink2file+'\nto path\n'+web_filePath+'\nis NOT readable. DELETING it. Check that file manually.'); # Millstone Hill was serving 72 byte gibberish files, warn and go by b/c there was no good data
        # END TRYING
        
        if( web_fileSize is not None ):
            web_dlSpeed = web_fileSize/(1024**2)/toc2tic; #MB/s, calc current web download speed
            if( FLG_verbose >= 2 ):
                print("\nINFO in get_got.py: Time to download: {:.2f} min, download speed: {:.2f} MB/s\n\n".format(toc2tic/60,web_dlSpeed)); # Extra space at end
            # END IF
        # END IF
        return web_dlSpeed
    elif( FLG_gotPage == 2 ):
        if( FLG_errorOn404 == True ):
            raise Exception('ERROR in get_got.py: File `'+weblink2file+'` returned "404: page not found", check the link validity.');
        else:
            if( FLG_verbose >= 1 ):
                print('WARNING in get_got.py: File `'+weblink2file+'` returned "404: page not found", check the link validity.');
            # END IF
        # END IF
    elif( FLG_gotPage == False ):
        if( FLG_ignoreFailure == False ):
            raise Exception('ERROR in get_got.py: Failed to download the file from link\n'+weblink2file+'\nto path\n'+web_filePath);
        else:
            if( FLG_verbose >= 1 ):
                print('WARNING in get_got.py: Failed to download the file from link\n'+weblink2file+'\nto path\n'+web_filePath);
            # END IF
        # END IF
    # END IF
# END DEF