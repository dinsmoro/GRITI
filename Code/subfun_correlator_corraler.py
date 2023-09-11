#  data1 stays static and data2 slides if FLG_correlator == 1
#  data2 stays static and data1 slides if FLG_correlator == 2
# tabulator saves length(data1_setNames) number of CSVs
# -NEG lags mean the sliding data was lagged (backward) in time, +POS lags mean the sliding data was pushed forward in time

import numpy as np
import joblib
from copy import deepcopy
from Code.subfun_filter import subfun_filter
from Code.subfun_timeMatch import subfun_timeMatch
from Code.subfun_correlator import subfun_correlator
from Code.subfun_strfind import strfind
# from Code.subfun_strstr import strstr
from Code.subfun_textNice import textNice

def subfun_correlator_corraler(
        data1, settings1, name1, data1_setNames, \
        data2, settings2, name2, data2_setNames, \
        dates, settings_plot, settings_paths, settings_config, \
        FLG_correlator, FLG_correlator_options, FLG_correlator_plot = False, \
        FLG_correlator_tabulator = True, FLG_tabulator_extraBit = '', FLG_enableText = True, \
        FLG_correlator_shiftDir = 'both', FLG_correlator_timeLim = 14400, FLG_correlator_timeBound = None, \
        filt1 = None,  filt2 = None, settings_spectra=None, reportDivisor=[60,'min']):

    dateRange_dayNum_full = dates['date range full dayNum'];
            
    if( not isinstance(data1_setNames,list) and not isinstance(data1_setNames,tuple) ):
        data1_setNames = [data1_setNames]; #make it a list
    #END IF
    if( not isinstance(data2_setNames,list) and not isinstance(data2_setNames,tuple) ):
        data2_setNames = [data2_setNames]; #make it a list
    #END IF
    data1_setNames = data1_setNames.copy(); #protect it from possible edits
    data2_setNames = data2_setNames.copy(); #protect it from possible edits
    
    if( settings1 == '!DIRECT' ):
        #---FOR THIS CASE data2 MUST BE A LIST OF [data2, timeUnique2] so things still work without the whole dict construct---
        #---name1 and data1_setNames should represent the <name1=name of overarching data type> and <data1_setNames=sub-set of data>--- 
        #(for real data it would be name1='OMNI', data1_setNames='SYM/H') 
        #(like for double keo data, name1='Delta-vTEC E USA' and data1_setNames='Delta-vTEC E USA' too so the name and tabulation files are named rightish)
        #denotes direct data/label were sent in, not a "standardized" dict
        data1 = {data1_setNames[0]:data1[0], 'time unique':data1[1], 'data rate':np.median(np.diff(data1[1]))}; #wrap it so it works
        settings1 = {'labels':{data1_setNames[0]:data1_setNames[0]},'data type':data1_setNames[0]}; #wrap it so it works
    #END IF
    
    if( settings2 == '!DIRECT' ):
        #---FOR THIS CASE data2 MUST BE A LIST OF [data2, timeUnique2] so things still work without the whole dict construct---
        #---name1 and data1_setNames should represent the <name1=name of overarching data type> and <data1_setNames=sub-set of data>--- 
        #(for real data it would be name1='OMNI', data1_setNames='SYM/H') 
        #(like for double keo data, name1='Delta-vTEC E USA' and data1_setNames='Delta-vTEC E USA' too so the name and tabulation files are named rightish)
        #denotes direct data/label were sent in, not a "standardized" dict
        data2 = {data2_setNames[0]:data2[0], 'time unique':data2[1], 'data rate':np.median(np.diff(data2[1]))}; #wrap it so it works
        settings2 = {'labels':{data2_setNames[0]:data2_setNames[0]},'data type':data2_setNames[0]}; #wrap it so it works
    #END IF
    
    if( FLG_tabulator_extraBit != '' ):
        FLG_tabulator_extraBit = '_'+FLG_tabulator_extraBit; #tack on a _ at the front
    #END IF
    
    #work through all names to detect sub-dicts
    if( FLG_enableText == True ):
        if( FLG_correlator == 1 ):
            print('\nWorking on Correlator '+name1+' n '+name2+' ('+name2+' slides w/ +POS lags going back in time b/c '+name2+' must be pushed fwd to align w/ '+name1+')');
        else:
            print('\nWorking on Correlator '+name2+' n '+name1+' ('+name1+' slides w/ +POS lags going back in time b/c '+name1+' must be pushed fwd to align w/ '+name2+')');
        #END IF
    #END IF
    
    FLG_fancyPlot = settings_plot['fancy plot'];
            
    def parallel_helper(data1, data1_alias, kr_data1, filt1, data1_setNames_d1, data1_adder, settings1, data2, kr_data2, filt2, data2_setNames_d2, settings2, \
                        kr_dataB_extended, FLG_correlator_options_d2, FLG_correlator_shiftDir, FLG_correlator_plot, settings_plot, settings_spectra, dateRange_zeroHr_day):
        
        data2_alias = data2; #for cheeky situations
        # FLG_data2_aliasFix = False; #default
        data2_adder = ''; #nothing usually
        
        #&| means sub-dict, $| means sub-array
        if( (strfind(data2_setNames_d2,'&|',opt=1) > 0) | (strfind(data2_setNames_d2,'$|',opt=1) > 0) ):
            data2_adder = '$\mathregular{_T}$';
            data2_setNames_subSegs = data2_setNames_d2.split('&|'); #get the sub-dicts
            data2_setNames_d2 = data2_setNames_subSegs[0]; #set the 1st thing to be the thing for the label to work
            for d2s in range(0,len(data2_setNames_subSegs)):
                if( strfind(data2_setNames_subSegs[d2s],'$|',opt=1) == 0 ): #check for sub-array for the sub-dict
                    data2_alias = data2_alias[data2_setNames_subSegs[d2s]]; #move down in the data structure
                else:
                    data2_setNames_subSegs_subVects = data2_setNames_subSegs[d2s].split('$|'); #split the sub-dict from the sub-vect
                    data2_alias = data2_alias[data2_setNames_subSegs_subVects[0]]; #move down in the data structure
                    if( data2_alias.shape[0] > data2_alias.shape[1] ):
                        data2_alias = data2_alias[:,int(data2_setNames_subSegs_subVects[1])]; #choose this array based on smaller array dim being the component
                    else:
                        data2_alias = data2_alias[int(data2_setNames_subSegs_subVects[1]),:]; #choose this array based on smaller array dim being the component
                    #END IF
                    if( data2_setNames_subSegs_subVects[1] == '0' ):
                        data2_adder = '$\mathregular{_X}$';
                    elif( data2_setNames_subSegs_subVects[1] == '1' ):
                        data2_adder = '$\mathregular{_Y}$';
                    elif( data2_setNames_subSegs_subVects[1] == '2' ):
                        data2_adder = '$\mathregular{_Z}$';
                    elif( data2_setNames_subSegs_subVects[1] == '3' ):
                        data2_adder = '$\mathregular{_T}$'; #(By^2 + Bz^2)^(1/2) - not implemented yet
                    #END IF
                #END IF
            #END FOR d2s
            data2_alias = {data2_setNames_d2:data2_alias}; #align it just right so it works with the code later
            # FLG_data2_aliasFix = True; #gotta fix the data2_alias later
        #END IF
        
        if( FLG_correlator == 1 ):
            sig1 = data1_alias[data1_setNames_d1][kr_data1]; #sig1 stays static
            time1 = data1['time unique'][kr_data1]-dateRange_zeroHr_day*86400;
            dataRate1 = data1['data rate'];
            if( (FLG_correlator_options_d2['mode'] != 'corr') ):
                #no need to limit if not corr because other modes cut out time frames and can use extra data
                sig2 = data2_alias[data2_setNames_d2][kr_dataB_extended]; #sig2 moves around
                time2 = data2['time unique'][kr_dataB_extended]-dates['date range zero hr dayNum'][1]*86400;
            else:
                sig2 = data2_alias[data2_setNames_d2][kr_data2]; #sig2 moves around
                time2 = data2['time unique'][kr_data2]-dates['date range zero hr dayNum'][1]*86400;
            #END IF
            dataRate2 = data2['data rate'];
            plotName = settings1['labels'][settings1['data type']]+data1_adder+ \
                ' & '+settings2['labels'][data2_setNames_d2]+data2_adder;
        else:
            sig1 = data2_alias[data2_setNames_d2][kr_data2]; #sig1 stays static
            time1 = data2['time unique'][kr_data2]-dateRange_zeroHr_day*86400;
            dataRate1 = data2['data rate'];
            if( (FLG_correlator_options_d2['mode'] != 'corr') ):
                #no need to limit if not corr because other modes cut out time frames and can use extra data
                sig2 = data1_alias[data1_setNames_d1][kr_dataB_extended]; #sig2 moves around
                time2 = data1['time unique'][kr_dataB_extended]-dates['date range zero hr dayNum'][1]*86400;
            else:
                sig2 = data1_alias[data1_setNames_d1][kr_data1]; #sig2 moves around
                time2 = data1['time unique'][kr_data1]-dates['date range zero hr dayNum'][1]*86400;
            #END IF
            dataRate2 = data1['data rate'];
            plotName = settings2['labels'][data2_setNames_d2]+data2_adder+ \
                ' & '+settings1['labels'][settings1['data type']]+data1_adder;
        #END IF
        # if( FLG_data2_aliasFix == True ): #in parallel don't need to undo the data2 aliasing
        #     data2_alias = data2; #for cheeky situations
        #     data2_adder = ''; #nothing usually
        # #END IF
        
        # #--- Filter if Needed Before Time-Match ---
        if( filt1 != None ):
            sig1 = subfun_filter( np.copy(sig1), filt1, dataTime = time1, dataRate = dataRate1, settings_spectra = settings_spectra, reduceWindow = 1, FLG_reportNaNs = False); #filter (or not)
        #END IF
        if( filt2 != None ):
            sig2 = subfun_filter( np.copy(sig2), filt2, dataTime = time2, dataRate = dataRate2, settings_spectra = settings_spectra, reduceWindow = 1, FLG_reportNaNs = False); #filter (or not)
        #END IF
        
        if( np.isclose(dataRate1,dataRate2) == False ):
            if( dataRate1 > dataRate2 ):
                if( np.all(kr_data2) == False ):
                    #attempt to keep all sig2 data available by extending time1 so that sig2 can be time matched properly but extend past the real sig1/time1
                    time1p = time1.copy(); #copy
                    if( time1[-1] < time2[-1] ):
                        time1p = np.concatenate( (time1p, np.arange(time1[-1]+dataRate1, time2[-1]+dataRate1, dataRate1)) ); #tack it on
                    #END IF
                    if( time1[0] > time2[0] ):
                        time1p = np.concatenate( (np.flip(np.arange(time1[0]-dataRate1, time2[0]-dataRate1, -dataRate1)), time1p) ); #tack it on
                    #END IF
                    sig2, time2 = subfun_timeMatch(sig2, time2, time1p, timeMatch_delta=dataRate1, FLG_removeNaNs=2, FLG_reportNaNs=False); #match the times to the same cadence
                else:
                    sig2, time2 = subfun_timeMatch(sig2, time2, time1, timeMatch_delta=dataRate1, FLG_removeNaNs=2, FLG_reportNaNs=False); #match the times to the same cadence
                #END IF
                dataRate2 = dataRate1; #set it
            else:
                if( (FLG_correlator_options_d2['mode'] != 'corr') & (np.all(kr_data2) == False) ):
                    #need to limit time2 here so sig1 can time match right if it wasn't limited earlier
                    sig1, time1 = subfun_timeMatch(sig1, time1, time2[kr_data2], timeMatch_delta=dataRate2, FLG_removeNaNs=2, FLG_reportNaNs=False); #match the times to the same cadence
                else:
                    sig1, time1 = subfun_timeMatch(sig1, time1, time2, timeMatch_delta=dataRate2, FLG_removeNaNs=2, FLG_reportNaNs=False); #match the times to the same cadence
                #END IF
                dataRate1 = dataRate2; #set it
            #END IF
        #END IF
        
        if( FLG_correlator_options_d2['mode'] == 'corr' ):
            corrRet_once = subfun_correlator(sig1, sig2, mode='corr', time1=time1, time2=time2, sig1_filt=None, sig2_filt=None, settings_spectra=settings_spectra, reportDivisor=reportDivisor, FLG_interpGaps=True); #calc some correlation
        elif( FLG_correlator_options_d2['mode'] == 'shift' ):
            corrRet_once = subfun_correlator(sig1, sig2, mode='shift', time1=time1, time2=time2, dataRate=dataRate1, sig1_filt=None, sig2_filt=None, \
                reportDivisor=reportDivisor, FLG_shiftDir=FLG_correlator_shiftDir, FLG_interpGaps=True, FLG_plot=FLG_correlator_plot, settings_plot=settings_plot, \
                settings_spectra=settings_spectra, plotName=plotName, FLG_enableText=False, FLG_fancyPlot=0); #calc some correlation
        elif( FLG_correlator_options_d2['mode'] == 'range' ):
            corrRet_once = subfun_correlator(sig1, sig2, mode='range', timeLimit=FLG_correlator_timeLim, time1=time1, time2=time2, dataRate=dataRate1, sig1_filt=None, sig2_filt=None, timeRange=FLG_correlator_options_d2['time range'], \
                reportDivisor=reportDivisor, FLG_shiftDir=FLG_correlator_shiftDir, FLG_interpGaps=True, FLG_plot=FLG_correlator_plot, settings_plot=settings_plot, \
                settings_spectra=settings_spectra, plotName=plotName, FLG_enableText=False, FLG_fancyPlot=0); #calc some correlation
        elif( FLG_correlator_options_d2['mode'] == 'interval' ):
            corrRet_once = subfun_correlator(sig1, sig2, mode='interval', timeLimit=FLG_correlator_timeLim, time1=time1, time2=time2, dataRate=dataRate1, sig1_filt=None, sig2_filt=None, timeInterval=FLG_correlator_options_d2['time interval'], intervalType=FLG_correlator_options_d2['interval type'], \
                reportDivisor=reportDivisor, FLG_shiftDir=FLG_correlator_shiftDir, FLG_interpGaps=True, FLG_plot=FLG_correlator_plot, settings_plot=settings_plot, \
                settings_spectra=settings_spectra, plotName=plotName, FLG_enableText=False, FLG_fancyPlot=0); #calc some correlation
        elif( FLG_correlator_options_d2['mode'] == 'interval manual' ):
            # this is a special call that builds time timeInterval off of 
            if( FLG_correlator_options_d2['method'] == 'req' ):
                if( FLG_correlator_options_d2['data type'] == 'OMNI' ):
                    timeTemp = data2['time unique']-dates['date range zero hr dayNum'][1]*86400;
                    if( FLG_correlator_options_d2['comparator'] == '>' ):
                        where_indexes = np.where( (data2[FLG_correlator_options_d2['sub data type']] > FLG_correlator_options_d2['comparator val']) == True )[0];
                    elif( FLG_correlator_options_d2['comparator'] == '>=' ):
                        where_indexes = np.where( (data2[FLG_correlator_options_d2['sub data type']] >= FLG_correlator_options_d2['comparator val']) == True )[0];
                    elif( FLG_correlator_options_d2['comparator'] == '<' ):
                        where_indexes = np.where( (data2[FLG_correlator_options_d2['sub data type']] < FLG_correlator_options_d2['comparator val']) == True )[0];
                    elif( FLG_correlator_options_d2['comparator'] == '<=' ):
                        where_indexes = np.where( (data2[FLG_correlator_options_d2['sub data type']] <= FLG_correlator_options_d2['comparator val']) == True )[0];
                    elif( FLG_correlator_options_d2['comparator'] == '==' ):
                        where_indexes = np.where( (data2[FLG_correlator_options_d2['sub data type']] == FLG_correlator_options_d2['comparator val']) == True )[0];
                    elif( FLG_correlator_options_d2['comparator'] == '!=' ):
                        where_indexes = np.where( (data2[FLG_correlator_options_d2['sub data type']] != FLG_correlator_options_d2['comparator val']) == True )[0];
                    #END IF
                    where_edges = np.where(np.insert(np.diff(where_indexes),0,0) != 1)[0];
                    timeMinEnforcer = FLG_correlator_options_d2['time enforcer']; #copy it out, gonna overwrite FLG_correlator_options
                    FLG_correlator_options_d2['time interval'] = np.zeros( (where_edges.size-1,2) ); #prep
                    for j in range(0,where_edges.size-1):
                        FLG_correlator_options_d2['time interval'][j,0] = timeTemp[where_indexes[where_edges[j]]];
                        FLG_correlator_options_d2['time interval'][j,1] = timeTemp[where_indexes[where_edges[j+1]-1]];
                    #END FOR j
                    kj = np.where(np.diff(FLG_correlator_options_d2['time interval'],axis=1) > timeMinEnforcer)[0]; #enforce the minimum time distance
                    FLG_correlator_options_d2['time interval'] = FLG_correlator_options_d2['time interval'][kj,:]; #get values that meet the time reqs
                elif( FLG_correlator_options_d2['data type'] == 'AMPERE' ):
                    timeTemp = data1['time unique']-dateRange_zeroHr_day*86400;
                    if( FLG_correlator_options_d2['comparator'] == '>' ):
                        where_indexes = np.where( (data1[FLG_correlator_options_d2['sub data type']] > FLG_correlator_options_d2['comparator val']) == True )[0];
                    elif( FLG_correlator_options_d2['comparator'] == '>=' ):
                        where_indexes = np.where( (data1[FLG_correlator_options_d2['sub data type']] >= FLG_correlator_options_d2['comparator val']) == True )[0];
                    elif( FLG_correlator_options_d2['comparator'] == '<' ):
                        where_indexes = np.where( (data1[FLG_correlator_options_d2['sub data type']] < FLG_correlator_options_d2['comparator val']) == True )[0];
                    elif( FLG_correlator_options_d2['comparator'] == '<=' ):
                        where_indexes = np.where( (data1[FLG_correlator_options_d2['sub data type']] <= FLG_correlator_options_d2['comparator val']) == True )[0];
                    elif( FLG_correlator_options_d2['comparator'] == '==' ):
                        where_indexes = np.where( (data1[FLG_correlator_options_d2['sub data type']] == FLG_correlator_options_d2['comparator val']) == True )[0];
                    elif( FLG_correlator_options_d2['comparator'] == '!=' ):
                        where_indexes = np.where( (data1[FLG_correlator_options_d2['sub data type']] != FLG_correlator_options_d2['comparator val']) == True )[0];
                    #END IF
                    where_edges = np.where(np.insert(np.diff(where_indexes),0,0) != 1)[0];
                    timeMinEnforcer = FLG_correlator_options_d2['time enforcer']; #copy it out, gonna overwrite FLG_correlator_options
                    FLG_correlator_options_d2['time interval'] = np.zeros( (where_edges.size-1,2) ); #prep
                    for j in range(0,where_edges.size-1):
                        FLG_correlator_options_d2['time interval'][j,0] = timeTemp[where_indexes[where_edges[j]]];
                        FLG_correlator_options_d2['time interval'][j,1] = timeTemp[where_indexes[where_edges[j+1]-1]];
                    #END FOR j
                    kj = np.where(np.diff(FLG_correlator_options_d2['time interval'],axis=1) > timeMinEnforcer)[0]; #enforce the minimum time distance
                    FLG_correlator_options_d2['time interval'] = FLG_correlator_options_d2['time interval'][kj,:]; #get values that meet the time reqs
                #END IF
            #END IF
            corrRet_once = subfun_correlator(sig1, sig2, mode='interval manual', timeLimit=FLG_correlator_timeLim, time1=time1, time2=time2, dataRate=dataRate1, sig1_filt=None, sig2_filt=None, timeInterval=FLG_correlator_options_d2['time interval'], \
                reportDivisor=reportDivisor, FLG_shiftDir=FLG_correlator_shiftDir, FLG_interpGaps=True, FLG_plot=FLG_correlator_plot, settings_plot=settings_plot, \
                settings_spectra=settings_spectra, plotName=plotName, FLG_enableText=False, FLG_fancyPlot=0); #calc some correlation
        else:
            print('WARNING in subfun_correlator_corraler: FLG_correlator - Unknown option "'+FLG_correlator_options_d2['mode']+'", skipping analysis for #'+str(d2));
        #END IF
        
        return corrRet_once
    #END DEF
    
    if( len(FLG_correlator_options) == len(data2_setNames) ):
        # data1[data1_setNames] = GRITI_AMPERE_integrator(data1, dates, settings1, plotLatRange, plotLongRange, AMPERE_integrateMethod, AMPERE_integrateMethod_val, AMPERE_integrateMethod_log=AMPERE_integrateMethod_log); #integrate the AMPERE data according to some settings
        corrRet = [[None for d2 in range(0,len(data2_setNames))] for d1 in range(0,len(data1_setNames))]; #preallocate
        if( np.any(FLG_correlator_timeBound == None) ):
            kr_data2 = (data2['time unique'] >= dateRange_dayNum_full[0,1]*86400) & (data2['time unique'] < (dateRange_dayNum_full[-1,1]+1)*86400); #get stuff outside the date range
            kr_data1 = (data1['time unique'] >= dateRange_dayNum_full[0,1]*86400) & (data1['time unique'] < (dateRange_dayNum_full[-1,1]+1)*86400); #get stuff outside the date range
            if( FLG_correlator == 1 ): #to cover extended time snipping that affects whatever ends up as sig2 (which is why they're different)
                kr_dataB_extended = (data2['time unique'] >= (dateRange_dayNum_full[0,1]-np.ceil(FLG_correlator_timeLim/86400).astype(np.int64))*86400) & (data2['time unique'] <= (dateRange_dayNum_full[-1,1]+1+np.ceil(FLG_correlator_timeLim/86400).astype(np.int64))*86400); #get stuff outside the date range
            else:
                kr_dataB_extended = (data1['time unique'] >= (dateRange_dayNum_full[0,1]-np.ceil(FLG_correlator_timeLim/86400).astype(np.int64))*86400) & (data1['time unique'] <= (dateRange_dayNum_full[-1,1]+1+np.ceil(FLG_correlator_timeLim/86400).astype(np.int64))*86400); #get stuff outside the date range
            #END IF
        else:
            kr_data2 = (data2['time unique'] >= FLG_correlator_timeBound[0]) & (data2['time unique'] < FLG_correlator_timeBound[-1]); #get stuff outside the date range
            kr_data1 = (data1['time unique'] >= FLG_correlator_timeBound[0]) & (data1['time unique'] < FLG_correlator_timeBound[-1]); #get stuff outside the date range
            if( FLG_correlator == 1 ): #to cover extended time snipping that affects whatever ends up as sig2 (which is why they're different)
                kr_dataB_extended = (data2['time unique'] >= (FLG_correlator_timeBound[0]-np.ceil(FLG_correlator_timeLim/86400).astype(np.int64)*86400)) & (data2['time unique'] <= (FLG_correlator_timeBound[-1]+(1+np.ceil(FLG_correlator_timeLim/86400).astype(np.int64))*86400)); #get stuff outside the date range
            else:
                kr_dataB_extended = (data1['time unique'] >= (FLG_correlator_timeBound[0]-np.ceil(FLG_correlator_timeLim/86400).astype(np.int64)*86400)) & (data1['time unique'] <= (FLG_correlator_timeBound[-1]+(1+np.ceil(FLG_correlator_timeLim/86400).astype(np.int64))*86400)); #get stuff outside the date range
            #END IF
        #END IF
        
        data1_alias = data1; #for cheeky situations  
        FLG_data1_aliasFix = False;
        data1_adder = ''; #nothing usually
        if( (FLG_correlator_plot == True) | (len(data2_setNames) == 1) ):
            data2_alias = data2; #for cheeky situations
            FLG_data2_aliasFix = False;
            data2_adder = ''; #nothing usually
        #END IF
        FLG_error = False; #notes an error
        for d1 in range(0, len(data1_setNames)):
            
            #&| means sub-dict, $| means sub-array
            if( (strfind(data1_setNames[d1],'&|',opt=1) > 0) | (strfind(data1_setNames[d1],'$|',opt=1) > 0) ):
                data1_adder = '$\mathregular{_F}$'; #default to F
                data1_setNames_subSegs = data1_setNames[d1].split('&|'); #get the sub-dicts
                data1_setNames[d1] = data1_setNames_subSegs[0]; #set the 1st thing to be the thing for the label to work
                for d1s in range(0,len(data1_setNames_subSegs)):
                    if( strfind(data1_setNames_subSegs[d1s],'$|',opt=1) == 0 ): #check for sub-array for the sub-dict
                        data1_alias = data1_alias[data1_setNames_subSegs[d1s]]; #move down in the data structure
                    else:
                        data1_setNames_subSegs_subVects = data1_setNames_subSegs[d1s].split('$|'); #split the sub-dict from the sub-vect
                        data1_alias = data1_alias[data1_setNames_subSegs_subVects[0]]; #move down in the data structure
                        if( data1_alias.shape[0] > data1_alias.shape[1] ):
                            data1_alias = data1_alias[:,int(data1_setNames_subSegs_subVects[1])]; #choose this array based on smaller array dim being the component
                        else:
                            data1_alias = data1_alias[int(data1_setNames_subSegs_subVects[1]),:]; #choose this array based on smaller array dim being the component
                        #END IF
                        if( data1_setNames_subSegs_subVects[1] == '0' ):
                            data1_adder = '$\mathregular{_X}$';
                        elif( data1_setNames_subSegs_subVects[1] == '1' ):
                            data1_adder = '$\mathregular{_Y}$';
                        elif( data1_setNames_subSegs_subVects[1] == '2' ):
                            data1_adder = '$\mathregular{_Z}$';
                        elif( data1_setNames_subSegs_subVects[1] == '3' ):
                            data1_adder = '$\mathregular{_T}$'; #(By^2 + Bz^2)^(1/2) - not implemented yet
                        #END IF
                    #END IF
                #END FOR d2s
                data1_alias = {data1_setNames[d1]:data1_alias}; #align it just right so it works with the code later
                FLG_data1_aliasFix = True; #gotta fix the data1_alias later
            #END IF
            
            if( (FLG_correlator_plot == False) & (len(data2_setNames) > 1) ): #no point to parallelize 1
                parallel_list = [None for d2 in range(0,len(data2_setNames))]; #Prep
                for d2 in range(0,len(data2_setNames)):
                    parallel_list[d2] = [data1, data1_alias, kr_data1, filt1, data1_setNames[d1], data1_adder, settings1, data2, kr_data2, filt2, data2_setNames[d2], settings2, \
                                        kr_dataB_extended, FLG_correlator_options[d2], FLG_correlator_shiftDir, FLG_correlator_plot, settings_plot, settings_spectra, dates['date range zero hr dayNum'][1]]; #load in the variables needed
                #END FOR d2
                
                with joblib.parallel_backend('loky'):
                    parallel_return = joblib.Parallel(n_jobs=settings_config['parallel num threads'])(joblib.delayed(parallel_helper)(h, i, j, k, l, m, n, o ,p, q, r, s, t, u, v, w, x, y, z) for h, i, j, k, l, m, n, o ,p, q, r, s, t, u, v, w, x, y, z in parallel_list); #will this not destroy the world?
                #END WITH
                
                for d2 in range(0,len(data2_setNames)):
                    corrRet[d1][d2] = deepcopy(parallel_return[d2]); #unpack
                    #--- loky won't let stuff print, so replicate the text ---
                    if( isinstance(corrRet[d1][d2],str) ):
                        #catch for an error that was returned
                        print(corrRet[d1][d2]);
                        FLG_error = True;
                    else:
                        if( FLG_enableText == True ):
                            if( corrRet[d1][d2]['plot name'] == '' ):
                                print('-- Correlation Calcing ---');
                            else:
                                print('-- On: '+corrRet[d1][d2]['plot name'].replace('$\mathregular{','').replace('}$','')+' ---');
                            #END IF
                            print(corrRet[d1][d2]['text results']);
                        #END IF
                    #END IF
                #END FOR d2
            else:
                data2_alias = data2; #for cheeky situations
                FLG_data2_aliasFix = False; #default
                data2_adder = ''; #nothing usually
                for d2 in range(0,len(data2_setNames)):
                    #&| means sub-dict, $| means sub-array
                    if( (strfind(data2_setNames[d2],'&|',opt=1) > 0) | (strfind(data2_setNames[d2],'$|',opt=1) > 0) ):
                        data2_adder = '$\mathregular{_F}$'; #default to F
                        data2_setNames_subSegs = data2_setNames[d2].split('&|'); #get the sub-dicts
                        data2_setNames[d2] = data2_setNames_subSegs[0]; #set the 1st thing to be the thing for the label to work
                        for d2s in range(0,len(data2_setNames_subSegs)):
                            if( strfind(data2_setNames_subSegs[d2s],'$|',opt=1) == 0 ): #check for sub-array for the sub-dict
                                data2_alias = data2_alias[data2_setNames_subSegs[d2s]]; #move down in the data structure
                            else:
                                data2_setNames_subSegs_subVects = data2_setNames_subSegs[d2s].split('$|'); #split the sub-dict from the sub-vect
                                data2_alias = data2_alias[data2_setNames_subSegs_subVects[0]]; #move down in the data structure
                                if( data2_alias.shape[0] > data2_alias.shape[1] ):
                                    data2_alias = data2_alias[:,int(data2_setNames_subSegs_subVects[1])]; #choose this array based on smaller array dim being the component
                                else:
                                    data2_alias = data2_alias[int(data2_setNames_subSegs_subVects[1]),:]; #choose this array based on smaller array dim being the component
                                #END IF
                                if( data2_setNames_subSegs_subVects[1] == '0' ):
                                    data2_adder = '$\mathregular{_X}$';
                                elif( data2_setNames_subSegs_subVects[1] == '1' ):
                                    data2_adder = '$\mathregular{_Y}$';
                                elif( data2_setNames_subSegs_subVects[1] == '2' ):
                                    data2_adder = '$\mathregular{_Z}$';
                                elif( data2_setNames_subSegs_subVects[1] == '3' ):
                                    data2_adder = '$\mathregular{_T}$'; #(By^2 + Bz^2)^(1/2) - not implemented yet
                                #END IF
                            #END IF
                        #END FOR d2s
                        data2_alias = {data2_setNames[d2]:data2_alias}; #align it just right so it works with the code later
                        FLG_data2_aliasFix = True; #gotta fix the data2_alias later
                    #END IF
                    
                    if( FLG_correlator == 1 ):
                        sig1 = data1_alias[data1_setNames[d1]][kr_data1]; #sig1 stays static
                        time1 = data1['time unique'][kr_data1]-dates['date range zero hr dayNum'][1]*86400;
                        dataRate1 = data1['data rate'];
                        if( (FLG_correlator_options[d2]['mode'] != 'corr') ):
                            #no need to limit if not corr because other modes cut out time frames and can use extra data
                            sig2 = data2_alias[data2_setNames[d2]][kr_dataB_extended]; #sig2 moves around
                            time2 = data2['time unique'][kr_dataB_extended]-dates['date range zero hr dayNum'][1]*86400;
                        else:
                            sig2 = data2_alias[data2_setNames[d2]][kr_data2]; #sig2 moves around
                            time2 = data2['time unique'][kr_data2]-dates['date range zero hr dayNum'][1]*86400;
                        #END IF
                        dataRate2 = data2['data rate'];
                        plotName = settings1['labels'][settings1['data type']]+data1_adder+ \
                            ' & '+settings2['labels'][data2_setNames[d2]]+data2_adder;
                    else:
                        sig1 = data2_alias[data2_setNames[d2]][kr_data2]; #sig1 stays static
                        time1 = data2['time unique'][kr_data2]-dates['date range zero hr dayNum'][1]*86400;
                        dataRate1 = data2['data rate'];
                        if( (FLG_correlator_options[d2]['mode'] != 'corr') ):
                            #no need to limit if not corr because other modes cut out time frames and can use extra data
                            sig2 = data1_alias[data1_setNames[d1]][kr_dataB_extended]; #sig2 moves around
                            time2 = data1['time unique'][kr_dataB_extended]-dates['date range zero hr dayNum'][1]*86400;
                        else:
                            sig2 = data1_alias[data1_setNames[d1]][kr_data1]; #sig2 moves around
                            time2 = data1['time unique'][kr_data1]-dates['date range zero hr dayNum'][1]*86400;
                        #END IF
                        dataRate2 = data1['data rate'];
                        plotName = settings2['labels'][data2_setNames[d2]]+data2_adder+ \
                            ' & '+settings1['labels'][settings1['data type']]+data1_adder;
                    #END IF
                    if( FLG_data2_aliasFix == True ):
                        data2_alias = data2; #for cheeky situations
                        data2_adder = ''; #nothing usually
                        FLG_data2_aliasFix = False; #turn off
                    #END IF
                    
                    #--- Filter if Needed Before Time-Match ---
                    if( filt1 != None ):
                        sig1 = subfun_filter( np.copy(sig1), filt1, dataTime = time1, dataRate = dataRate1, settings_spectra = settings_spectra, reduceWindow = 1, FLG_reportNaNs = False); #filter (or not)
                    #END IF
                    if( filt2 != None ):
                        sig2 = subfun_filter( np.copy(sig2), filt2, dataTime = time2, dataRate = dataRate2, settings_spectra = settings_spectra, reduceWindow = 1, FLG_reportNaNs = False); #filter (or not)
                    #END IF
                    
                    if( np.isclose(dataRate1,dataRate2) == False ):
                        if( dataRate1 > dataRate2 ):
                            if( np.all(kr_data2) == False ):
                                #attempt to keep all sig2 data available by extending time1 so that sig2 can be time matched properly but extend past the real sig1/time1
                                time1p = time1.copy(); #copy
                                if( time1[-1] < time2[-1] ):
                                    time1p = np.concatenate( (time1p, np.arange(time1[-1]+dataRate1, time2[-1]+dataRate1, dataRate1)) ); #tack it on
                                #END IF
                                if( time1[0] > time2[0] ):
                                    time1p = np.concatenate( (np.flip(np.arange(time1[0]-dataRate1, time2[0]-dataRate1, -dataRate1)), time1p) ); #tack it on
                                #END IF
                                sig2, time2 = subfun_timeMatch(sig2, time2, time1p, timeMatch_delta=dataRate1, FLG_removeNaNs=2, FLG_reportNaNs=False); #match the times to the same cadence
                            else:
                                sig2, time2 = subfun_timeMatch(sig2, time2, time1, timeMatch_delta=dataRate1, FLG_removeNaNs=2, FLG_reportNaNs=False); #match the times to the same cadence
                            #END IF
                            dataRate2 = dataRate1; #set it
                        else:
                            if( (FLG_correlator_options[d2]['mode'] != 'corr') & (np.all(kr_data2) == False) ):
                                #need to limit time2 here so sig1 can time match right if it wasn't limited earlier
                                sig1, time1 = subfun_timeMatch(sig1, time1, time2[kr_data2], timeMatch_delta=dataRate2, FLG_removeNaNs=2, FLG_reportNaNs=False); #match the times to the same cadence
                            else:
                                sig1, time1 = subfun_timeMatch(sig1, time1, time2, timeMatch_delta=dataRate2, FLG_removeNaNs=2, FLG_reportNaNs=False); #match the times to the same cadence
                            #END IF
                            dataRate1 = dataRate2; #set it
                        #END IF
                    #END IF
                    
                    if( FLG_correlator_options[d2]['mode'] == 'corr' ):
                        corrRet[d1][d2] = subfun_correlator(sig1, sig2, mode='corr', time1=time1, time2=time2, sig1_filt=None, sig2_filt=None, settings_spectra=settings_spectra, reportDivisor=reportDivisor, FLG_interpGaps=True); #calc some correlation
                    elif( FLG_correlator_options[d2]['mode'] == 'shift' ):
                        corrRet[d1][d2] = subfun_correlator(sig1, sig2, mode='shift', time1=time1, time2=time2, dataRate=dataRate1, sig1_filt=None, sig2_filt=None, reportDivisor=reportDivisor, \
                            FLG_shiftDir=FLG_correlator_shiftDir, FLG_interpGaps=True, FLG_plot=FLG_correlator_plot, settings_plot=settings_plot, settings_paths=settings_paths, \
                            settings_spectra=settings_spectra, plotName=plotName, FLG_enableText=FLG_enableText, FLG_fancyPlot=FLG_fancyPlot); #calc some correlation
                    elif( FLG_correlator_options[d2]['mode'] == 'range' ):
                        corrRet[d1][d2] = subfun_correlator(sig1, sig2, mode='range', timeLimit=FLG_correlator_timeLim, time1=time1, time2=time2, dataRate=dataRate1, sig1_filt=None, sig2_filt=None, timeRange=FLG_correlator_options[d2]['time range'], \
                            reportDivisor=reportDivisor, FLG_shiftDir=FLG_correlator_shiftDir, FLG_interpGaps=True, FLG_plot=FLG_correlator_plot, settings_plot=settings_plot, settings_paths=settings_paths, \
                            settings_spectra=settings_spectra, plotName=plotName, FLG_enableText=FLG_enableText, FLG_fancyPlot=FLG_fancyPlot); #calc some correlation
                    elif( FLG_correlator_options[d2]['mode'] == 'interval' ):
                        corrRet[d1][d2] = subfun_correlator(sig1, sig2, mode='interval', timeLimit=FLG_correlator_timeLim, time1=time1, time2=time2, dataRate=dataRate1, sig1_filt=None, sig2_filt=None, timeInterval=FLG_correlator_options[d2]['time interval'], intervalType=FLG_correlator_options[d2]['interval type'], \
                            reportDivisor=reportDivisor, FLG_shiftDir=FLG_correlator_shiftDir, FLG_interpGaps=True, FLG_plot=FLG_correlator_plot, settings_plot=settings_plot, settings_paths=settings_paths, \
                            settings_spectra=settings_spectra, plotName=plotName, FLG_enableText=FLG_enableText, FLG_fancyPlot=FLG_fancyPlot); #calc some correlation
                    elif( FLG_correlator_options[d2]['mode'] == 'interval manual' ):
                        # this is a special call that builds time timeInterval off of 
                        if( FLG_correlator_options[d2]['method'] == 'req' ):
                            if( FLG_correlator_options[d2]['data type'] == 'OMNI' ):
                                timeTemp = data2['time unique']-dates['date range zero hr dayNum'][1]*86400;
                                if( FLG_correlator_options[d2]['comparator'] == '>' ):
                                    where_indexes = np.where( (data2[FLG_correlator_options[d2]['sub data type']] > FLG_correlator_options[d2]['comparator val']) == True )[0];
                                elif( FLG_correlator_options[d2]['comparator'] == '>=' ):
                                    where_indexes = np.where( (data2[FLG_correlator_options[d2]['sub data type']] >= FLG_correlator_options[d2]['comparator val']) == True )[0];
                                elif( FLG_correlator_options[d2]['comparator'] == '<' ):
                                    where_indexes = np.where( (data2[FLG_correlator_options[d2]['sub data type']] < FLG_correlator_options[d2]['comparator val']) == True )[0];
                                elif( FLG_correlator_options[d2]['comparator'] == '<=' ):
                                    where_indexes = np.where( (data2[FLG_correlator_options[d2]['sub data type']] <= FLG_correlator_options[d2]['comparator val']) == True )[0];
                                elif( FLG_correlator_options[d2]['comparator'] == '==' ):
                                    where_indexes = np.where( (data2[FLG_correlator_options[d2]['sub data type']] == FLG_correlator_options[d2]['comparator val']) == True )[0];
                                elif( FLG_correlator_options[d2]['comparator'] == '!=' ):
                                    where_indexes = np.where( (data2[FLG_correlator_options[d2]['sub data type']] != FLG_correlator_options[d2]['comparator val']) == True )[0];
                                #END IF
                                where_edges = np.where(np.insert(np.diff(where_indexes),0,0) != 1)[0];
                                timeMinEnforcer = FLG_correlator_options[d2]['time enforcer']; #copy it out, gonna overwrite FLG_correlator_options
                                FLG_correlator_options[d2]['time interval'] = np.zeros( (where_edges.size-1,2) ); #prep
                                for j in range(0,where_edges.size-1):
                                    FLG_correlator_options[d2]['time interval'][j,0] = timeTemp[where_indexes[where_edges[j]]];
                                    FLG_correlator_options[d2]['time interval'][j,1] = timeTemp[where_indexes[where_edges[j+1]-1]];
                                #END FOR j
                                kj = np.where(np.diff(FLG_correlator_options[d2]['time interval'],axis=1) > timeMinEnforcer)[0]; #enforce the minimum time distance
                                FLG_correlator_options[d2]['time interval'] = FLG_correlator_options[d2]['time interval'][kj,:]; #get values that meet the time reqs
                            elif( FLG_correlator_options[d2]['data type'] == 'AMPERE' ):
                                timeTemp = data1['time unique']-dates['date range zero hr dayNum'][1]*86400;
                                if( FLG_correlator_options[d2]['comparator'] == '>' ):
                                    where_indexes = np.where( (data1[FLG_correlator_options[d2]['sub data type']] > FLG_correlator_options[d2]['comparator val']) == True )[0];
                                elif( FLG_correlator_options[d2]['comparator'] == '>=' ):
                                    where_indexes = np.where( (data1[FLG_correlator_options[d2]['sub data type']] >= FLG_correlator_options[d2]['comparator val']) == True )[0];
                                elif( FLG_correlator_options[d2]['comparator'] == '<' ):
                                    where_indexes = np.where( (data1[FLG_correlator_options[d2]['sub data type']] < FLG_correlator_options[d2]['comparator val']) == True )[0];
                                elif( FLG_correlator_options[d2]['comparator'] == '<=' ):
                                    where_indexes = np.where( (data1[FLG_correlator_options[d2]['sub data type']] <= FLG_correlator_options[d2]['comparator val']) == True )[0];
                                elif( FLG_correlator_options[d2]['comparator'] == '==' ):
                                    where_indexes = np.where( (data1[FLG_correlator_options[d2]['sub data type']] == FLG_correlator_options[d2]['comparator val']) == True )[0];
                                elif( FLG_correlator_options[d2]['comparator'] == '!=' ):
                                    where_indexes = np.where( (data1[FLG_correlator_options[d2]['sub data type']] != FLG_correlator_options[d2]['comparator val']) == True )[0];
                                #END IF
                                where_edges = np.where(np.insert(np.diff(where_indexes),0,0) != 1)[0];
                                timeMinEnforcer = FLG_correlator_options[d2]['time enforcer']; #copy it out, gonna overwrite FLG_correlator_options
                                FLG_correlator_options[d2]['time interval'] = np.zeros( (where_edges.size-1,2) ); #prep
                                for j in range(0,where_edges.size-1):
                                    FLG_correlator_options[d2]['time interval'][j,0] = timeTemp[where_indexes[where_edges[j]]];
                                    FLG_correlator_options[d2]['time interval'][j,1] = timeTemp[where_indexes[where_edges[j+1]-1]];
                                #END FOR j
                                kj = np.where(np.diff(FLG_correlator_options[d2]['time interval'],axis=1) > timeMinEnforcer)[0]; #enforce the minimum time distance
                                FLG_correlator_options[d2]['time interval'] = FLG_correlator_options[d2]['time interval'][kj,:]; #get values that meet the time reqs
                            #END IF
                        #END IF
                        corrRet[d1][d2] = subfun_correlator(sig1, sig2, mode='interval manual', timeLimit=FLG_correlator_timeLim, time1=time1, time2=time2, dataRate=dataRate1, sig1_filt=None, sig2_filt=None, timeInterval=FLG_correlator_options[d2]['time interval'], \
                            reportDivisor=reportDivisor, FLG_shiftDir=FLG_correlator_shiftDir, FLG_interpGaps=True, FLG_plot=FLG_correlator_plot, settings_plot=settings_plot, settings_paths=settings_paths, \
                            settings_spectra=settings_spectra, plotName=plotName, FLG_enableText=FLG_enableText, FLG_fancyPlot=FLG_fancyPlot); #calc some correlation
                    else:
                        print('WARNING in subfun_correlator_corraler: FLG_correlator - Unknown option "'+FLG_correlator_options[d2]['mode']+'", skipping analysis for #'+str(d2));
                    #END IF
                #END FOR d2
            #END IF
            if( FLG_data1_aliasFix == True ):
                data1_alias = data1; #for cheeky situations 
                data1_adder = ''; #nothing usually
                FLG_data1_aliasFix = False; #turn off
            #END IF
        #END FOR d1
    else:
        print('WARNING in subfun_correlator_corraler: FLG_correlator - Length of FLG_correlator_options ('+str(len(FLG_correlator_options))+') does not equal length of data2_setNames ('+str(len(data2_setNames))+')');
    #END IF
    if( (FLG_correlator_tabulator == True) & ('corrRet' in locals()) & (not FLG_error) ):
        if( (FLG_correlator_options.count(FLG_correlator_options[0]) == len(FLG_correlator_options)) & (FLG_correlator_options[0]['mode'] == 'range') ):
            import pandas as pd
            import os
            if( not os.path.isdir(os.path.join(settings_paths['cwd'],'Tabulations')) ):
                os.makedirs(os.path.join(settings_paths['cwd'],'Tabulations'));
            #END IF
            for d1 in range(0,len(data1_setNames)): #run through all data1_setNames and save individually
                #gather data into vector
                corr_fileName = os.path.join(settings_paths['cwd'],'Tabulations','correlator_tabulator_'+name1+'-'+settings1['labels'][data1_setNames[d1]].replace(' ','_')+FLG_tabulator_extraBit+'.csv');
                corr_colName = np.array(FLG_correlator_options[0]['time range'])/3600; #prep, conv to hr
                corr_colNameStr = textNice(corr_colName).replace('[','').replace(']','')+' hr';
                corr_rowNames = data2_setNames.copy(); #prep
                #--rename rowNames as needed--
                for d2 in range(0,len(corr_rowNames)): #run through all data1_setNames and save individually
                    data2_setNames_d2 = corr_rowNames[d2]; #get it out
                    #&| means sub-dict, $| means sub-array
                    if( (strfind(data2_setNames_d2,'&|',opt=1) > 0) | (strfind(data2_setNames_d2,'$|',opt=1) > 0) ):
                        data2_setNames_subSegs = data2_setNames_d2.split('&|'); #get the sub-dicts
                        data2_setNames_d2 = data2_setNames_subSegs[0]; #set the 1st thing to be the thing for the label to work
                        for d2s in range(0,len(data2_setNames_subSegs)):
                            if( strfind(data2_setNames_subSegs[d2s],'$|',opt=1) != 0 ): #check for sub-array for the sub-dict
                                data2_setNames_subSegs_subVects = data2_setNames_subSegs[d2s].split('$|'); #split the sub-dict from the sub-vect
                                if( data2_setNames_subSegs_subVects[1] == '0' ):
                                    data2_setNames_d2 = data2_setNames_d2+'_X';
                                elif( data2_setNames_subSegs_subVects[1] == '1' ):
                                    data2_setNames_d2 = data2_setNames_d2+'_Y';
                                elif( data2_setNames_subSegs_subVects[1] == '2' ):
                                    data2_setNames_d2 = data2_setNames_d2+'_Z';
                                elif( data2_setNames_subSegs_subVects[1] == '3' ):
                                    data2_setNames_d2 = data2_setNames_d2+'_T';
                                #END IF
                            #END IF
                        #END FOR d2s
                        corr_rowNames[d2] = data2_setNames_d2; #put it back (prob not needed since python memory sharing stuff but explicit is good)
                    #END IF
                #END FOR d2
                
                corr_holdr = ['' for d2 in range(0,len(corrRet[d1]))]; #prep
                for d2 in range(0,len(corrRet[d1])):
                    if( np.any(corrRet[d1][d2]['corr'] != None) ):
                        corr_holdr[d2] = textNice(np.round(corrRet[d1][d2]['corr'].item(),3)).replace('0.','.')+' @ '+textNice(corrRet[d1][d2]['time shift'].item()/60);#print resupts
                    #END IF
                #END FOR d2  
                FLG_changed = False; #keeps from writing if not needed
                if( os.path.isfile(corr_fileName) ):
                    #read it in for editz
                    corr_csv = pd.read_csv(corr_fileName,index_col=0);
                    corr_csvColumns = list(corr_csv.columns);
                    corr_csvRows = list(corr_csv.index);
                    if( corr_colNameStr in corr_csvColumns ):
                        for d2 in range(0,len(corr_rowNames)):
                            if( corr_rowNames[d2] in corr_csvRows ):
                                if( (corr_csv.loc[corr_rowNames[d2],corr_colNameStr] != corr_holdr[d2]) & (corr_holdr[d2] != '') ):
                                    corr_csv.loc[corr_rowNames[d2],corr_colNameStr] = corr_holdr[d2]; #load in latest info
                                    FLG_changed = True;
                                elif( (corr_csv.loc[corr_rowNames[d2],corr_colNameStr] != corr_holdr[d2]) & (corr_holdr[d2] == '') ):
                                    print('WARNING in subfun_correlator_corraler: FLG_correlator_tabulator column name '+corr_colNameStr+' & row name '+corr_rowNames[d2]+' not recorded due to corr coeff being None.');
                                #END IF
                            else:
                                if( corr_holdr[d2] != '' ):
                                    corr_csv.loc[corr_rowNames[d2],corr_colNameStr] = corr_holdr[d2]; #load in latest info
                                    FLG_changed = True;
                                else:
                                    print('WARNING in subfun_correlator_corraler: FLG_correlator_tabulator column name '+corr_colNameStr+' & row name '+corr_rowNames[d2]+' not recorded due to corr coeff being None.');
                                #END IF
                            #END IF
                        #END FOR d2
                    else:
                        if( not all(j == '' for j in corr_holdr) ):
                            corr_csv = pd.concat((corr_csv, pd.DataFrame(corr_holdr, corr_rowNames, [corr_colNameStr])),axis=1); #write it in in one go
                            FLG_changed = True;
                        else:
                            print('WARNING in subfun_correlator_corraler: FLG_correlator_tabulator column name '+corr_colNameStr+' & row names '+str(corr_rowNames)+' not recorded due to corr coeff being None.');
                        #END IF
                    #END IF
                else:
                    if( not all(j == '' for j in corr_holdr) ):
                        #create
                        corr_csv = pd.DataFrame(corr_holdr, corr_rowNames, [corr_colNameStr]); #write it in in one go
                        FLG_changed = True;
                    else:
                        print('WARNING in subfun_correlator_corraler: FLG_correlator_tabulator column name '+corr_colNameStr+' & row names '+str(corr_rowNames)+' not recorded due to corr coeff being None.');
                    #END IF
                #END IF
                #write to csv
                if( FLG_changed == True ):
                    corr_csv.to_csv(corr_fileName);
                #END IF
            #END FOR d1
        else:
            print('WARNING in subfun_correlator_corraler: FLG_correlator_tabulator is ON but FLG_correlator_options are NOT identical or mode is NOT range. Printing FLG_correlator_options:\n'+str(FLG_correlator_options)); #report issue
        #END IF
    #END IF
    
    if( len(data1_setNames) == 1 ):
        corrRet = corrRet[0]; #remove the excess list
    #END IF
    
    return corrRet
#END DEF