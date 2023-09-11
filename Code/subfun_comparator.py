#excludes basedon rules

#rule_holder uses a format of 
# {'ref data path':['TEC'], 'data path':[['SuperMAG','SMUs'],['OMNI','Bz GSM']], 'comparison':['elevated|auto,pos only & rate|auto & nan', 'less than|0 & nan'], 'time offset':[122*60, 122*60]}
#data path is how to get to the data you want from the data dict
#comparison is the comparison to do, & adds checks that are done in order, | adds options for the check, , separates the multiple options
#time offset offsets the time so it works good

#dates needed for sunrise/sunset but otherwise it's cool

import numpy as np
from decimal import Decimal
from Code.subfun_sunAlsoRises import sunAlsoRises
from Code.subfun_timeMatch import subfun_timeMatch
from Code.subfun_strfind import strfind
from Code.subfun_date_to_dayNum import subfun_date_to_dayNum

def subfun_comparator(rule_holder, data, dates=None):    
    
    if( not isinstance(rule_holder['ref data path'], (list, tuple)) ): #test if single path with ['dataSet','dataType'] instead of [['dataSet','dataType'],...]
        rule_holder['ref data path'] = [rule_holder['ref data path']]; #wrap it in a list so stuff works
    #END IF
    if( not isinstance(rule_holder['data path'][0], (list, tuple)) ): #test if single path with ['dataSet','dataType'] instead of [['dataSet','dataType'],...]
        rule_holder['data path'] = [rule_holder['data path']]; #wrap it in a list so stuff works
    #END IF
    if( not isinstance(rule_holder['comparison'], (list, tuple)) ): #test if single comparison with 'do compare' instead of ['do compare',...]
        rule_holder['comparison'] = [rule_holder['comparison']]; #wrap it in a list so stuff works
    #END IF
    if( np.isscalar(rule_holder['time offset']) ): #test if scalar
        rule_holder['time offset'] = [rule_holder['time offset']]; #wrap it in a list so stuff works
    #END IF
    
    #get the data reference to use
    data_ref = data; #alias
    for jj in range(0, len(rule_holder['ref data path'])):
        data_ref = data_ref[rule_holder['ref data path'][jj]]; #drill down
    #END FOR i
    
    kCompare = np.zeros(data_ref['time unique'].size, dtype=np.bool_); #preallocate
    for dd in range(0, len(rule_holder['data path'])):
        data2compare = data; #alias
        for jj in range(0,len(rule_holder['data path'][dd])):
            data2compare = data2compare[rule_holder['data path'][dd][jj]]; #drill down
        #END FOR jj
        
        if( data[rule_holder['data path'][dd][0]]['data rate'] != data_ref['data rate'] ):
            #resize!
            data2compare = subfun_timeMatch(data2compare, data[rule_holder['data path'][dd][0]]['time unique']+rule_holder['time offset'][dd], data_ref['time unique'], timeMatch_delta=data_ref['data rate'], FLG_removeNaNs=0, FLG_reportNaNs=False, FLG_useSum=0)[0];
        #END IF
        
        kCompare_method = rule_holder['comparison'][dd].split('&'); #this allows multiple filtering methods to be chained together
        for jk in range(0,len(kCompare_method)):
            
            
            if( strfind(kCompare_method[jk],'|',opt=1) > 0 ): #check for filt info
                kCompare_method_specifics = kCompare_method[jk][kCompare_method[jk].find('|')+1:].split(','); #get any filt specifics
                kCompare_method_specifics = [strang.lower().replace(' ','') for strang in kCompare_method_specifics]; #reduce and standardize
                kCompare_method[jk] = kCompare_method[jk][:kCompare_method[jk].find('|')]; #remove the filt info
            else:
                kCompare_method_specifics = None; #no specifics
            #END IF
            kCompare_method_curr = kCompare_method[jk].lower().replace('-','').replace(' ',''); #get the current method w/o specifics and properly simplified
            
            #--- value comparisons require '<|0' for 0 to be the comparison value ---
            if( ('lessthan' == kCompare_method_curr) or ('<' == kCompare_method_curr) ):
                kCompare_method_specifics[0] = Decimal(kCompare_method_specifics[0]);
                if( (Decimal(kCompare_method_specifics[0]) % 1) == 0 ): #convert to a number to use
                    kCompare_method_specifics[0] = np.int64(kCompare_method_specifics[0]);
                else:
                    kCompare_method_specifics[0] = np.float64(kCompare_method_specifics[0]);
                #END IF
                kCompare = kCompare | (data2compare < kCompare_method_specifics[0]);
                
            elif( ('lessthanorequal' == kCompare_method_curr) or ('<=' == kCompare_method_curr) ):
                kCompare_method_specifics[0] = Decimal(kCompare_method_specifics[0]);
                if( (Decimal(kCompare_method_specifics[0]) % 1) == 0 ): #convert to a number to use
                    kCompare_method_specifics[0] = np.int64(kCompare_method_specifics[0]);
                else:
                    kCompare_method_specifics[0] = np.float64(kCompare_method_specifics[0]);
                #END IF
                kCompare = kCompare | (data2compare <= kCompare_method_specifics[0]);
                
            elif( ('greaterthan' == kCompare_method_curr) or ('>' == kCompare_method_curr) ):
                kCompare_method_specifics[0] = Decimal(kCompare_method_specifics[0]);
                if( (Decimal(kCompare_method_specifics[0]) % 1) == 0 ): #convert to a number to use
                    kCompare_method_specifics[0] = np.int64(kCompare_method_specifics[0]);
                else:
                    kCompare_method_specifics[0] = np.float64(kCompare_method_specifics[0]);
                #END IF
                kCompare = kCompare | (data2compare > kCompare_method_specifics[0]);
                
            elif( ('greaterthanorequal' == kCompare_method_curr) or ('>=' == kCompare_method_curr) ):
                kCompare_method_specifics[0] = Decimal(kCompare_method_specifics[0]);
                if( (Decimal(kCompare_method_specifics[0]) % 1) == 0 ): #convert to a number to use
                    kCompare_method_specifics[0] = np.int64(kCompare_method_specifics[0]);
                else:
                    kCompare_method_specifics[0] = np.float64(kCompare_method_specifics[0]);
                #END IF
                kCompare = kCompare | (data2compare >= kCompare_method_specifics[0]);
                
            #--- elevated requires a distance from the mean value (or auto) (kinda equiv to # of stdev from mean) in 1st spot with optional pos/neg only decorator, e.g., 'elevated|auto,neg only' ---
            elif( 'elevated' in kCompare_method_curr ):               
                kCompare_distFromMedian_noSqr = data2compare - np.nanmedian(data2compare);
                kCompare_distFromMedian = (kCompare_distFromMedian_noSqr)**2; #distance from the median (no sqrt, just leave it at pwr2)
                kCompare_medianDistFromMedian = np.nanmedian(kCompare_distFromMedian); #median of the distance from the median
                
                if( (kCompare_method_specifics[0] == 'auto') | (kCompare_method_specifics[0].replace('.','',1).isdigit() == False) ):
                    #auto creates a comparison value via the ratio between mean and median (e.g., spikes should cause elevated mean and 
                    kCompare_curr = kCompare_distFromMedian > kCompare_medianDistFromMedian*np.nanmean(kCompare_distFromMedian)/kCompare_medianDistFromMedian; #comparison value here is a scalar multiplier, akin to # of standard deviations
                else:
                    kCompare_method_specifics[0] = Decimal(kCompare_method_specifics[0]);
                    if( (Decimal(kCompare_method_specifics[0]) % 1) == 0 ): #convert to a number to use
                        kCompare_method_specifics[0] = np.int64(kCompare_method_specifics[0]);
                    else:
                        kCompare_method_specifics[0] = np.float64(kCompare_method_specifics[0]);
                    #END IF
                    kCompare_curr = kCompare_distFromMedian > kCompare_medianDistFromMedian*kCompare_method_specifics[0]; #comparison value here is a scalar multiplier, akin to # of standard deviations
                #END IF
                
                if( ('posonly' in kCompare_method_specifics[1:]) | ('positiveonly' in kCompare_method_specifics[1:]) | ('+only' in kCompare_method_specifics[1:]) ):
                    kCompare_curr = kCompare_curr & (kCompare_distFromMedian_noSqr > 0); #extra comparator for pos only
                elif( ('negonly' in kCompare_method_specifics[1:]) | ('negativeonly' in kCompare_method_specifics[1:]) | ('-only' in kCompare_method_specifics[1:]) ):
                    kCompare_curr = kCompare_curr & (kCompare_distFromMedian_noSqr < 0); #extra comparator for neg only
                #END IF
                
                kCompare = kCompare | kCompare_curr; #combine in
                
            #--- rate requires a distance from the mean value (or auto) (kinda equiv to # of stdev from mean) in 1st spot with optional pos/neg only decorator, e.g., 'rate|7,+ only' ---
            elif( 'rate' in kCompare_method_curr ):
                kCompare_rate = np.insert(np.diff(data2compare)/data_ref['data rate'], 0, 0); #get the rate, insert 0 at 0 so same size as data2compare
                kCompare_rate_distFromMedian_noSqr = kCompare_rate - np.nanmedian(kCompare_rate);
                kCompare_rate_distFromMedian = (kCompare_rate_distFromMedian_noSqr)**2; #distance from the median (no sqrt, just leave it at pwr2)
                kCompare_rate_medianDistFromMedian = np.nanmedian(kCompare_rate_distFromMedian); #median of the distance from the median
                
                if( (kCompare_method_specifics[0] == 'auto') | (kCompare_method_specifics[0].replace('.','',1).isdigit() == False) ):
                    #auto creates a comparison value via the ratio between mean and median (e.g., spikes should cause elevated mean and 
                    kCompare_curr = kCompare_rate_distFromMedian > kCompare_rate_medianDistFromMedian*np.nanmean(kCompare_rate_distFromMedian)/kCompare_rate_medianDistFromMedian; #comparison value here is a scalar multiplier, akin to # of standard deviations
                else:
                    kCompare_method_specifics[0] = Decimal(kCompare_method_specifics[0]);
                    if( (Decimal(kCompare_method_specifics[0]) % 1) == 0 ): #convert to a number to use
                        kCompare_method_specifics[0] = np.float64(kCompare_method_specifics[0]);
                    else:
                        kCompare_method_specifics[0] = np.int64(kCompare_method_specifics[0]);
                    #END IF
                    kCompare_curr = kCompare_rate_distFromMedian > kCompare_rate_medianDistFromMedian*kCompare_method_specifics[0]; #comparison value here is a scalar multiplier, akin to # of standard deviations
                #END IF
                
                if( ('posonly' in kCompare_method_specifics[1:]) | ('positiveonly' in kCompare_method_specifics[1:]) | ('+only' in kCompare_method_specifics[1:]) ):
                    kCompare_curr = kCompare_curr & (kCompare_rate_distFromMedian_noSqr > 0); #extra comparator for pos only
                elif( ('negonly' in kCompare_method_specifics[1:]) | ('negativeonly' in kCompare_method_specifics[1:]) | ('-only' in kCompare_method_specifics[1:]) ):
                    kCompare_curr = kCompare_curr & (kCompare_rate_distFromMedian_noSqr < 0); #extra comparator for neg only
                #END IF
                
                kCompare = kCompare | kCompare_curr; #combine in
                
            elif( 'nan' in kCompare_method_curr ):
                kCompare = kCompare | np.isnan(data2compare); #tack on nan check\
          
            elif( ('sunrisesunset' in kCompare_method_curr) | ('sunrisensunset' in kCompare_method_curr) | ('sunriseandsunset' in kCompare_method_curr) ):
                #reqs 'sunrise|lat#,long#,timeBeforeToExclude,timeAfterToExclude' - note these can be tacked onto any other analysis, consider adding a "dummy" data if needed since these don't need a data loaded (other than the ref)
                kCompare_method_specifics[0] = str2num(kCompare_method_specifics[0]); #convert em
                kCompare_method_specifics[1] = str2num(kCompare_method_specifics[1]);
                kCompare_method_specifics[2] = str2num(kCompare_method_specifics[2]);
                kCompare_method_specifics[3] = str2num(kCompare_method_specifics[3]);
                [sunRise, sunSet, dateRange_fullPad] = sunAlsoRises(dates['date range full'],kCompare_method_specifics[0],kCompare_method_specifics[1]); #get the sunrise and sunset times for the days at the loc  
                
                dateRange_dayNum_fullPad = subfun_date_to_dayNum(dateRange_fullPad); #convert to dayNum
                sunRise = (sunRise + dateRange_dayNum_fullPad[:,1] - dates['date range zero hr dayNum'][1])*86400; #sec, center around zero hr and convert ot hrs
                sunSet = (sunSet + dateRange_dayNum_fullPad[:,1] - dates['date range zero hr dayNum'][1])*86400; #sec, center around zero hr and convert ot hrs
                sunRise = sunRise[1:]; #remove 1st
                sunSet = sunSet[:-1]; #remove last
                
                for i in range(0, sunRise.size):
                    kCompare = kCompare | ((data_ref['time unique aligned'] < sunRise[i]+kCompare_method_specifics[3]) & (data_ref['time unique aligned'] > sunRise[i]-kCompare_method_specifics[2])); #combine in
                    kCompare = kCompare | ((data_ref['time unique aligned'] < sunSet[i]+kCompare_method_specifics[3]) & (data_ref['time unique aligned'] > sunSet[i]-kCompare_method_specifics[2])); #combine in
                #END FOR i
                
            elif( 'sunrise' in kCompare_method_curr ):
                #reqs 'sunrise|lat#,long#,timeBeforeToExclude,timeAfterToExclude' - note these can be tacked onto any other analysis, consider adding a "dummy" data if needed since these don't need a data loaded (other than the ref)
                kCompare_method_specifics[0] = str2num(kCompare_method_specifics[0]); #convert em
                kCompare_method_specifics[1] = str2num(kCompare_method_specifics[1]);
                kCompare_method_specifics[2] = str2num(kCompare_method_specifics[2]);
                kCompare_method_specifics[3] = str2num(kCompare_method_specifics[3]);
                [sunRise, _, dateRange_fullPad] = sunAlsoRises(dates['date range full'],kCompare_method_specifics[0],kCompare_method_specifics[1]); #get the sunrise and sunset times for the days at the loc  
                
                dateRange_dayNum_fullPad = subfun_date_to_dayNum(dateRange_fullPad); #convert to dayNum
                sunRise = (sunRise + dateRange_dayNum_fullPad[:,1] - dates['date range zero hr dayNum'][1])*3600; #sec, center around zero hr and convert ot hrs
                sunRise = sunRise[1:]; #remove 1st
                
                for i in range(0, sunRise.size):
                    kCompare = kCompare | ((data_ref['time unique aligned'] < sunRise[i]+kCompare_method_specifics[3]) & (data_ref['time unique aligned'] > sunRise[i]-kCompare_method_specifics[2])); #combine in
                #END FOR i
                
            elif( 'sunset' in kCompare_method_curr ):
                #reqs 'sunrise|lat#,long#,timeBeforeToExclude,timeAfterToExclude' - note these can be tacked onto any other analysis, consider adding a "dummy" data if needed since these don't need a data loaded (other than the ref)
                kCompare_method_specifics[0] = str2num(kCompare_method_specifics[0]); #convert em
                kCompare_method_specifics[1] = str2num(kCompare_method_specifics[1]);
                kCompare_method_specifics[2] = str2num(kCompare_method_specifics[2]);
                kCompare_method_specifics[3] = str2num(kCompare_method_specifics[3]);
                [_, sunSet, dateRange_fullPad] = sunAlsoRises(dates['date range full'],kCompare_method_specifics[0],kCompare_method_specifics[1]); #get the sunrise and sunset times for the days at the loc  
                
                dateRange_dayNum_fullPad = subfun_date_to_dayNum(dateRange_fullPad); #convert to dayNum
                sunSet = (sunSet + dateRange_dayNum_fullPad[:,1] - dates['date range zero hr dayNum'][1])*3600; #sec, center around zero hr and convert ot hrs
                sunSet = sunSet[:-1]; #remove last
                
                for i in range(0, sunRise.size):
                    kCompare = kCompare | ((data_ref['time unique aligned'] < sunSet[i]+kCompare_method_specifics[3]) & (data_ref['time unique aligned'] > sunSet[i]-kCompare_method_specifics[2])); #combine in
                #END FOR i
                
            else:
                print('WARNING in subfun_exclude: Bad comparison requested of "'+rule_holder['comparison'][dd]+'", not doing any comparing for that one.');
            #END IF
            
        #END FOR jk
    #END FOR dd
    
    return kCompare
#END DEF

def str2num(strang):
    if( (Decimal(strang) % 1) == 0 ): #convert to a number to use
        numb = np.int64(strang);
    else:
        numb = np.float64(strang);
    #END IF
    return numb
#END DEF