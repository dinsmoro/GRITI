#GOAL: Zero-Phase Filter GNSS TEC based on 
#RD on 9/18/2018
#
#INPUT: 
#IN: DATA TO BE FILTERED **NOTE THAT DATA SHOULD HAVE A DAY BEFORE AND AFTER OF DATA TO FILTER CORRECTLY ON DAY EDGES**
#IN: FILTER FREQS - 
#**NOTE THAT FILTER FREQ IS ASSUMED TO BE A PERIOD IN minutes BY DEFAULT**
#low & high get 1 filter freq
#band gets 2 filter freqs
#IN: FILTER TYPE - 
#0 = low pass filter
#1 = high pass filter [DEFAULT]
#2 = band-pass filter, expects 2 numbers for the filter inputs

#OUTPUT: 
#OUT: Filtered data

#options!: 
#0 = filter assuming input is a period (minutes) [DEFAULT]
#1 = filter input is a freqeuncy (Hz, 1/sec)
#

#def GRITI_TEC_subfun_filter(bigData, filtFreq, filtType = 1, filtPeriodOrFreq = 0):



