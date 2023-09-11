"""
Makes numbers look nice
"""
from numpy import ndarray

def textNice(number):
    if( isinstance(number,str) ):
        text = textNice_cleanser(number); #already a string
    else:
        if( isinstance(number,ndarray) | isinstance(number,list) | isinstance(number,tuple) ):
            text = '';
            for i in range(0,len(number)):
                text = text + textNice_cleanser(str(number[i])) + ', ';
            #END FOR i
            text = text.rstrip(', '); #remove the last one
        else:
            text = textNice_cleanser(str(number)); #convert number to a string
        #END IF
    #END IF
    
    return text #success
#END DEF

def textNice_cleanser(text):
    
    #=====decimal-required stuff=====
    decLoc = text.find('.');
    if( decLoc > -1 ):
        #---remove trailing zeros---
        if( text[-1] == '0' ):
            text = text.rstrip('0'); #only strip 0 if there's a decimal (prevents removing like 10 -> 1)
        #END IF
        
        #---remove repeating values---
        if( text[-1].isdigit() & text[-2].isdigit() & (len(text[decLoc+1:-1]) > 2) ):
            lastNsecondLast = text[-1] == text[-2];
            gr8trThan4 = int(text[-1]) > 4;
            if( ((text[-1] == str(int(text[-2])+1)) & gr8trThan4) | lastNsecondLast ): #catch for possible repeats
                if( lastNsecondLast ):
                    if( not gr8trThan4 ):
                        ender = ''; #nothing to add on since something like .54444 should become .544
                    else:
                        ender = str(int(text[-1])+1); #weird string that was like 0.76666666666 that should become 0.7667 so this helps
                    #END IF
                else:
                    ender = text[-1];
                #END IF
                cntr = -2;
                # cntr_same = 1;
                FLG_same = True;
                decBit = text[decLoc+1:-1];
                decBit_len = len(decBit);
                sameChar = text[-2];
                while( FLG_same ):
                    if( (sameChar == decBit[cntr]) & (decBit_len > abs(cntr)) ):
                        # cntr_same += 1; #increment count
                        cntr -= 1; #decriment
                    else:
                        FLG_same = False; #leaving
                    #END IF
                #END WHILE
                if( cntr < -3 ): #using negative counting to enable 1 cntr to do it all (want to ensure at least 3 repeating values so that it can be condensed to 2 repeating values)
                    text = text[:decLoc+1]+decBit[:cntr+3]+ender; #rebuild it
                #END IF
            #END IF
        #END IF

        text = text.rstrip('.').replace('. ',' ').replace('.]',']').replace('.0]',']').replace('.0,',','); #always want to strip a . on the end of a number b/c useless info if 10. for plotting
    #END IF
    
    #=====comma-required stuff=====
    commaLoc = text.find(',');
    if( commaLoc == -1 ):
        text = text.lstrip('[').rstrip(']'); #removes the situation [33] where 1 number in an array is printed
    #END IF
    
    return text
#END DEF