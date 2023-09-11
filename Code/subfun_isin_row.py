import numpy as np

def isin_row(a, b): #inspired by https://stackoverflow.com/a/67213552/2403531
    #calcs to see if each row of a is in b, yields numpy array that is the same length as a
    z = set(map(tuple, b)); #it is abs essential for performance this is outside of the call below
    return np.asarray([row in z for row in map(tuple, a)])
#END DEF