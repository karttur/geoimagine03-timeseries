'''
Created on 17 Feb 2021

@author: thomasgumbricht
'''

# Standard library imports

# Third party imports

import numpy as np
 
# Package application imports

from geoimagine.ktnumba import InterpolateLinearNaNNumba

class TimeSeriesNumba:
    ''' translator between Time Series and Numba 
    '''
    
    def __init__(self):
        '''Empty call to access the functions
        '''
        pass
    
    def _FillAlongAxis(self, ts, validfraction=0.5):
        '''Linear interpolation of NaN, calls Numba function
        '''
        
        if np.all(np.isnan(ts)):
                        
            return self.dstNullDarr[self.activeComp]
        
        if np.isnan(np.sum(ts)):
            
            #non_nans = (~np.isnan(ts)).sum()
            non_nans = np.count_nonzero(~np.isnan(ts))
            
            nans = np.count_nonzero(np.isnan(ts))
            
            if float(non_nans)/ts.shape[0] < validfraction:
                                
                return self.dstNullDarr[self.activeComp]
            
            avg = np.nanmean(ts)
            
            if np.isnan(ts[0]):
                
                ts[0] = avg
                
            if np.isnan(ts[ts.shape[0]-1]):
                
                ts[ts.shape[0]-1] = avg
                
            ts = InterpolateLinearNaNNumba(ts)
            
            return ts
        
        else:
            
            return ts