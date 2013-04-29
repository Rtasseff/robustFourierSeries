"""
Robust statistical models.

This code was taken directly from scikits statsmodels.  
The robustFourierSeries code is dependent on this module in 
statsmodels, but to avoid needing the whole program I simply 
pulled out the relevent parts. 
All of the code from statsmodels was originaly published under
the  Modified (3-clause) BSD license at http://opensource.org/licenses/BSD-3-Clause.


"""
import norms
from .scale import mad, stand_mad, Huber, HuberScale, hubers_scale

#from scikits.statsmodels import NoseWrapper as Tester
#test = Tester().test
