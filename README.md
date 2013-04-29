robustFourierSeries
===================

Python 2.7 module for robust estimation of a Fourier Series, and methods for estimating periodic signals in data.


This module contains methods to calculate
the periodogram for time series data by
estimating the coefficents of a discreat 
Fourier Series.  The g-statistic is also
estimated and tested to identify the periodic 
component that is statistically significant.
This uses the robust non-uniform sampling 
algorithm described in Ahdesmai 2007
(Robust regression for periodicity detection in non-uniformly 
sampled time-course gene expression data) 
and the genral pereto distribution estimates
from Knijnenburg 2009
(Fewer permutations, more accurate P-values)
Also uses the robust regression form 
scikits statsmodels
(http://statsmodels.sourceforge.net/install.html)

The basic method is estSigGStat; however,
many other methods are avalible.

Please excuse the spelling and poor commenting :)

     Copyright (C) 2003-2012 Institute for Systems Biology
                             Seattle, Washington, USA.
 
     This library is free software; you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public
     License as published by the Free Software Foundation; either
     version 2.1 of the License, or (at your option) any later version.
 
     This library is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
     Lesser General Public License for more details.
 
     You should have received a copy of the GNU Lesser General Public
     License along with this library; if not, write to the Free Software
     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA

