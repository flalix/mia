'''
Created on 10/02/2015

@author: flavio

in ~/Downloads/PyCogent-1.5.3$
python setup.py build

paper: http://genomebiology.com/content/pdf/gb-2007-8-8-r171.pdf
help: http://pycogent.org/
source: http://sourceforge.net/projects/pycogent/?source=typ_redirect

'''
import math, Bio

from cogent.evolve.models import GTR

a = math.sqrt(2)
sub_mod = GTR(with_rate=True, distribution='gamma')
print sub_mod

    