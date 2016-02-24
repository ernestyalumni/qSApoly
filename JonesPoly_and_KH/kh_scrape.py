## kh_scrape.py
## From the output of knotkit for kh (Khovanov homology), we "scrape" or 
## "clean" text file outputs into useful formats for sympy and Sage Math
################################################################################# 
## Copyleft 2015, Ernest Yeung <ernestyalumni@gmail.com>                            
## 20160222
##                                                                               
## This program, along with all its code, is free software; you can redistribute 
## it and/or modify it under the terms of the GNU General Public License as 
## published by the Free Software Foundation; either version 2 of the License, or   
## (at your option) any later version.                                        
##     
## linkedin     : ernestyalumni                                                    
## wordpress    : ernestyalumni                                                    
#################################################################################  

def scrape_file(filename="kh_T0201"):
    """
    scrape_file = scrape_file(filename="kh_T0201")
    
    Processes or "cleans" the output from Knotkit (kk) and "khplainout" option.
    Returns a dictionary with the filename for 'kind' of knot it is, the 
    polynomial (which represents the topological invariant) as a string, under 
    'poly' and the rank, under 'rank' of the invariant.  

    INPUT/Parameters
    ----------------
    filename='kh_T0201' : <string>
    File name of file with the output from the kk program 
    """
    with open(filename,"r") as khfile:
        khstr=khfile.readlines()
        khstr=[string.replace('\n','').strip() for string in khstr]
    
    # clean the polynomial string in khstr[0]
    polystr=[term.strip() for term in khstr[0].split('+')]
    
    terms=[]
    for term in polystr:
        if term.find('x2') == -1: # check that 'x2' substring is at least there
            terms.append(term)
        elif term.find('x2') == 0 or term.find('x2') == 1: # no 'x1' terms case
            terms.append(term)
        else:
            x2ind = term.index('x2') # index or "location" in string of 'x2'
            termf = term[0:x2ind]+'*'+term[x2ind:]
            terms.append(termf)
    polystr = ' + '.join(terms)

    knot = {'kind':filename,'poly':polystr,'rank':int(khstr[1])}
    return knot

def scrape_bat(filename="khTknotsless14"):
    """
    scrape_bat = (filename='khTknotsless14')
    "Scrape" or "clean" a batch of Khovanov Homology computed Jones polynomials from 
    knotkit to Sage Math.

    Returns a list of text (strings) that are the polynomials in the original file:
    remember that you'll have to keep track manually of which polynomial corresponds
    to which Torus knot.
    """
    f = open("khTknotsless14",'r')
    khbat = [line.replace('\n','') for line in f.readlines()]
    
    cleanedbat=[]
    for line in khbat:
        polystr=[term.strip() for term in line.split('+')]

        terms=[]
        for term in polystr:
            if term.find('x2') == -1: # check that 'x2' substring is at least there
                terms.append(term)
            elif term.find('x2') == 0 or term.find('x2') == 1: # no 'x1' terms case
                terms.append(term)
            else:
                x2ind = term.index('x2') # index or "location" in string of 'x2'
                termf = term[0:x2ind]+'*'+term[x2ind:]
                terms.append(termf)
        polystr = ' + '.join(terms)
        cleanedbat.append(polystr)
    return cleanedbat
    

#####################################################################################
##### EXAMPLEs of USAGE in Sage Math
#####################################################################################
"""
from kh_scrape import scrape_file
Qkh.<t,q> = PolynomialRing(RationalField(),2)
T0203stf = scrape_file('kh_T0203')
T0203poly = sage_eval( T0203stf['poly'], locals={'x1':t,'x2':q} )
T0203poly.substitute(t=-1) # for some reason, subs doesn't work in this case

from kh_scrape import scrape_bat
Tknots14sage=[sage_eval(line,locals={'x1':t,'x2':q}) for line in Tknots14]

"""
