## Jonespoly_and_kh_sage.py
## Jones polynomial and Khovanov Homology implemented in Sage Math
############################################################################ 
## Copyleft 2016, Ernest Yeung <ernestyalumni@gmail.com>                            
## 20160223
##                                                                               
## This program, along with all its code, is free software; you can redistribute 
## it and/or modify it under the terms of the GNU General Public License as 
## published by the Free Software Foundation; either version 2 of the License, or   
## (at your option) any later version.                                        
##     
## This program is distributed in the hope that it will be useful,               
## but WITHOUT ANY WARRANTY; without even the implied warranty of              
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                 
## GNU General Public License for more details.                             
##                                                         
## linkedin     : ernestyalumni                                                    
## wordpress    : ernestyalumni                                                    
############################################################################ 

import snappy

# list of (p,q) coprime numbers for all torus knots of less than 14 crossings
TKNOTSLIST14 = [(2,1),(2,3),(2,5),(2,7),(2,9),(2,11),(2,13),(3,4),(3,5)]

TKNOTS14SNAP = {}

for K in TKNOTSLIST14:
    TKNOTS14SNAP[K] = snappy.Link('T'+str(K))

# Then one can enter as a tuple for a key to obtain the torus knot desired
# e.g. TKNOTS14SNAP[(2,3)].jones_poly() # -q^4 + q^3 + q

try:
    from kh_scrape import scrape_bat
    Tknots14 = scrape_bat('khTknotsless14')
    Tknots14sage=[sage_eval(line,locals={'x1':t,'x2':q}) for line in Tknots14] # from Khovanov Homology computation

    # normalize to the unknot
    normTknots14sage = [K/Tknots14sage[0] for K in Tknots14sage]


"""
Examples of USAGE:
------------------

from kh_scrape import scrape_bat
Tknots14 = scrape_bat('khTknotsless14')
Tknots14sage = [ sage_eval(line, locals={'x1':t,'x2':q}) for line in Tknots14 ]

# to get a LaTeX printout of all the categorified Jones polynomials
for K in Tknots14sage:
    print latex(K)

normTknots14sage = [K/Tknots14sage[0] for K in Tknots14sage]

# Compare Jones polynomials from SnaPy to Jones polynomials from decategorified Khovanov Homology:
x = var('x')
for i in range(1,len(TKNOTS14SNAP)):
    SnaPyvskh = TKNOTS14SNAP[TKNOTSLIST14[i]].jones_poly().subs(q=x) == normTknots14sage[i].subs(t=-1).subs(q=sqrt(x))
    print bool( SnaPyvskh )

print "If all True, then J(K;t=-1,q) = J(K;q)!"

# Explore what one can do with the categorified Jones polynomial with command dir(): for instance, for the T(2,3) trefoil knot,
dir(Tknots14sage[1])
# one sees the gradient module and we can compute the partials (derivatives) immediately:
Tknots14sage[1].gradient()
# [3*t^2*q^9 + 3*t^2*q^7 + 2*t*q^7 + 2*t*q^5,
#  9*t^3*q^8 + 7*t^3*q^6 + 7*t^2*q^6 + 5*t^2*q^4 + 3*q^2 + 1]

# Also, explore the Newton Polytope with newton_polytope:

Tknots14sage[1].newton_polytope()
# A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices

# this command should print a pretty picture of the Newton Polytope
Tknots14sage[1].newton_polytope().show()
# Launched png viewer for Graphics object consisting of 6 graphics primitives


"""
