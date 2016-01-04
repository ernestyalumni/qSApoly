## hypergeo.sage
## This is my implementation of hypergeometric polynomials in Sage Math
## The main reference that I'll liberally copy from, and set notation to, is from  
##   An algorithmic proof theory for hypergeometric (ordinary and "q") 
##   multisum/integral identities. Invent. Math., 108 (1992). 575-633
##   by H.S. Wilf and D. Zeilberger
################################################################################
## Copyleft 2016, Ernest Yeung <ernestyalumni@gmail.com>
## 20160104
## 
## This program, along with all its code, is free software; you can redistribute 
## it and/or modify it under the terms of the GNU General Public License as 
## published by the Free Software Foundation; either version 2 of the License, 
## or (at your option) any later version. 
## 
## This program is distributed in the hope that it will be useful, 
## but WITHOUT ANY WARRANTY; without even the implied warranty of 
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
## GNU General Public License for more details.
## 
## linkedin     : ernestyalumni 
## twitter      : ernestyalumni 
## wordpress    : ernestyalumni 
################################################################################

## cf. Wilf and Zeilberger (1992)

R2.<x,y> = PolynomialRing(QQ,2); R
f2 = (x+y)**2
f3 = (x+y)**3
f2.factor()
f2.factor(x)
f3.factor()

binomialpoly = lambda n: (x+y)**n
binomialpolyexpand = lambda n: sum([binomial(n,k)*x**k*y**(n-k) for k in range(n+1)])

# check of binomial theorem as of pp. 577 Wilf and Zeilberger (1992)
binomialpoly(1) == binomialpolyexpand(1)
binomialpoly(2) == binomialpolyexpand(2)
binomialpoly(3) == binomialpolyexpand(3)
binomialpoly(4) == binomialpolyexpand(4)

# rising factorial
a_n = lambda a, n: prod([ (a+k) for k in range(n) ])
q = var("q")
aq_n = lambda a,n: prod([ (1-a*q**k) for k in range(n)])
a = var("a")

