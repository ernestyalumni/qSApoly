## hypergeo.py
## This is my implementation of hypergeometric polynomials in sympy
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
import sympy
from sympy import binomial, factorial, floor, rf, Product, Symbol, symbols, Sum # rf RisingFactorial

from sympy.utilities.lambdify import lambdify, implemented_function
from sympy import Function
from sympy import Rational as Rat # treat integers as "sympy" integers

from sympy.abc import alpha, beta

from sympy.polys.orthopolys import (
    hermite_poly,
    laguerre_poly,
    legendre_poly,
    jacobi_poly
)

## cf. Wilf and Zeilberger (1992)

x, y = symbols('x y')
f2 = (x+y)**2
f2.expand()
f3 = (x+y)**3
f3.expand()

n = Symbol('n',positive=True, integer=True)
k = Symbol('k',integer=True)
#binomialpoly = implemented_function(Function('(x+y)^n'), lambda n: (x+y)**n)
#binomialpoly_n = lambdify(n, binomialpoly(n))
binomialpoly = (x+y)**n
binomialpolyexpand = Sum( binomial(n,k)*x**k*y**(n-k),(k,0,n))

binomialpoly.subs(n,1)==binomialpolyexpand.subs(n,1).doit() # True
binomialpoly.subs(n,2).expand()==binomialpolyexpand.subs(n,2).doit() # True
binomialpoly.subs(n,3).expand()==binomialpolyexpand.subs(n,3).doit() # True
binomialpoly.subs(n,4).expand()==binomialpolyexpand.subs(n,4).doit() # True

# rising factorial
a,q = symbols('a q')
#a_n = Product(a+k,(k,0,n-1))
a_n  = rf(a,n)
aq_n = Product(Rat(1)-a*q**k,(k,0,n-1))
# e.g. aq_n.subs(n,3).doit()

# Hermite polynomials

HermiteF = factorial(n)*(Rat(-1))**k*(Rat(2)*x)**(n-Rat(2)*k)/ ( factorial(n-Rat(2)*k)*factorial(k))
HermiteP = Sum( HermiteF,(k,0,floor(n/2)) )
hermite_poly(1,x) == HermiteP.subs(n,1).doit() # True
hermite_poly(2,x) == HermiteP.subs(n,2).doit() # True
hermite_poly(3,x) == HermiteP.subs(n,3).doit() # True
hermite_poly(4,x) == HermiteP.subs(n,4).doit() # True

(HermiteF.subs(n,n+1)/HermiteF).simplify() # 2*x*(n + 1)/(-2*k + n + 1)
(HermiteF.subs(k,k+1)/HermiteF).simplify() # -(2*k - n)*(2*k - n + 1)/(4*x**2*(k + 1))

LaguerreF = binomial(n+alpha,n-k)*(-x)**k/factorial(k)
LaguerreP = Sum( LaguerreF, (k,0,n))

for alphavar in range(4):
    for nvar in range(4):
        print alphavar, nvar, LaguerreP.subs(n,nvar).subs(alpha,alphavar).doit(), LaguerreP.subs(n,nvar).subs(alpha,alphavar).doit() == laguerre_poly(nvar,x,alpha=alphavar) # True

(LaguerreF.subs(n,n+1)/LaguerreF).simplify() # (alpha + n + 1)/(-k + n + 1)
(LaguerreF.subs(k,k+1)/LaguerreF).simplify() # x*(k - n)/((k + 1)*(alpha + k + 1))

LegendreF = Rat(1)/(Rat(2)**n)*(binomial(n,k))**2*(x-Rat(1))**k*(x+Rat(1))**(n-k)
LegendreP = Sum( LegendreF,(k,0,n))

for nvar in range(4):
    legendre_poly(nvar,x) == LegendreP.subs(n,nvar).doit().expand()  # True

(LegendreF.subs(n,n+1)/LegendreF).simplify() # (n + 1)**2*(x + 1)/(2*(k - n - 1)**2)
(LegendreF.subs(k,k+1)/LegendreF).simplify() # (k - n)**2*(x - 1)/((k + 1)**2*(x + 1))

JacobiF = Rat(1)/(Rat(2)**n)*binomial(n+alpha,n-k)*binomial(n+beta,k)*(x-Rat(1))**k*(x+Rat(1))**(n-k)
JacobiP = Sum( JacobiF, (k,0,n))
for alphvar in range(5):
    for betvar in range(5):
        for nvar in range(5):
            jacobi_poly(nvar,alphvar,betvar,x)==JacobiP.subs(n,nvar).subs(alpha,alphvar).subs(beta,betvar).doit().expand() # True

(JacobiF.subs(n,n+1)/JacobiF).simplify() # (x + 1)*(alpha + n + 1)*(beta + n + 1)/(2*(k - n - 1)*(-beta + k - n - 1))
(JacobiF.subs(k,k+1)/JacobiF).simplify() # (k - n)*(x - 1)*(-beta + k - n)/((k + 1)*(x + 1)*(alpha + k + 1))


