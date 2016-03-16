/* Indefinite sum of rational function using Gosper's algorithm.
 * GP implementation by R. Stephan, 2004.
 * 2004-07-12: added Eh,gvar by Michael Somos
 * 2004-07-13: removed poldisp code
 *
 * Given a rational function f(x)=P(x)/Q(x), computes rational
 * functions s(x),t(x), such that f(x) = s(x+1)-s(x)+t(x), with the
 * denominator of t minimal. Only univariate case is supported.
 * 
 * See end of file for bibliography. The method was described by 
 * Paule and Pirastu [Pir92,Pau95,Pir96]. See Petkovsek/Wilf/Zeilberger's
 * book [PWZ97] for an extensive description of Gosper's algorithm
 * and more references and ideas.
 */

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
/* variable of polynomial/rational function otherwise default x */
{gvar(x)=if(type(x)=="t_POL",variable(x),
         if(type(x)=="t_RFRAC",variable(denominator(x)),'x))}
/* increment expression variable by h */
Eh(x,v,h=1)= if(!v,v=gvar(x));if(!v,x,subst(x,v,v+h))
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

gsum(rf,x)=
{ local(f,p,wp,h,g,u,v,D,ad,bd,m,rhs,t,num1,num2);

  if(!x,x = gvar(rf));
  /* handle whole part of fraction */
  wp = numerator(rf) \ denominator(rf);
  f = rf-wp;
  if (wp,
    h = 0; g = poldegree(wp); p = bernvec(g\2+1);
    while (g > 0,
      v = vector(g+2); v[1] = 1; v[2] = -(g+1) / 2;
      forstep(k = 3, g+1, 2, v[k] = binomial(g+1,k-1) * p[ (k+1)\2 ]);
      h += pollead(wp) / (g+1) * Pol(v,x);
      wp -= pollead(wp) * x^g; g = poldegree(wp)
    );
    wp = h + wp*x;
  );
  
  /* determine optimal denominators */
  p = denominator(f);
  v = gosper (p/Eh(p),x);
  u = p*v[1]/Eh(v[3],x,-1);
  v = qplus (v[2],x);
  
  /* set up numerator equation, solve it */
  D = lcm (p, lcm (u, lcm (v, Eh(u,x))));
  ad = poldegree(u)-1;
  bd = poldegree(v)-1;
  m = matrix (ad+bd+2,ad+bd+2);
  rhs = vectorv (ad+bd+2);
  t = (D/p)*numerator(f);
  for (k=1, ad+bd+2, 
    rhs[k] = -polcoeff(t,k-1));
  t = D/v;
  for (k=0, bd,
    for (l=0, poldegree(t),
      m[k+l+1,ad+k+2] = polcoeff(t,l)));
  t = D/u;
  for (k=0, ad,
    for (l=0, poldegree(t),
      m[k+l+1,k+1] = -polcoeff(t,l)));
  t = D/Eh(u,x);
  for (k=0, ad,
    for (l=0, k,
      bin = binomial(k,l);
      for (ll=0, poldegree(t),
        m[l+ll+1,k+1] += bin*polcoeff(t,ll))));
  t = matsolve (m, rhs);

  /* build up solution */
  num1 = 0; num2 = 0;
  for (k=0, ad, num1 += -t[k+1] * x^k);
  for (k=0, bd, num2 += -t[k+ad+2] * x^k);
/*print(if(subst(wp+num1/u,x,x+1)-wp-num1/u+num2/v-rf==0,"OK","FAIL"));*/
  [wp+num1/u, num2/v]
}

qplus(q,x)=
{ local(d,g,h,qp);

  d = disp (q,q,x);
  g = q;
  qp = 1;
  forstep (j = d, 0, -1,
    h = part_of (Eh(g,x,-j), g);
    g = g/h;
    qp = qp * Eh(h,x,j);
  ); qp
}

gosper(pol,x)=
{ local(facq,facr,l,p,q,r,tp,v,lc,d);

  p = 1;
  q = numerator(pol);
  r = denominator(pol);
  v = dispset (q,r,x);
  for (k=1, length(v), 
    d = gcd (q, Eh(r,x,v[k]));
    q = q/d;
    r = r / Eh(d,x,-v[k]);
    p = p * prod(kk=1, v[k], Eh(d,x,-kk));
    );
  v = vector(3);
  v[1] = p;
  v[2] = q;
  v[3] = r;
  v
}

part_of(p1,q1)=
{local(P,qq,p,q);

  p=p1; q=q1;
  if (poldegree(p)==0||poldegree(q)==0, return(1));
  P = 1;
  while (poldegree(gcd(p,q)), 
    qq = q/gcd(p,q);
    P = P*gcd(p,q);
    q = qq
  ); P
}

dispset(p,q,x)=
{ local(d,l,r);

  r = polresultant (Eh(q,x,hh), p/gcd(p,p'), x);
  r = r/content(r);
  r = r/gcd(r,hh^poldegree(r));
  d = disp (p,q,x);
  if (d<=0, return(listcreate(0)));
  l = listcreate (d);
  for (k=1, d-1,
    if (subst (r, hh, k) == 0,
      listput (l, k)
    )
  );
  listput (l, d);
  l
}

disp(p,q,x)=
{ local(d,r,g1,g2,m,n);

  g1 = GFF(p,x);
  g2 = GFF(q,x);
  m = prod(k=1,length(g1),g1[k]);
  n = Eh(prod(k=1,length(g2),Eh(g2[k],x,1-k)), x, hh);
  r = polresultant (m/gcd(m,m'), n/gcd(n,n'), x);
  r = r/content(r);
  r = r/gcd(r,hh^poldegree(r));
  d = divisors(abs(polcoeff(r,0)));
  forstep (k=length(d),1,-1,
    if (subst(r,hh,d[k])==0,
      return(d[k])
    )
  );
  0
}

GFF(pol,x)=
{
  local(i,g,pr,p);

  i = 1;
  g = gcd(pol, Eh(pol,x));
  pr = pol / Eh(g,x,-1);
  p = vector(1000);
  p[i] = gcd (pol/g, pr);
  while (poldegree(g),
    i = i + 1;
    pol = g;
    pr = pr / p[i-1];
    g = Eh (pol/pr,x);
    p[i] = gcd (pol/g, pr);
    );
  vector (i, n, if(n==1,g*p[1],p[n]))
}

/*******************************************************************

\begin{thebibliography}{A99}
\bibitem[Gra94]{b-gkp} 
R. L. Graham, D. E. Knuth and O. Patashnik, \emph{
Concrete Mathematics}, 2nd ed., Addison-Wesley, 1994
\bibitem[Pau95]{b-paule}
P. Paule, \emph{Greatest factorial factorization and symbolic summation},
J. Symb. Comp., 20 (1995) 235-268.
\bibitem[PWZ97]{b-pwz}
M. Petkov\v sek, H. S. Wilf and D. Zeilberger, $A=B$, A.~K.~Peters,
Wellesley, MA, 1996
\bibitem[Pir92]{b-pirm}
R. Pirastu, \emph{Algorithmen zur Summation rationaler
Funktionen}, Master's Thesis, Erlangen 1992
\bibitem[Pir96]{b-pirt}
R. Pirastu, \emph{On combinatorial identities: symbolic summation and
umbral calculus}, PhD Thesis, Linz 1996
\end{thebibliography}
%\verb+http://www.risc.uni-linz.ac.at/research/combinat/risc/publications/+

*********************************************************************/
/* Examples: */

gsum((-x^2-3*x-3)/(x^4+2*x^3-3*x^2-4*x+2))
gsum((4*x+5)/(2*(4*x+1)*(2*x+3)))
gsum(x*(x+3)/(x-1)/(x+4))
gsum((x+1)^2*(x+3)/(x*(x+2)^2))
gsum((-3*x-5)/(x^5+23/3*x^4+199/9*x^3+793/27*x^2+1400/81*x+784/243))
gsum((x^2-3*x+1)/((x-1)^2*x^3*(x+3)*(x*x+1)*(x^2+4*x+5)^2))
gsum(((x-3)*(x-2)^2*(x+2)*(x+5)^2)/((x-4)*(x+1)^3*(x+3)^2))
gsum(1/(x^3*(x+2)^2*(x+3)*(x^2+1)*(x^2+4*x+5)))


