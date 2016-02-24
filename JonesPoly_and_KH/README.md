# JonesPoly_and_KH
## The Jones Polynomial and Khovanov Homology

Remember that detailed and verbose instructions and descriptions are in `qsuperApoly.pdf` (`qsuperApoly.tex` is the TeX file) in the `LateX_and_pdf` directory of this `qSApoly` repository.  Otherwise, the *gist* is this:

* Be sure to have `knotkit` installed and `SnaPy` installed in `Sage Math`.  
* If you want to *confirm for yourself* (i.e. recreate the "categorified" Jones polynomial, computed from Khovanov Homology), start with files
  - `kk.cpp` (thank you to Prof. N. Dunfield for this modified version of `kk.cpp` which you copy over the original `kk.cpp` of `knotkit`, and then run `make`, and run the executable `kk` in order to have the output of the Khovanov Homology calculation be a text (string)
  - `kk_kh_batch.py` : run this to make all the "categorified" Jones polynomials for Torus knots of less than 14 crossings
* Otherwise, the categorified Jones polynomials from Khovanov Homology *J(K;t,q)* using `knotkit` has been saved as text files in the following files (and this repository serves as a place to store these polynomials and agree upon their conventions):
  - `kh_T0201`
  - `kh_T0203`
  - `kh_T0205`
  - `kh_T0207`
  - `kh_T0209`
  - `kh_T0211`
  - `kh_T0304`
  - `kh_T0305`
  - `khTknotsless14` (this file has all the categorified Jones polynomials from Khovanov Homology *J(K;t,q)* for torus knots of less than 14 crossings)

EY: 20160224 `kh_T0213` is missing because my MacBook Air stalled

So if we agree upon these polynomials, then
* Open up *Sage Math* and run `kh_scrape.py` in it by importing in the 2 functions in the file, 
  - `scrape_file` (to scrape a file containing a single polynomial, like `kh_T0203`
  - `scrape_bat` (to scrape a file containing many polynomials, like `khTknotsless14`
* Then take a look at `Jonespoly_and_kh_sage.py` which you can load in Sage Math and follow the examples in there to confirm that the categorified Jones polynomials *J(K;t,q)* indeed reproduce the Jones polynomials from *SnaPy*, by setting t=-1
* It would be nice to explore Sage Math modules (functions) that are available for the polynomials, such as `newton_polytope` 

### Example Output
In Sage Math, in the working directory including the files in this directory `JonesPoly_and_KH` of this repository,
```
sage: from kh_scrape import scrape_bat
sage: Tknots14 = scrape_bat('khTknotsless14')
sage: Tknots14sage = [ sage_eval(line, locals={'x1':t,'x2':q}) for line in Tknots14 ]
sage: latex(Tknots14sage[1]) # J(T(2,3);t,q)
t^{3} q^{9} + t^{3} q^{7} + t^{2} q^{7} + t^{2} q^{5} + q^{3} + q
sage: latex(Tknots14sage[2]) # J(T(2,5);t,q)
t^{5} q^{15} + t^{5} q^{13} + t^{4} q^{13} + t^{4} q^{11} + t^{3} q^{11} + t^{3} q^{9} + t^{2} q^{9} + t^{2} q^{7} + q^{5} + q^{3}
```