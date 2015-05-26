[nls2](http://cran.r-project.org/web/packages/nls2/index.html) is an [R](http://www.r-project.org) package that adds the "brute-force" and related algorithms as well as multiple starting values to the R [nls](http://finzi.psych.upenn.edu/R/library/stats/html/nls.html) function.  `nls2` is free software licensed under the GPL and available from [CRAN](http://cran.r-project.org/web/packages/nls2/index.html).  It provides a function, `nls2`, which is a superset of the R `nls` function which it, in turn, calls.

## News ##

**We will now be using the newly created [sqldf discussion group](http://groups.google.com/group/sqldf) for discussion of nls2 as well.**

March 8, 2013.  nls2 0.2 is now on [CRAN](http://cran.r-project.org/web/packages/nls2/index.html).  There are two new algorithms -- `"plinear-brute"` and `"plinear-random"`.  They perform brute force or random search on the non-linear parameters setting the linear parameters to the optimum values given the non-linear ones.  There are also minor documentation updates and improvements.

August 22, 2010. nls2 version 0.1-3 uploaded to CRAN.  See [NEWS](http://cran.r-project.org/web/packages/nls2/NEWS) file.

## Installation ##
As with any R package on CRAN, `nls2` can be installed from within R using
```
install.packages("nls2")
```

## Help ##
After installation, the `nls2` help page can be accessed from within R via
```
library(nls2)
help(nls2)
```
or it can be accessed online [here](http://cran.r-project.org/web/packages/nls2/nls2.pdf).

## Examples ##
There are examples:
  * at the end of the `nls2` [help page](http://cran.r-project.org/web/packages/nls2/nls2.pdf).
  * posted on r-help:  see this [question](http://tolstoy.newcastle.edu.au/R/e4/help/08/05/12190.html) and its [answer](http://tolstoy.newcastle.edu.au/R/e4/help/08/05/12355.html).
  * on [pages 28 and 29](http://books.google.com/books?id=9Aq5k0hZLykC&pg=PA28&lpg=PA28&dq=Ritz+Streibig+nls2&source=bl&ots=j7rL2VKbj6&sig=2fec6I2i9Hrw_Pf0PWouVf95d7U#PPA29,M1) of the book [Nonlinear Regression with R](http://books.google.ca/books?id=o1KtNAEACAAJ&dq=nonlinear+regression+with+R&hl=en&sa=X&ei=ru05UYaVJ6qG2gXDrIDgBQ&redir_esc=y) by Ritz & Streibig.

## Citation ##
The citation for this package can be obtained by issuing this command from within R:
```
citation("nls2")
```

(Note that this package is unrelated to the [software of the same name](http://w3.jouy.inra.fr/unites/miaj/public/AB/nls2/welcome.html) associated with the book Statistical Tools for Nonlinear Regression by Huet et al.)