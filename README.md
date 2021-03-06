----

This repository is no longer supported. It has been merged with the regular one from the original team of DiceKriging: https://github.com/DiceKrigingClub/DiceKriging , and all features are translated there.

You can also access the website dedicated to DiceKriging: https://dicekrigingclub.github.io/www/

----


Linux & OSX [![Build Status](https://travis-ci.org/IRSN/DiceKriging.png)](https://travis-ci.org/IRSN/DiceKriging)
Windows [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/DiceKrigingClub/DiceKriging-l4eun?branch=master&svg=true)](https://ci.appveyor.com/project/DiceKrigingClub/DiceKriging-l4eun)

[![codecov](https://codecov.io/gh/IRSN/DiceKriging/branch/master/graph/badge.svg)](https://codecov.io/gh/IRSN/DiceKriging)

# DiceKriging: Kriging methods for computer experiments

This repository is a fork of regular DiceKriging sources (available at http://cran.r-project.org/web/packages/DiceKriging).
It contains some fixes and supplement for testing purpose:

 * control random seed for genoud loglik optimization (kmEstimate.R)
 * allow arbitrary optim.method function to be passed (expected to work like optim)
 * add scaling derivatives to be used within DiceOptim/EGO (covStruct_Scaling.R, covStruct_Scaling_Affine.R)
 * big improvement of scaling code for cpu efficiency

Installation
------------

You can install the latest version of the code:

  * using the `devtools` R package:
```
# Install devtools, if you haven't already.
install.packages("devtools")

devtools::install_github("DiceKriging", "IRSN")
```
  * using pre-built packages:
    * Windows: https://github.com/IRSN/DiceKriging/releases/download/windows/DiceKriging_1.5.55.zip
    * Linux: https://github.com/IRSN/DiceKriging/releases/download/linux/DiceKriging_1.5.55.tar.gz
    * OSX: https://github.com/IRSN/DiceKriging/releases/download/osx/DiceKriging_1.5.55.tgz

![Analytics](https://ga-beacon.appspot.com/UA-109580-20/DiceKriging)
