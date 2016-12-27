# Experimental extension of polymake for toric geometry

## Purpose
This extension for [polymake] (https://www.polymake.org) contains my private
code for experimenting with toric geometry. It is not well documented and many
methods might do different things from what you expect. If you have any
questions, please send me an email.


## Interfacing Singular
There is sample code for building a proof-of-concept type interface to
[Singular] (http://www.singular.uni-kl.de/) in the file
[apps/fulton/rules/singular.rules]
(https://github.com/lkastner/pmExperiments/blob/master/apps/fulton/rules/singular.rules).
This code uses the singular_eval method from the polymake core to execute
Singular commands. The polymake objects get equipped with string properties,
containing the variable names of the Singular objects.

## Cyclic quotient singularities
The file [apps/fulton/rules/cyclic_quotient.rules]
(https://github.com/lkastner/pmExperiments/blob/master/apps/fulton/rules/cyclic_quotient.rules)
contains all algorithms developed in my [thesis]
(http://www.diss.fu-berlin.de/diss/receive/FUDISS_thesis_000000101520) that
were too specific to go into the polymake core. Some of these may be different
than from the print version, since polymake has had numerous updates making
small changes necessary.

## Differential graded Lie algebras
The file [apps/fulton/rules/dgla.rules]
(https://github.com/lkastner/pmExperiments/blob/master/apps/fulton/rules/dgla.rules)
contains code for dealing with differential graded Lie algebras of toric
varieties. I am not really sure what this does, please ask [Klaus Altmann]
(http://www.math.fu-berlin.de/altmann/) for the mathematical details. You can
ask me, if you know the mathematical details and have questions about the
implementation. 

## Deformation theory
There are several files dealing with toric deformation theory:

* apps/fulton/rules/t1.rules
* apps/fulton/rules/deformations.rules
* apps/fulton/rules/downgrading.rules
* apps/fulton/rules/subvarietyInTV\*.rules

You can determine the age of these files using git or estimating from the
cleanliness of the code.

