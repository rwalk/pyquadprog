pyquadprog
========

This module provides a Python binding for the Fortran quadprog library of Berwin Turlach. It follows closely the popular R binding:
http://cran.r-project.org/web/packages/quadprog/index.html.

The module can be used to solve quadratic programs of the form
    
    minimize in x:    (-d^Tx + 1/2 x^T D x)
    subject to:    A^T >=b

where the matrix D is symmetric positive definite.


