#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
This module solves quadratic programs of the form
    
    minimize in x:    (-d^Tx + 1/2 x^T D x)
    subject to:    A^T >=b

where the matrix D is symmetric positive definite using the methods
of Goldfarb and Idnani (1982, 1983).  The solvers are compiled Fortran routines
implemented by Berwin A. Turlach originally for the S-language.  These solvers
were ported to the popular R package "quadprog" by Andreas Weingessel.

The module design follows closely the R implementation.

"""
import numpy as np
from pyquadprog import aind, qpgen1, qpgen1_inplace, qpgen2, qpgen2_inplace

class Solution:
    '''A solution class for the result of a call to a quadprog solver.
    
    Attributes:

        solution: solution to the quadratic program

        crval: scaler, value of the objetive function at the solution

        unconstrained solution:  nx1, global minimizer

        iact: qx1, constraints active in the final iteration

        nact: scalar, number of constraints active in the final iteration

        iterations: 2x1, fist component is number of iterations; 
            the second the number of constraints deleted after they became active'''

    def __init__(self, n, q):
        self.iact = np.zeros(q)
        self.nact = 0
        self.solution = np.zeros(n)
        self.unconstrained = np.zeros(n)
        self.Lagrangian = np.zeros(q)
        self.crval = 0
        r = min(n,q)
        self.work = np.zeros(int((2*n+r*(r+5.0)/2+2*q+1)))
        self.iterations = np.zeros(2)
    
    def __str__(self):
        s = "Solution = " + self.solution.__str__()
        return s

def solveCompactFormQP(Dmat, dvec, Amat, Aind, bvec, meq=0, factorized=False, destroyOK=False):
    ''' Solve a quadratic program using the active set methods of Goldfarb and Idnani (1982, 1983)
    
    This methods solves a quadratic program of the form:
    
    minimize in x:    (-d^Tx + 1/2 x^T D x)
    subject to:    A^T >=b

    where A is stored in a memory efficient form using the matrices Amat and Aind.

    Arguments:

    Dmat -- matrix D appearing in the quadratic form; when factorized = True, Dmat = R^{-1}
         where D=R^T
    
    dvec -- vector d appearing in the quadratic form
    
    Amat -- matrix with amat[k,i] equal to the k-th non-zero element in column i of A.

    Aind -- matrix with:
        Aind[i,i] equal to the number of non-zero elements in column i of A
        Aind[k,i] for k>1 equal to j if the (k-1)st non-zero element of A[:,i]
        is A(i,j).

    bvec -- the vector of constants in the inquality constraint

    meq  -- the first meq constraints of the inequality constraint will be treated as equalities

    factorized -- when factorized is true, Dmat=R^{-1} where D = R^TR.
    
    destroyOK -- if True, Amat, dvec, Aind, bvec will be destoryed on exit but computation
    may be performed more efficiently.
    
    '''

    n = Dmat.shape[0]
    anrow, q = Amat.shape
    if bvec is None:
        bvec = np.zeros(q)
    if n != Dmat.shape[1]:
        raise ValueError("Dmat is not symmetric.")
    if n != len(dvec):
        raise ValueError("Dmat is not compatible with dvec.")
    if meq>q or meq<0:
        raise ValueError("meq is invalid.")
    if anrow+1 != Aind.shape[0] or q!= Aind.shape[1] or q!= len(bvec):
        raise ValueError("Amat, Aind, and bvec are incompatible!")

    # check the indices for the input format.
    Aindok = aind(Aind, q, n)
    if Aindok==False:
        raise ValueError("Aind contains illegal indexes.")

    sol = Solution(n,q)
    err = factorized
    
    if destroyOK:
        qpgen1_inplace(Dmat,
            dvec,
            Amat,
            sol.solution,
            sol.crval,
            Amat,
            Aind,
            bvec,
            q,
            meq,
            sol.iact,
            sol.nact,
            sol.iterations,
            sol.work,
            err)        
    else:
        sol.solution, sol.crval, sol.iact, sol.nact, sol.iter, sol.work = qpgen1( 
         Dmat,
         dvec,
         Amat,
         Aind,
         bvec,
         meq,
         err)
    if err == 1:
        ValueError("constraints are inconsistent, no solution.")
    if err == 2:
        ValueError("matrix D in quadratic function is not positive definite!")
    return(sol)

def solveQP(Dmat, dvec, Amat, bvec=None, meq=0, factorized=False, destroyOK=False):
    ''' Solve a quadratic program using the active set methods of Goldfarb and Idnani (1982, 1983)
    
    This methods solves a quadratic program of the form:
    
    minimize in x:    (-d^Tx + 1/2 x^T D x)
    subject to:    A^T >=b

    Arguments:

    Dmat -- The matrix D appearing in the quadratic form; when factorized = True, Dmat = R^{-1}
         where D=R^T
    
    dvec -- vector appearing in the quadratic form
    
    Amat -- matrix A appearing the constraint inquality

    bvec -- the vector b of constants in the constraint inequality

    meq  -- the first meq (integer) constraints of the inequality constraint will be treated as equalities

    factorized -- when factorized is true, Dmat=R^{-1} where D = R^TR.

    destroyOK -- if True, Amat, dvec, Aind, bvec will be destoryed on exit but computation
    may be performed more efficiently.'''    
        
    n = Dmat.shape[0]
    q = Amat.shape[1]
    if bvec is None:
        bvec = np.zeros(q)
    if n != Dmat.shape[1]:
        raise ValueError("Dmat is not symmetric.")
    if n != len(dvec):
        raise ValueError("Dmat is not compatible with dvec.")
    if n != Amat.shape[0]:
        raise ValueError("Amat and dvec are not compatible.")
    if q != len(bvec):
        raise ValueError("Amat and bvec are incompatible.")
    if meq>q or meq<0:
        raise ValueError("meq is invalid.")  
    err = factorized
    
    sol = Solution(n,q)

    if destroyOK:
        qpgen2_inplace(
            Dmat,
            dvec,
            Amat,
            bvec,
            meq,
            err
            )
    else:
        sol.solution, sol.crval, sol.iact, sol.nact, sol.iter, sol.work = qpgen2( 
            Dmat,
            dvec,
            Amat,
            bvec,
            meq,
            err
            )
    if err == 1:
        ValueError("constraints are inconsistent, no solution.")
    if err == 2:
        ValueError("matrix D in quadratic function is not positive definite!")
    return(sol)

if __name__ == "__main__":
    pass