#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Test routines for the quadprog package.  Excepted where noted, all examples are drawn 
from the R quadprog package documentation and test suite.
"""
import numpy as np
import quadprog
import unittest

class TestQuadprog(unittest.TestCase):
    
    def setUp(self):
        print(dir(quadprog))
        pass
    
    def test_solveQP_basic(self):       
        '''Solve a simple quadratic program.
         
         Minimize in x: -(0 5 0) %*% x + 1/2 x^T x
         Subject to:      A^T b >= b0
        
         with b0 = (-8,2,0)^T
        
         and      (-4  2  0) 
              A = (-3  1 -2)
                  ( 0  0  1)
        
        '''
        expected = [0.4761905, 1.0476190, 2.0952381]
        Dmat = np.identity(3)
        dvec = np.array([0,5,0])
        Amat = np.array([[-4, 2, 0],[-3, 1, -2], [0, 0, 1]])
        bvec = np.array([-8,2,0])
        sol = quadprog.solveQP(Dmat,dvec,Amat,bvec)
        print(self.test_solveQP_basic.__doc__ + '\nExpected: ' + expected.__str__())
        np.testing.assert_almost_equal(sol.solution, np.array(expected))

    def test_solveCompactFormQP_basic(self):       
        '''Solve a simple quadratic progam using the compact storage format for the constraint data.

         Minimize in x: -(0 5 0) %*% x + 1/2 x^T x
         Subject to:      A^T b >= b0
        
         with b0 = (-8,2,0)^T
        
         and      (-4  2  0) 
              A = (-3  1 -2)
                  ( 0  0  1)
         using a compact form of A.
        
         
        '''
        expected = [0.4761905, 1.0476190, 2.0952381]
        Dmat = np.identity(3)
        dvec = np.array([0,5,0])
        Aind = np.array([[2,2,2], [1,1,2], [2,2,3]])
        Amat = np.array([[-4,2,-2],[-3,1,1]])
        bvec = np.array([-8,2,0])
        sol = quadprog.solveCompactFormQP(Dmat, dvec, Amat, Aind, bvec)
        print(self.test_solveCompactFormQP_basic.__doc__+ '\nExpected: ' + expected.__str__())
        np.testing.assert_almost_equal(sol.solution, np.array(expected))

if __name__ == "__main__":
    unittest.main()

