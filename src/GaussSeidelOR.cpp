#include <iostream>
#include <stdio.h>
#include "Matrix.h"
#include "Vector.h"
#include "Solver.h"
#include "IterativeSolver.h"
#include "GaussSeidelOR.h"
#include <cstdlib>
#include <ctime>
#include <cmath>


  



void GaussSeidelOR :: solve(Vector &x, Vector x0, int omega)
{
  this->printGeneralSolverHeader();
  this->printIterativeSolverHeader();
  
  double norm = 1.0;
  double s, norm0, absnorm;
  int size, iter = 0;
  size = this->A.giveNumberOfRows();
  Vector r0(size), x_new(size), r(size);
  for(int i = 1; i <= size; i++){
    s = 0;
    for(int j = 1; j <= size; j++){
      s += A.at(i,j) * x0.at(j);
    }
    r0.at(i) = b.at(i) - s;
    x.at(i) = x0.at(i);
  }
  
  x_new.setToZero();
  norm = r0.computeEuclideanNorm();
  absnorm = norm;
  norm0 = norm;
  while (((norm > tol) || (absnorm > tol)) && (iter < maxiter)){
    iter++;
    
    for(int i = 1; i <= size; i++){
      s=0;
      for(int j = 1; j <= i-1; j++){
	s += A.at(i,j) * x_new.at(j);
      }
      for(int j = i+1; j <= size; j++){
	s += A.at(i,j) * x.at(j);
      }
      
      x_new.at(i) = omega * (b.at(i) - s) / A.at(i,i) + (1 - omega) * x.at(i);
    }
    
    for(int i = 1; i <= size; i++){
      s=0;
      for(int j = 1; j <= size; j++){
	s += A.at(i,j) * x_new.at(j);
      }
      r.at(i) = b.at(i) - s;
      x.at(i) = x_new.at(i);	
    }
    absnorm = r.computeEuclideanNorm();
    norm = r.computeEuclideanNorm() / norm0;
    this->printIterativeSolverInfo(iter, absnorm, norm);
  }
 
}
