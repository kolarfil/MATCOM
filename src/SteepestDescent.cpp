#include <iostream>
#include <stdio.h>
#include "Matrix.h"
#include "Vector.h"
#include "Solver.h"
#include "IterativeSolver.h"
#include "SteepestDescent.h"
#include <cstdlib>
#include <ctime>
#include <cmath>


  



void SteepestDescent :: solve(Vector &x, Vector x0)
{

  this->printGeneralSolverHeader();
  this->printIterativeSolverHeader();
  double norm = 1.0;
  double s, norm0, absnorm, rho, alpha;
  int size, iter = 0;
  size = this->A.giveNumberOfRows();
  x.setToZero();
  Vector r(size);
  Vector q(size);
 
  for(int i = 1; i <= size; i++){
    s = 0;
    for(int j = 1; j <= size; j++){
      s += A.at(i,j) * x0.at(j);
    }
    r.at(i) = b.at(i) - s;
    x.at(i) = x0.at(i);
  }
  norm = r.computeEuclideanNorm();
  absnorm = norm;
  norm0 = norm;
  
  while (((norm > tol) || (absnorm > tol)) && (iter < maxiter)){
    iter++; 
    double* res = r.values();
    rho = dotProduct(res,res,size);
    for(int i = 1; i <= size; i++){
      s = 0;
      for(int j = 1; j <= size; j++){
	s += A.at(i,j) * r.at(j);
      }
      q.at(i) = s;
    }
    alpha = rho / dotProduct(res,q.values(), size);
    for(int i = 1; i <= size; i++){
      x.at(i) += alpha * r.at(i);
    }
    for(int i = 1; i <= size; i++){
      r.at(i) -= alpha * q.at(i);   
    }
    absnorm = r.computeEuclideanNorm();
    norm = r.computeEuclideanNorm() / norm0; 
    this->printIterativeSolverInfo(iter, absnorm, norm);
  }

}
