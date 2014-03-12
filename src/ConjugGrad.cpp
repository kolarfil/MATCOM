#include <iostream>
#include <stdio.h>
#include "Matrix.h"
#include "Vector.h"
#include "Solver.h"
#include "IterativeSolver.h"
#include "ConjugGrad.h"
#include <cstdlib>
#include <ctime>
#include <cmath>


void ConjugGrad :: solve(Vector &x, Vector x0)
{
 

  this->printGeneralSolverHeader();
  this->printIterativeSolverHeader();
  double norm = 1.0;
  double s, norm0, absnorm, rho, rho1, alpha, beta;
  int size, iter = 0;
  size = this->A.giveNumberOfRows();
  x.setToZero();
  Vector r(size);
  Vector q(size);
  Vector p(size);
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
    if (iter == 1){
      for(int i = 1; i <= size; i++){
	p.at(i) = r.at(i);
      }
    } else {
      for(int i = 1; i <= size; i++){
	beta = rho/rho1;
	p.at(i) = r.at(i) + beta * p.at(i);
      }
    }
    for(int i = 1; i <= size; i++){
      s = 0;
      for(int j = 1; j <= size; j++){
	s += A.at(i,j) * p.at(j);
      }
      q.at(i) = s;
    }
    alpha = rho / dotProduct(p.values(),q.values(), size);
    for(int i = 1; i <= size; i++){
      x.at(i) += alpha * p.at(i);
    }
    for(int i = 1; i <= size; i++){
      r.at(i) -= alpha * q.at(i);   
    }
    rho1 = rho;
    absnorm = r.computeEuclideanNorm();
    norm = r.computeEuclideanNorm() / norm0; 
    this->printIterativeSolverInfo(iter, absnorm, norm);
  }


}
