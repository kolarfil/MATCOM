#include <iostream>
#include <stdio.h>
#include "Matrix.h"
#include "Vector.h"
#include "Solver.h"
#include "IterativeSolver.h"
#include <cstdlib>
#include <ctime>
#include <cmath>



  

void IterativeSolver :: printIterativeSolverHeader()
{
  string type;
  type = this->getTypeOfSolver();
  printf("Type of solver: ");
  std::cout << type << std::endl;
  printf("Iteration             Error           Relative Error       \n");

  printf("___________________________________________________________\n");
  
}

void IterativeSolver :: printIterativeSolverInfo(int iter, double Error, double relError)
{

  
  printf("%d                 %10.5e           %10.5e                 \n", iter, Error, relError);


}

void IterativeSolver :: getDiagonalPreconditioningMatrix()
{
  
  int size = this->A.giveNumberOfRows();
  Matrix* Preconditioner = new Matrix;
  Preconditioner->resize(size, size);
  P = Preconditioner;
  for ( int i = 1; i<=size; i++){
    for ( int j = 1; j<=size; j++){ 
      if ( i == j){
	P->at(i,j) = A.at(i,j);
      } else{
	P->at(i,j) = 0;
      }
    }
  }
}

void IterativeSolver :: computePreconditionedResiduum(Vector &preres, Vector &res, int precType)
{

  if(precType == 1){
    this->computeDiagonalPreconditionedResiduum(preres, res);
  } else{
    printf("unknown preconditioner type");
    //break;
  }
}


void IterativeSolver :: computeDiagonalPreconditionedResiduum(Vector &z, Vector &r)
{
  
  if(P){
    for ( int i = 1; i<=this->P->giveNumberOfRows(); i++){
      z.at(i) = r.at(i) / P->at(i,i);
    }
  } else{
    this->getDiagonalPreconditioningMatrix();
      for ( int i = 1; i<=this->P->giveNumberOfRows(); i++){
	z.at(i) = r.at(i) / P->at(i,i);
      }
    
  }
}


