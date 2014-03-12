#include <iostream>
#include <stdio.h>
#include "Matrix.h"
#include "Vector.h"
#include "Solver.h"
#include <cstdlib>
#include <ctime>
#include <cmath>


  


void Solver :: printGeneralSolverHeader()
{
  printf("___________________________________________________________\n");
  printf("               MATCOM - Linear system solver               \n");
  printf("              Copyright (C) 2012 Filip Kolarik             \n");
  printf("___________________________________________________________\n");
  printf("                                                           \n");
  printf("Number of equations: %d\n",this->A.giveNumberOfRows() );
  //printf("Iteration" );
}




Vector Solver :: solveBySimpleIterations()
{
  double norm = 1.0;
  double norm0, s;
  Matrix I;
  int size;
  size = this->A.giveNumberOfRows();
  I.resize(size, size );
  I.beUnitMatrix();
  Vector x(size); 
  Vector x_new(size); 
  Vector r(size);
  Vector r0(size);
  for(int i = 1; i <= size; i++){
    for(int j = 1; j <= size; j++){
      r0.at(i) = b.at(i);
    }
  }
  norm = r0.computeEuclideanNorm();
  norm0 = norm;

  while (norm > 0.0001){
   
    for(int i = 1; i<=size;i++){
      for(int j = 1; j<=size;j++){ 
	x_new.at(i) = ( A.at(i,j) + I.at(i,j) ) * x.at(j) - b.at(i);
      }
      x.at(i) = x_new.at(i);
    }
    
    for(int i = 1; i<=size;i++){
      s = 0;
      for(int j = 1; j <= size; j++){
	s += A.at(i,j) * x_new.at(j);
      }
      r.at(i) = b.at(i) - s; 
    }
    norm = r.computeEuclideanNorm() / norm0; 
  
  }
  return x;
}
 

void Solver :: solveByJacobiOR(Vector &x, Vector x0, int maxiter, double tol, int omega)
{
  
  double norm = 1.0;
  double s, norm0;
  int size, iter = 0;
  size = this->A.giveNumberOfRows();
  Vector r0, x_new, r;
  x_new.resize(size);
  r.resize(size);
  r0.resize(size);
  //x.resize(size);
  for(int i = 1; i <= size; i++){
    s = 0;
    for(int j = 1; j <= size; j++){
      s += A.at(i,j) * x0.at(j);
    }
  r0.at(i) = b.at(i) - s;
  }
  
  x.setToZero();
  norm = r0.computeEuclideanNorm();

  norm0 = norm;
  while ((norm > tol) && (iter < maxiter)){
    iter++;
    
    for(int i = 1; i <= size; i++){
      s=0;
      for(int j = 1; j <= i-1; j++){
	s += A.at(i,j) * x.at(j);
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
    
    norm = r.computeEuclideanNorm() / norm0; 
  }
 
}



void Solver :: solveByGaussSeidelOR(Vector &x, Vector x0, int maxiter, double tol, int omega)
{
   
  double norm = 1.0;
  double s, norm0;
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
  norm0 = norm;
  while ((norm > tol) && (iter < maxiter)){
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

    norm = r.computeEuclideanNorm() / norm0;
  }
 
}

void Solver :: solveBySteepestDescent(Vector &x, Vector x0, int maxiter, double tol)
{
 double norm = 1.0;
 double s, norm0, rho, alpha;
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
  norm0 = norm;
  
  while ((norm > tol) && (iter < maxiter)){
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
    norm = r.computeEuclideanNorm() / norm0; ;
  }

}

void Solver :: solveByConjGrad(Vector &x, Vector x0, int maxiter, double tol)
{
  /*function [x,relres,iter,flag]=conjgrad(A,b,x,P,nmax,tol)
    %CONJGRAD Conjugate gradient method
    %  [X,RELRES,ITER,FLAG]=CONJGRAD(A,B,X0,NMAX,TOL,OMEGA) attempts 
    %  to solve the system A*X=B with the conjugate gradient method. TOL 
    %  specifies the tolerance of the method. NMAX specifies the maximum number of 
    %  iterations. X0 specifies the initial guess. P is a preconditioner. RELRES is the 
    %  relative residual. If FLAG is 1 RELRES > TOL. ITER is the 
    %  iteration number at which X is computed.
    flag=0; iter=0; bnrm2=norm(b);
    if bnrm2==0, bnrm2=1; end
    r=b-A*x; relres=norm(r)/bnrm2;
    if relres<tol, return, end
    for iter = 1:nmax 
    z=P\r; rho=r'*z;
    if iter>1
    beta=rho/rho1;
    p=z+beta*p;
    else
    p=z;
    end
    q=A*p;              
    alpha=rho/(p'*q);
    x=x+alpha*p;
    r=r-alpha*q;
    relres=norm(r)/bnrm2;
    if relres<=tol, break, end
    rho1 = rho;
    end
    if relres>tol, flag = 1; end
    return
  */

 double norm = 1.0;
 double s, norm0, rho, alpha, beta;
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
  norm0 = norm;
  
  while ((norm > tol) && (iter < maxiter)){
    iter++; 
    double* res = r.values();
    rho = dotProduct(res,res,size);
    if (iter == 1){
      for(int i = 1; i <= size; i++){
	p.at(i) = r.at(i);
      }
    } else {
      for(int i = 1; i <= size; i++){
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
    norm = r.computeEuclideanNorm() / norm0; ;
  }


}





Vector Solver :: solveByPrecondSteepestDescent(Vector x0, int maxiter, double tol)
{
  /*[n,m]=size(A);
    if n ~= m, error('Only squared systems'); end
    flag = 0; iter = 0; bnrm2 = norm( b );
    if bnrm2==0, bnrm2 = 1; end
    r=b-A*x;  relres=norm(r)/bnrm2;
    if relres<tol, return, end
    for iter=1:nmax
    z=P\r;  
    rho=r'*z;
    q=A*z;
    alpha=rho/(z'*q);
    x=x+alpha*z;
    r=r-alpha*q;
    relres=norm(r)/bnrm2;
    if relres<=tol, break, end
    end
    if relres>tol, flag = 1; end
    return
  */



}


Vector Solver :: solveByPrecondConjGrad(Vector x0, int maxiter, double tol)
{
  /*function [x,relres,iter,flag]=conjgrad(A,b,x,P,nmax,tol)
    %CONJGRAD Conjugate gradient method
    %  [X,RELRES,ITER,FLAG]=CONJGRAD(A,B,X0,NMAX,TOL,OMEGA) attempts 
    %  to solve the system A*X=B with the conjugate gradient method. TOL 
    %  specifies the tolerance of the method. NMAX specifies the maximum number of 
    %  iterations. X0 specifies the initial guess. P is a preconditioner. RELRES is the 
    %  relative residual. If FLAG is 1 RELRES > TOL. ITER is the 
    %  iteration number at which X is computed.
    flag=0; iter=0; bnrm2=norm(b);
    if bnrm2==0, bnrm2=1; end
    r=b-A*x; relres=norm(r)/bnrm2;
    if relres<tol, return, end
    for iter = 1:nmax 
    z=P\r; rho=r'*z;
    if iter>1
    beta=rho/rho1;
    p=z+beta*p;
    else
    p=z;
    end
    q=A*p;              
    alpha=rho/(p'*q);
    x=x+alpha*p;
    r=r-alpha*q;
    relres=norm(r)/bnrm2;
    if relres<=tol, break, end
    rho1 = rho;
    end
    if relres>tol, flag = 1; end
    return
  */


}
