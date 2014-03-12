#include <iostream>
#include <stdio.h>
#include "Matrix.h"
#include "Vector.h"
#include "Solver.h"
#include "IterativeSolver.h"
#include "JacobiOR.h"
#include "GaussSeidelOR.h"
#include "SteepestDescent.h"
#include "PrecSteepestDescent.h"
#include "ConjugGrad.h"
#include "PrecConjugGrad.h"

int main()
{
  //using namespace std;
  
 
  
  Matrix M;
  FILE * f = fopen("posdefmat32.dat", "rt");
  int size = M.getSizeOfFile(f);
  fclose(f);
  f = fopen("posdefmat32.dat", "rt");
  M.setMatrix(f, size);
  M.printMatrix();
  
  
  Vector b;
  FILE * f2 = fopen("vec_posdefmat32.dat", "rt");
  int size2 = b.getSizeOfFile(f2);
  fclose(f2);
  f2 = fopen("vec_posdefmat32.dat", "rt");
  b.setVector(f2, size2);
  b.printVector();
  
  /*
 Matrix M;
  FILE * f = fopen("mat64.dat", "rt");
  int size = M.getSizeOfFile(f);
  fclose(f);
  f = fopen("mat64.dat", "rt");
  M.setMatrix(f, size);
  M.printMatrix();
  
  
  Vector b;
  FILE * f2 = fopen("vec64.dat", "rt");
  int size2 = b.getSizeOfFile(f2);
  fclose(f2);
  f2 = fopen("vec64.dat", "rt");
  b.setVector(f2, size2);
  b.printVector();
  */
  /*
  Vector c;
  FILE * f3 = fopen("vec2.dat", "rt");
  int size3 = c.getSizeOfFile(f3);
  fclose(f3);
  f3 = fopen("vec2.dat", "rt");
  c.setVector(f3, size3);
  c.printVector();
 
  Vector V;
  V.resize(5);
  
  
  Vector d;
  //double sum = dotProduct(c,b);
  d = b+c;
  d.printVector();

  Vector r1, r2;

  r1=d*M;
  r2 = M*d;
  r1.printVector();
  r2.printVector();
  
  Matrix A;
  FILE * f4 = fopen("mat2.dat", "rt");
  int size4 = A.getSizeOfFile(f4);
  fclose(f4);
  f4 = fopen("mat2.dat", "rt");
  A.setMatrix(f4, size4);
  A.printMatrix();
  
  r2 = (A+M)*d;
  r2.printVector();
  Matrix K;
  K = A;
  K.printMatrix();
  Matrix H;
  //H = M*K;
  //H.printMatrix();
  K = M+A;
  M.printMatrix();
  K.printMatrix();
  H = M*K;
  H.printMatrix();
  //K[2][2] = 100;

  M.printMatrix();
  d.printVector();
  //K(2,2)=555;
  K.printMatrix();
  */
   
  /*
  Solver s(M, b);
  //  Vector solutionSimple = s.solveBySimpleIterations();
  Vector solutionJOR(3);
  s.solveByJacobiOR(solutionJOR, x0, 100, 0.0001, 1);
  solutionJOR.printVector();
  
  x0.setToZero();
  Vector solutionSOR(3);
  s.solveByGaussSeidelOR(solutionSOR, x0, 100, 0.0001, 1);
  solutionSOR.printVector();
  
  x0.setToZero();
  Vector solutionSteepestDescent(3);
  s.solveBySteepestDescent(solutionSteepestDescent, x0,100, 0.0001);
  solutionSteepestDescent.printVector();

  x0.setToZero();
  Vector solutionConjGrad(3);
  s.solveByConjGrad(solutionConjGrad, x0,100, 0.0001);
  solutionConjGrad.printVector();
  */
  int eqn = M.giveNumberOfRows();
  double tol = 0.00001;
  int iter = 1000;
  Vector x0(eqn);
  x0.setToZero();
  /*
  JacobiOR s(M,b, tol, iter);
  Vector solutionJOR(eqn);
  s.solve(solutionJOR, x0, 1);
  solutionJOR.printVector();
  
  GaussSeidelOR s2(M,b, tol, 100000);
  Vector solutionSOR(eqn);
  s2.solve(solutionSOR, x0, 1);
  solutionSOR.printVector();

  SteepestDescent s3(M,b, tol, iter);
  Vector steepest(eqn);
  s3.solve(steepest, x0);
  steepest.printVector();
  */
  ConjugGrad s4(M,b, tol, iter);
  Vector conjug(eqn);
  s4.solve(conjug, x0);
  conjug.printVector();
  /*
  PrecSteepestDescent s5(M,b, tol, iter);
  Vector Precsteepest(eqn);
  s5.solve(Precsteepest, x0, 1);
  Precsteepest.printVector();
  */
  PrecConjugGrad s6(M,b, tol, iter);
  Vector Precconjug(eqn);
  s6.solve(Precconjug, x0, 1);
  Precconjug.printVector();
  int konec = 2;

  Matrix K;
  int neq = 100;
  K.beRandomMatrix(3);
  K.printMatrix();
  K.beSymmetricPart();
  K.printMatrix();
  Matrix L;
  L.bePositiveDefiniteMatrix(neq);
  //L.printMatrix();
  Vector c;
  c.beTestRhsVector(L);
  x0.resize(neq);
  x0.setToZero();
  //x0.setValues(100,1);
  //x0.setValues(1000,2);
  //x0.setValues(10000,3);
  //x0.setValues(100000,4);
  //PrecConjugGrad s7(L,c, tol, iter);
  ConjugGrad s7(L,c, tol, iter);
  Vector Con(neq);
  
  s7.solve(Con,x0);
  Con.printVector();
  PrecConjugGrad s8(L,c, tol, iter);
  Vector PreCon(neq);
  s8.solve(PreCon,x0,1);
  PreCon.printVector();
  return 0;



}
