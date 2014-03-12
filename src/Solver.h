class Vector;
class Matrix;

class Solver
{
 protected:
  
  Matrix A;
  Vector b;

 public:
  Solver(Matrix &mat, Vector &vec)
    {
      A = mat;
      b = vec;
    }
  
  
  Solver(){}

  ~Solver(){}
  

  void printGeneralSolverHeader();

  
  Vector solveBySimpleIterations();
  void solveByJacobiOR(Vector &solution, Vector x0, int maxiter, double tol, int omega);
  void solveByGaussSeidelOR(Vector &solution, Vector x0, int maxiter, double tol, int omega);
  void solveBySteepestDescent(Vector &solution, Vector x0, int maxiter, double tol);
  void solveByConjGrad(Vector &solution, Vector x0, int maxiter, double tol);
  Vector solveByPrecondSteepestDescent(Vector x0, int maxiter, double tol);
  Vector solveByPrecondConjGrad(Vector x0, int maxiter, double tol);


};
