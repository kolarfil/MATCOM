

class Vector;
class Matrix;
  
using namespace std;

class IterativeSolver : public Solver
{
 protected:
  double tol;
  int maxiter;
  Matrix* P; //Preconditioner

 public:
 IterativeSolver(Matrix &mat, Vector &vec, double tolerance, int maxit):
  Solver(mat, vec)
    {
      tol = tolerance;
      maxiter = maxit;
      P = NULL;
    }
  
  ~IterativeSolver()
    {
      /* if(P){
	delete P;
	}*/
    }
  
  void printIterativeSolverHeader();
  void printIterativeSolverInfo(int iter, double Error, double relError);
  virtual string getTypeOfSolver() = 0;
  
  void getDiagonalPreconditioningMatrix();
  void computePreconditionedResiduum(Vector &preres, Vector &res, int precType);
  void computeDiagonalPreconditionedResiduum(Vector &preres, Vector &res);
};
