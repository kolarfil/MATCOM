class Vector;
class Matrix;

//#include "IterativeSolver.h"


class PrecSteepestDescent : public IterativeSolver
{
  //protected:
  

 public:
 PrecSteepestDescent(Matrix &mat, Vector &vec, double tolerance, int maxit) : 
  IterativeSolver(mat, vec,tolerance, maxit)
    {
      
    }
  
  


  ~PrecSteepestDescent(){}
  void solve(Vector &x, Vector x0, int precType);

  virtual string getTypeOfSolver(){return "Preconditioned Steepest Descent method";}
};
