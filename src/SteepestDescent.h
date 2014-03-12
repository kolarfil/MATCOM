class Vector;
class Matrix;

//#include "IterativeSolver.h"


class SteepestDescent : public IterativeSolver
{
  //protected:
  

 public:
  SteepestDescent(Matrix &mat, Vector &vec, double tolerance, int maxit) : 
  IterativeSolver(mat, vec,tolerance, maxit)
    {
      
    }
  
  


  ~SteepestDescent(){}
  void solve(Vector &x, Vector x0);

  virtual string getTypeOfSolver(){return "Steepest Descent method";}
};
