class Vector;
class Matrix;

//#include "IterativeSolver.h"


class JacobiOR : public IterativeSolver
{
  //protected:
  

 public:
  JacobiOR(Matrix &mat, Vector &vec, double tolerance, int maxit) : 
  IterativeSolver(mat, vec,tolerance, maxit)
    {
      
    }
  
  


  ~JacobiOR(){}
  void solve(Vector &x, Vector x0, int omega);
  virtual string getTypeOfSolver(){return "Jacobi over relaxation method";}

};
