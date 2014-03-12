class Vector;
class Matrix;

//#include "IterativeSolver.h"


class GaussSeidelOR : public IterativeSolver
{
  //protected:
  

 public:
  GaussSeidelOR(Matrix &mat, Vector &vec, double tolerance, int maxit) : 
  IterativeSolver(mat, vec,tolerance, maxit)
    {
      
    }
  
  


  ~GaussSeidelOR(){}
  void solve(Vector &x, Vector x0, int omega);

  virtual string getTypeOfSolver(){return "Gauss-Seidel over relaxation method";}
};
