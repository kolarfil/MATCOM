class Vector;
class Matrix;

//#include "IterativeSolver.h"


class PrecConjugGrad : public IterativeSolver
{
  //protected:
  

 public:
  PrecConjugGrad(Matrix &mat, Vector &vec, double tolerance, int maxit) : 
  IterativeSolver(mat, vec,tolerance, maxit)
    {
      
    }
  
  


  ~PrecConjugGrad(){}
  void solve(Vector &x, Vector x0, int precType);
  virtual string getTypeOfSolver(){return "Preconditioned Conjugate Gradient method";}

};
