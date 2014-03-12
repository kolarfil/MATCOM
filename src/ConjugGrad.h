class Vector;
class Matrix;

//#include "IterativeSolver.h"


class ConjugGrad : public IterativeSolver
{
  //protected:
  

 public:
  ConjugGrad(Matrix &mat, Vector &vec, double tolerance, int maxit) : 
  IterativeSolver(mat, vec,tolerance, maxit)
    {
      
    }
  
  


  ~ConjugGrad(){}
  void solve(Vector &x, Vector x0);
  virtual string getTypeOfSolver(){return "Conjugate Gradient method";}

};
