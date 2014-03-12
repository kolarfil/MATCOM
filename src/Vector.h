class Matrix;



class Vector
{
 protected:
  int nRows;
  double *val;

 public:
  Vector(int n)
    {
      nRows = n;
      val = new double [nRows];
      
    }
  
  Vector(){}

  ~Vector()
    {
      //  if(val){delete []val;}
    }

  double* values(){return val;}


  void setValues(double h, int i);

  Vector operator+(const Vector &b) const;
  Vector operator-(Vector &b);
  Vector operator*(Matrix &m);
  int giveNumberOfRows(){return nRows;}

  void printVector();
  double at(int i)const {return val[ i -1 ];}
  double& at(int i) {return val[ i -1 ];}
  void resize(int rows);
  void setToZero();
  void setToOnes();
  void beTestRhsVector(Matrix &A);
  int getSizeOfFile(FILE *f);
  void setVector(FILE *f, int size);
 
  double computeEuclideanNorm();
  
  
};

double dotProduct(double *a, double *b, int i);
