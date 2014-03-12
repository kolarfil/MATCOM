class Vector;


class Matrix
{
 protected:
  
  int nRows;
  int nColumns;
  double **A;

 public:
  Matrix(int n, int m)
    {
      nRows = n;
      nColumns = m;
      
      A = new double *[nRows];
      for(int i = 0; i< nRows; i++){
	A[i] = new double [nColumns];
      }
    }
  
  Matrix(){}
    

  ~Matrix()
    {
      /* if(A){ */
      /* 	for (int i;i<nRows;i++){ */
      /* 	  delete []A[i]; */
      /* 	} */
      /* 	delete []A; */
      /* } */
    }
  
  int giveNumberOfRows() { return nRows; }
  
  int giveNumberOfColumns() { return nColumns; }
  void resize(int rows, int columns);
  //Matrix* resize(int rows, int columns);
  double at( int i, int j) const { return A[i-1][j-1]; }
  double& at( int i, int j) { return A[i-1][j-1]; }
 
  
  void printMatrix();
  void beRandomMatrix(int n);
  void beSymmetricPart();
  void bePositiveDefiniteMatrix(int n);
  
  void setMatrix(FILE *f, int size);
  int getSizeOfFile(FILE *f);
  void beUnitMatrix();

  Matrix operator+(const Matrix &m) const;
  Vector operator*(Vector &v);
  Matrix operator*(Matrix &m);
};
