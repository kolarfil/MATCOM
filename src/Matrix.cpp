
#include <iostream>
#include <stdio.h>
#include "Matrix.h"
#include "Vector.h"
#include <cstdlib>
#include <ctime>
#include <cmath>

/*//copy constructor
Matrix :: Matrix(const Matrix& m)
{


}*/


void Matrix :: resize(int rows, int columns)
//Matrix* Matrix :: resize(int rows, int columns)
{


  this->nRows = rows;
  this->nColumns = columns;

  this->A = new double *[nRows];
  for(int i = 0; i< nRows; i++){
    A[i] = new double [nColumns];
  }
  
}



Matrix Matrix :: operator+(const Matrix &m) const
//Matrix Matrix :: operator+(Matrix &m)
{
  
  //int size = m.giveNumberOfRows();
  Matrix result;
  result.resize( this->nRows, this->nColumns); 
  for ( int i = 1; i<=nRows; i++){
    for ( int j = 1; j<=nColumns; j++){
      result.at(i,j) = this->at(i,j) + m.at(i,j);

    }
    
    
  } 
  return result;
}



Matrix Matrix :: operator*(Matrix &m)
{
  
  Matrix result;
  result.resize( this->nRows, m.nColumns); 
  double sum=0.0;
  for ( int i = 1; i<=this->nRows; i++){
    for ( int j = 1; j<=m.nColumns; j++){
      for ( int k = 1; k <= m.nRows; k++){
	sum += this->at(i,k) * m.at(k,j);
      }      
      result.A[i-1][j-1] = sum;
      sum = 0.0;
    } 
    //return *this;
  }
  return result;
}


Vector Matrix :: operator*(Vector &v)
{
  
  Vector result;
  result.resize( this->nColumns);
  int size = v.giveNumberOfRows();
  double sum=0.0;
  for ( int i = 1; i<=this->nRows; i++){
    for ( int j = 1; j<=size; j++){
      
      sum += this->at(i,j) * v.at(j);
    }
    result.setValues(sum, i);
    sum = 0.0;
    
    //return *this;
  }
  return result;
  
}




void Matrix :: printMatrix() 
// Prints the receiver on screen.
{
    int i, j;

  int nRows = this->giveNumberOfRows();
  int nColumns = this->giveNumberOfColumns();
    
        for ( i = 1; i <= nRows; i++ ) {
            for ( j = 1; j <= nColumns; j++ ) {
                printf( "%10.3e  ", this->at(i, j) );
            }

            printf("\n");
        }
	
}


void Matrix :: beRandomMatrix(int n)
{
  srand(time(NULL));
  this->resize(n,n);
  int i, j;
  double value;
  int number, number2;
  int nRows = this->giveNumberOfRows();
  int nColumns = this->giveNumberOfColumns();
    for ( i = 1; i <= nRows; i++ ) {
      for ( j = 1; j <= nColumns; j++){
	
	number = rand();
	number2 = number % 1000;
	this->at(i,j) = number2 + pow(number, 0.01);
	      
	    }
    }
    //return A;  
}


void Matrix :: beSymmetricPart()
{

 int nRows = this->giveNumberOfRows();
  int nColumns = this->giveNumberOfColumns();
    for ( int i = 1; i <= nRows; i++ ) {
      for ( int j = i; j <= nColumns; j++){
	
	this->at(i,j) = (0.5) * ( this->at(i,j) + this->at(j,i));
	this->at(j,i) =  this->at(i,j);    
	    }
    }

}

void Matrix :: bePositiveDefiniteMatrix(int n)
{
  Matrix tmp;
  tmp.resize(n,n);
  double s=0;
  this->beRandomMatrix(n);
  this->beSymmetricPart();
    for ( int i = 1; i <= n; i++ ){
      for ( int j = i; j <= n; j++){
	s = 0;
	for ( int k = 1; k <= n; k++ ) {
	  s += this->at(i,k) * this->at(k,j);
	}
	tmp.at(i,j) = s;
	tmp.at(j,i) = s;    
      }
    }
    for ( int i = 1; i <= n; i++ ){
      for ( int j = 1; j <= n; j++){
	this->at(i,j) = tmp.at(i,j);
      }
    }
}

void Matrix :: setMatrix(FILE *f, int size)
{

 
  double * help_vector;
  //help_vector = new double [9];
  int i = 0;
  int check;
  int k = 0;
  
  help_vector = new double [size];
   
 while (check != EOF){
   
    check = fscanf(f, "%lf", &help_vector[i]);
    //check = fscanf(f, "%lf", &pokus);
    i++;
 }

 int n = sqrt(size+1);
 this->resize(n,n);
 
 for (int i = 1; i <= n; i++ ) {
   for (int j = 1; j <= n; j++ ) {
     
     //this->A[i-1][j-1] = help_vector[k];
     this->at(i,j) = help_vector[k];
     k++;
   }
 }
 fclose(f);

 delete help_vector;
}


void Matrix :: beUnitMatrix()
{
  for ( int i = 1; i<=this->nRows; i++){
   for ( int j = 1; j<=this->nColumns; j++){ 
     if ( i == j){
       this->at(i,j) = 1;
     } else{
       this->at(i,j) = 0;
     }
   }
  }
}




int Matrix :: getSizeOfFile(FILE *f)
{
  double check, pokus;
  int i = 0;
  while (check != EOF){
   
    //check = fscanf(f, "%lf", &help_vector[i]);
    check = fscanf(f, "%lf", &pokus);
   i++;
 }
  return --i;

}
