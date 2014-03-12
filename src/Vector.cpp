
#include <iostream>
#include <stdio.h>
#include "Vector.h"
#include "Matrix.h"
#include <cstdlib>
#include <ctime>
#include <cmath>

void Vector :: resize(int rows)
{


  this->nRows = rows;
 
  this->val = new double [nRows];

  //return this;
  //return;
}


void Vector :: setValues(double h, int i)
{
  this->at(i) = h;
}


void Vector :: setToZero()
{
  for(int j = 1; j <= nRows; j++){
    this->at(j) = 0;   
  }
  
}

void Vector :: setToOnes()
{
  for(int j = 1; j <= nRows; j++){
    this->at(j) = 1;   
  }
  
}




Vector Vector :: operator+(const Vector &b) const
{
  
  Vector result;
  result.resize(this->nRows);

  for ( int i = 1; i<=nRows; i++){

    result.at(i) = this->at(i) + b.at(i);
  }

  return result;
} 

Vector Vector :: operator-(Vector &b)
{
  
  Vector result;
  result.resize(this->nRows);

  for ( int i = 1; i<=nRows; i++){

    result.at(i) = this->at(i) - b.at(i);
  }

  return result;
} 


Vector Vector :: operator*(Matrix &m)
{
  
  Vector result;
  int size = m.giveNumberOfRows();
  result.resize(size);

  double sum=0.0;
  for ( int i = 1; i<=size; i++){
    for ( int j = 1; j<=this->nRows; j++){
      
      sum += m.at(i,j) * this->at(j);
    }
    result.at(i) = sum;
    sum = 0.0;
    
    //return *this;
  }
  return result;
  
}



void Vector :: printVector()
{
  int i;
  int nRows = this->giveNumberOfRows();
  
  for ( i = 1; i <= nRows; i++ ) {
           
    printf( "%10.3e  ", this->at(i));
       

    printf("\n");
  }

}
int Vector :: getSizeOfFile(FILE *f)
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

void Vector :: setVector(FILE *f, int size)
{
 
  int i = 1;
  int check;
 
  this->resize(size);
  
  while (check != EOF){
    
    check = fscanf(f, "%lf", &this->at(i));
    //check = fscanf(f, "%lf", &pokus);
    i++;
  }
  
  
  fclose(f);
}

double dotProduct(double *a, double *b, int i)
{
  double sum = 0.0;
  while (i--){
    sum += (*a++) * (*b++);
  }  
return sum;

}


double Vector :: computeEuclideanNorm()
{
  double sum = 0.0;
  
    sum = dotProduct(val,val, nRows);
  

  return sqrt(sum);
}


void Vector :: beTestRhsVector(Matrix &A)
{
  int n = A.giveNumberOfRows();
  Vector tmp;
  tmp.resize(n);
  tmp.setToOnes();
  this->resize(n);
  for ( int i = 1; i<=n; i++){
    for ( int j = 1; j<=n; j++){
      this->at(i) += A.at(i,j)*tmp.at(j);
    }
  }

    
}
