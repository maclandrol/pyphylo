#ifndef matrixoperations_c
#define matrixoperations_c

/* Copyright June 2004: Rachel Bevan
   Email: rachel@mcb.mcgill.ca
   NOTE: NOT Done commenting
*/
/* MATRIX STUFF */
#include"matrixOperations.h"


void matrixMultiplyDiagA(Array2D<double> &result,const Array2D<double> &diag, const Array2D<double> &A){

  Array2D<double> intermediate(A.dim1(), A.dim2(),0.0);
  int i,j;
  int dim1=A.dim1(), dim2=A.dim2();
  
  for(i=0;i<dim1;i++){
    for(j=0;j<dim2;j++){
      intermediate[i][j]=A[i][j]*diag[i][0];
    }
  }
  result=intermediate;
  return;
}


void matrixMultiplyBTransDiag(Array2D<double> &result,const Array2D<double> &BTrans, const Array2D<double> &diag1){

  Array2D<double> intermediate(BTrans.dim1(), BTrans.dim2(),0.0);
  int i,j;
  int dim1=BTrans.dim1(), dim2=BTrans.dim2();
  
  for(i=0;i<dim1;i++){
    for(j=0;j<dim2;j++){
      intermediate[i][j]=BTrans[i][j]*diag1[j][0];
    }
  }
  result=intermediate;
  return;
}

void matrixMultiply(Array2D<double> & result, const Array2D<double> &matrix1, const Array2D<double> &matrix2){
  int dim1=matrix1.dim1(), dim2=matrix1.dim2(), dim3=matrix2.dim2();
  int i,j,k;
  double val1, val2;
  
  
  if(dim2==matrix2.dim1()){
    Array2D<double> result1(dim1,dim3,0.0);
    for(i=0;i<dim1;i++){
      for(j=0;j<dim2;j++){
	for(k=0;k<dim3;k++){
	  val1=matrix1[i][j];
	  val2=matrix2[j][k];
	  if((val1!=0.0)&&(val2!=0.0))
	    result1[i][k]=result1[i][k]+val1*val2;
	}
      }
    }
    result=result1;
  }
  else
    fprintf(stderr, "Cannot multiply matrices! Dimensions are wrong!\n");

  return;

}



void matrixDifference(Array2D<double> &result, const Array2D<double> &matrix1, const Array2D<double> &matrix2){
  int dim1=matrix1.dim1(), dim2=matrix1.dim2();
  int i,j;
    
  if((dim1==matrix2.dim1())&&(dim2==matrix2.dim2())){
    Array2D<double> result1(dim1,dim2,0.0);
    for(i=0;i<dim1;i++){
      for(j=0;j<dim2;j++){
	result1[i][j]=matrix1[i][j]-matrix2[i][j];
      }
    }
    result=result1;
  }
  else
    fprintf(stderr, "Cannot subtract matrices! Different dimensions!\n");

}

void matrixAddition(Array2D<double> &result, const Array2D<double> &matrix1, const Array2D<double> &matrix2){
  int dim1=matrix1.dim1(), dim2=matrix1.dim2();
  int i,j;
    
  if((dim1==matrix2.dim1())&&(dim2==matrix2.dim2())){
    Array2D<double> result1(dim1,dim2,0.0);
    for(i=0;i<dim1;i++){
      for(j=0;j<dim2;j++){
	result1[i][j]=matrix1[i][j]+matrix2[i][j];
      }
    }
    result=result1;
  }
  else
    fprintf(stderr, "Cannot subtract matrices! Different dimensions!\n");

}

void vectorMultiply(double *result, const Array2D<double> &vec1, const Array2D<double> &vec2){

  int dim1=vec1.dim2(),dim2=vec2.dim1();
  int i;
  (*result)=0.0;
  //fprintf(stderr, "dim1 %i dim2 %i dim1 %i dim2 %i \n", vec1.dim1(), vec1.dim2(),vec2.dim1(),vec2.dim2());
  if((dim1==dim2)&&(vec1.dim1()==1)){
    for(i=0;i<dim1;i++){
      (*result)=(*result)+vec1[0][i]*vec2[i][0];
      
    }
  }
}


#endif
