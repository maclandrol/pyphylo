/* Copyright November 2005: Rachel Bevan
   Email: rachel@mcb.mcgill.ca
   NOTE: NOT Done commenting
*/

#ifndef matrixoperations_h
#define matrixoperations_h

#include <stdlib.h>
#include <stddef.h>
#include "tnt.h"
#include "jama_qr.h"
#include"global.h"
//using namespace TNT;
using namespace JAMA;

void matrixMultiplyBTransDiag(Array2D<double> &result, const Array2D<double> &BTrans, const Array2D<double> &diag1);
void matrixMultiplyDiagA(Array2D<double> &result,const Array2D<double> &diag, const Array2D<double> &A);
void matrixMultiply(Array2D<double> & result,const Array2D<double> &matrix1, const Array2D<double> &matrix2);
void matrixDifference(Array2D<double> &result, const Array2D<double> &matrix1, const Array2D<double> &matrix2);
void vectorMultiply(double *result, const Array2D<double> &vec1, const Array2D<double> &vec2);
void matrixAddition(Array2D<double> &result, const Array2D<double> &matrix1, const Array2D<double> &matrix2);
#endif
