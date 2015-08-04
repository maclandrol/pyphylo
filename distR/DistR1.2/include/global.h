/*******************************************************
filename:       global.H
author: 	David Bryant (bryant@math.mcgill.ca)
created: 	12 Aug 1999
copyright:	David Bryant 1999. 

This contains all the external include files for the
project, all of the global macros and constants, and
the routine error_message(" ") which prints error message and aborts.


********************************************************/

#ifndef GLOBAL_H
#define GLOBAL_H

//Standard and STL include files
#include<iostream>
#include<stdio.h>
#include<string>
#include<unistd.h>
#include<vector>
#include<math.h>
#include<fstream>
#include<cstdlib>

using namespace std;


//Constants
#define BIG_NUMBER 1000000 //Used in all maximisation loops
#define BIG_FLOAT 1.0e300
#define INF 1.0e100
//Macros
#define ABS(g) ((g>0) ? (g) : -(g) )
#define MAX(A,B) ((A) > (B) ? (A):(B))
#define MIN(A,B) ((A) < (B) ? (A):(B))

typedef vector<vector<double> > dbl_array; /* STL implementation of a matrix */

//Utility routines
void error_message(const string& message); /* prints error message and aborts. */

#endif
