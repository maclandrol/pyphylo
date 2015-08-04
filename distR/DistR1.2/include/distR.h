/* Copyright November 2005: Rachel Bevan
   Email: rachel@mcb.mcgill.ca
   NOTE: NOT Done commenting
*/

#ifndef calcdistances_h
#define calcdistances_h

#include<stdio.h>
#include<stddef.h>
#include<stdlib.h>
#include<math.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include "tnt.h"
#include "phylogeny.h"
#include "wrapper.h"
#include "jama_qr.h"
#include"global.h"
#include"distanceList.h"
//using namespace TNT;
using namespace JAMA;

double LeastSquaresError(const Array2D<double> &X, const Array2D<double> &truedistances, distanceList *mylist,char **alltaxa, int numtaxa,int numtrees);
int compareSpeciesLists(char **list1, int length1, char**list2, int length2);
double ** calculateVarianceMatrixfromAlignment(double **currdistances, char **myalign, char ** alignspecies, char **distspecies, int numspecies);
double ** calculateDistanceMatrixfromTree(Treenode *mytree, int numspecies);
void calculateTreeRates(distanceList *distmats, int numtrees, char **alltaxa, int numtaxa, char *dirName);
int CalculateSequenceIndex( char *taxaname, char **alignmentnames,int numspecies);
double bulmerVarianceCalculation(char * seq1, char *seq2,double distance);
int CorrectID(char * taxaname, char **alltaxa, int numtaxa);
void createDistanceMatrix(const Array2D<double> &solutions,const Array2D<int> &indexconversions,int numtaxa,const Array1D<int> &removedflags, Array2D<double> &distances );
#endif
