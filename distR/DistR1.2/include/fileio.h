/* Copyright November 2005: Rachel Bevan
   Email: rachel@mcb.mcgill.ca
   NOTE: NOT Done commenting
*/

#ifndef fileio_h
#define fileio_h

#include "tnt.h"
#include "jama_qr.h"
#include<stdio.h>
#include<stddef.h>
#include<stdlib.h>
#include"wrapper.h"
#include "phylogeny.h"
#include <string.h>
#include "distanceList.h"
using namespace JAMA;

#define MAXSTRSIZE 200

int readNexusFormat(char *infile, int readnexusdistances,double ***distances, int *numspecies, char ***speciesnames,char ***alignment,int readnexusalignment);
distanceList *readInputFromFiles(FILE *listofdistances, FILE *listofaligns, int *numberofgenes, char ***alltaxa, int *numtaxa,int distanceformat, int alignmentformat);
Treenode **getTreesFromFiles(FILE *, FILE *,char ****listoftaxanames, int ** listoftreesizes, int *numtrees, char ***alltaxa, int *numtaxa);
char **readPhylipAlignment(FILE *currentalignment, int *numspecies, char ***speciesnames);
char **readFastaAlignment(FILE *currentalignment, int *numspecies, char ***speciesnames);
char *getStringFromFile (FILE *datafile, char **title);
void printNexusFile(const Array2D<double> &solutions,const Array2D<int> &indexconversion, char **alltaxa, int numtaxa, const Array1D<int> &removedflags,double coeff, char *dirName);
void printRatesFile(const Array2D<double> &rates,distanceList *distlist,  char *dirName);
void printMatrixToFile(const Array2D<double> &coeffs_subset);
void printVectorToFile(const Array1D<double> &constants_subset);
void printMatrixtoMatlabformat(const Array2D<double> & currmat,FILE *outfile);


#endif
