/* Copyright November 2005: Rachel Bevan
   Email: rachel@mcb.mcgill.ca
   NOTE: NOT Done commenting
*/

#ifndef distancelist_h
#define distancelist_h

#include<stdio.h>
#include<stddef.h>
//#include"tnt_array2d.h"
//#include"tnt_array2d_utils.h"

//using namespace TNT;
//using namespace JAMA;

typedef struct distnode{
  double **distances;
  double **variances;
  char **alignmentnames;
  int numspecies;
  char **speciesnames;
  distnode *next;
  char *genename;
  double totaldistance;
}distanceList;

distanceList * newDistanceListnode(int numspecies,char **speciesnames);
int setDistancesInDistanceList(distanceList *curr, double **dists);
double ** getDistancesFromDistanceList(distanceList *curr);
int setVariancesInDistanceList(distanceList *curr, double **vars);
double ** getVariancesFromDistanceList(distanceList *curr);
int setNextNodeInDistanceList(distanceList *curr, distanceList *next);
distanceList * getNextNodeFromDistanceList(distanceList *curr);
int setNumSpeciesInDistanceList(distanceList *curr, int numspecies);
int getNumSpeciesFromDistanceList(distanceList *curr);
int setSpeciesNamesInDistanceList(distanceList *curr, char **speciesnames);
char **getSpeciesNamesFromDistanceList(distanceList *curr);
int setAlignmentNamesInDistanceList(distanceList *curr, char **alignnames,int numspecies);
char **getAlignmentNamesFromDistanceList(distanceList *curr);
int setGeneNameInDistanceList(distanceList *curr, char* genename);
char *getGeneNameFromDistanceList(distanceList *curr);
int setTotalDistanceInDistanceList(distanceList *curr);
double getTotalDistanceFromDistanceList(distanceList *curr);
char **calculateTotalListOfSpeciesFromDistanceList(distanceList *head,int *numspecies);

#endif
