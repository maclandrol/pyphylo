/* Copyright November 2005: Rachel Bevan
   Email: rachel@mcb.mcgill.ca
   NOTE: NOT Done commenting
*/


#ifndef distancelist_cpp
#define distancelist_cpp

#include"distanceList.h"
#include"wrapper.h"
#include<stddef.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

distanceList * newDistanceListnode(int numspecies,char **speciesnames){

  distanceList *newnode;
  int i;

  newnode=(distanceList *)Calloc(1,sizeof(distanceList));
  newnode->next=NULL;
  newnode->numspecies=numspecies;
  newnode->speciesnames=(char **)Calloc(numspecies,sizeof(char *));
 
  for(i=0;i<numspecies;i++){
    newnode->speciesnames[i]=(char*)Calloc((strlen(speciesnames[i])+1),sizeof(char));
    strcpy(newnode->speciesnames[i],speciesnames[i]);
   
  }
  newnode->genename=NULL;
  newnode->totaldistance=0.0;
  return newnode;
}

/* need to copy here - not working - causes a seg fault */
int setDistancesInDistanceList(distanceList *curr,double **dists){ 
  if(curr!=NULL){
    //for(i=0;i<curr->numspecies;i++){
    //  for(j=0;j<curr->numspecies;j++){
    //  
    // }
    //}
    curr->distances=dists;
    return 0;
  }
  else
    return -1;
}
double ** getDistancesFromDistanceList(distanceList *curr){  
  if(curr!=NULL){
   return curr->distances;
  }
  else
    return NULL;
}

int setVariancesInDistanceList(distanceList *curr, double **vars){

   if(curr!=NULL){
     curr->variances=vars;
    return 0;
  }
  else
    return -1;

}

double ** getVariancesFromDistanceList(distanceList *curr){

  if(curr!=NULL){
    return curr->variances;
  }
  else
    return NULL;

}


int setNextNodeInDistanceList(distanceList *curr, distanceList *next){  

  if(curr!=NULL){
    curr->next=next;
    return 0;
  }
  else
    return -1;
}

distanceList * getNextNodeFromDistanceList(distanceList *curr){  

  if(curr!=NULL)
    return curr->next;
  else
    return NULL;
}


int setNumSpeciesInDistanceList(distanceList *curr, int numspecies){

  if(curr!=NULL){
    curr->numspecies=numspecies;
    return 0;
  }
  else
    return -1;
  
}
int getNumSpeciesFromDistanceList(distanceList *curr){

  if(curr!=NULL)
    return curr->numspecies;
  else
    return -1;


}

int setSpeciesNamesInDistanceList(distanceList *curr, char **speciesnames){

  if(curr!=NULL){
    curr->speciesnames=speciesnames;
    return 0;
  }
  else
    return -1;

};
char **getSpeciesNamesFromDistanceList(distanceList *curr){

  if(curr!=NULL)
    return curr->speciesnames;
  else
    return NULL;

};

int setAlignmentNamesInDistanceList(distanceList *curr, char **alignnames,int numspecies){
  int i;

  if(curr!=NULL){
    curr->alignmentnames=(char **)Calloc(numspecies,sizeof(char *));
    for(i=0;i<numspecies;i++){
      curr->alignmentnames[i]=(char *)Calloc((strlen(alignnames[i])+1),sizeof(char));
      strcpy(curr->alignmentnames[i],alignnames[i]);

    }
    return 0;
  }
  else
    return -1;

};
char **getAlignmentNamesFromDistanceList(distanceList *curr){

   if(curr!=NULL)
    return curr->alignmentnames;
  else
    return NULL;


};
int setGeneNameInDistanceList(distanceList *curr, char* genename){

  if(curr!=NULL){
    curr->genename=(char*)Calloc(strlen(genename),sizeof(char*));
    strcpy(curr->genename,genename);
    return 0;
  }
  else
    return -1;

};

char *getGeneNameFromDistanceList(distanceList *curr){
  
  if(curr!=NULL){
    return curr->genename;
  }
  else
    return NULL;
};

int setTotalDistanceInDistanceList(distanceList *curr){
  int i,j;
  double **distances;
  double totaldistance=0.0;
  if(curr!=NULL){
    distances=curr->distances;
    if(distances!=NULL){
      for(i=0;i<curr->numspecies;i++){
	for(j=i+1;j<curr->numspecies;j++){
	  totaldistance=totaldistance+curr->distances[i][j];
	}
      }
      curr->totaldistance=totaldistance;
      return 0;
    }
    else
      return -1;
  }
  else
    return -1;

}

double getTotalDistanceFromDistanceList(distanceList *curr){

   if(curr!=NULL){
    return curr->totaldistance;
  }
  else
    return -1.0;
}

char **calculateTotalListOfSpeciesFromDistanceList(distanceList *head, int *numspecies){

  char **all_species=NULL;
  distanceList *curr=head;
  int currspecies;
  int *speciesfound;
  int totalspecies=0,newspecies=0;
  char **species_names;
  int species_index;
  int strlength,i,j;
  int speciesnotfound=0;
  int numgenesnotfound=0;
  int numgenes=0;

  while(curr!=NULL){
    speciesnotfound=0;
    currspecies=getNumSpeciesFromDistanceList(curr);
    newspecies=currspecies;
    speciesfound=(int *)Calloc(currspecies,sizeof(int));
    species_names=getSpeciesNamesFromDistanceList(curr);
    for(i=0;i<currspecies;i++){
      speciesfound[i]=0;
      for(j=0;j<totalspecies;j++){
	if(strcmp(species_names[i],all_species[j])==0){
    fprintf(stderr, "%s vs %s\n", species_names[i], all_species[j] );

	  speciesfound[i]=1;
	  j=totalspecies;
	  newspecies--;
	}
	
      }
      if(speciesfound[i]==0){
	speciesnotfound++;
      }
    }
    //fprintf(stderr, "speciesnotfound %i currspecies %i\n",speciesnotfound,currspecies);
    numgenes++;
    if(speciesnotfound==currspecies)
      numgenesnotfound++;

    species_index=totalspecies;
    totalspecies=totalspecies+newspecies;
    if(all_species==NULL)
      all_species=(char **)Calloc(totalspecies,sizeof(char *));
    else
      all_species=(char **)Realloc(all_species,(sizeof(char *)*(totalspecies-newspecies)),(sizeof(char *)*totalspecies));
    j=0;
    if(newspecies>0){
      for(i=species_index;i<totalspecies;i++){
	while(speciesfound[j]==1)
	  j++;
	strlength=strlen(species_names[j]);
	all_species[i]=(char*)Calloc((strlength+1),sizeof(char));
	strcpy(all_species[i],species_names[j]);
	j++;
      }
    }

    fprintf(stderr, "numgenesnotfound %i numgenes %i\n",numgenesnotfound,numgenes);
   
    free(speciesfound);
    curr=curr->next;
  }

  if(numgenes==numgenesnotfound){
    fprintf(stderr, "Error: Absolutely no overlap in species sets!\nExiting program...\n");
    exit(1);
  }
  *numspecies=totalspecies;
  
  return all_species;
}

#endif
