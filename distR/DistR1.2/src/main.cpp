/* Copyright November 2005: Rachel Bevan
   Email: rachel@mcb.mcgill.ca
*/
#include<stdio.h>
#include<stddef.h>
#include<stdlib.h>
#include <sys/stat.h>  
#include"wrapper.h"
#include<math.h>
#include <time.h>
#include <unistd.h>
#include "phylogeny.h"
#include <string.h>
#include "distR.h"
#include"matrixOperations.h"
#include"fileio.h"
#include"distanceList.h"
#include <sys/time.h>

void helpFull ();
void readCommandLine(int argc, char **argv);
char myGetOpt(int *currargc, int argc, char **argv, char *options);
void FreeData(distanceList *distmats, int numtrees, char **alltaxa, int numtaxa);
char myoptarg[MAXSTRSIZE];

int main(int argc, char **argv){

  /* Read in arguments from the commandline */
  readCommandLine(argc, argv);
  return 0;
}



void helpFull ()
{
  char *usage =  "Usage: [options] \n"
"Options:\n"
"     -h        This help screen\n"
"     -t File   Specify file which contains a list of tree files in NEWICK format\n"
"               Must use either this option, or the -d, or -b options to specify\n"
"               trees\n"
"                     \n"
"     -d File   Specify file which contains a list of distance matrix files in\n"
"               NEXUS format\n"
"               Must use either this option, or the -t, or -b options to specify\n"
"               trees.\n"
"                     \n"
"     -p File   Specify file which contains a list of alignment files in PHYLIP\n"
"               format - both interleaved and sequential are acceptable\n"
"               Must use either this option, or the -n, or -b options to specify\n"
"               alignments.\n"
"                      \n"
"     -n File   Specify file which contains a list of alignment files in NEXUS\n"
"                format\n"
"               Must use either this option, or the -p, or -b options to specify\n"
"               alignments.\n"
"                      \n"
"               Each line in the list of alignments file, contains the file name\n"
"               of the alignment corresponding to the tree/distance matrix in the\n"
"               appropriate file\n"
"     -b File   Specify file which contains a list of NEXUS files that have both\n"
"               a distance matrix and an alignment\n"
"               This option allows for specification of both tree distances and\n"
"               alignment in same file.  Note:  I have removed this option due to\n"
"               bugs in the code with more complex nexus files.  If you wish to \n"
"               use this option simply uncomment the code in the switch statement\n" 
"               under option b.\n"
"                      \n"
"     -o  Dir   Specify the output directory\n"
"                      \n"
"     Please note that if no alignment files are provided the program defaults to\n"
"     PHYLIP formatted alignment files which the file extension (i.e. the last '.*')\n"
"     of the tree/distance files specified by -t/-d option are changed to .phy\n"
"                      \n"
"     Also note that the user MUST specify either a list of trees in newick format\n"
"     with the -t option OR a list of distance matrices in NEXUS format with the -d\n"
"     option\n"
    "                      \n";

  fprintf(stderr, "%s\n", usage);
  
  return;
}


void readCommandLine(int argc, char **argv){

  char *progname = NULL;         /* The name of the executable program */
  char c;                        /* For reading options */
  char *dirName = NULL;          /* Output directory name */
  int commandcount = 1;          /* Keep track of the number of command-line args */
  int numtrees=0;                /* Number of trees specified by list of trees file */
  char **alltaxa=NULL;           /* Names of all the taxa over all the trees */
  int numtaxa=0;                 /* Number of taxa over all the trees */
  FILE *listofalignsfile=NULL;   /* File to specify the names of the alignment files */
  FILE *listofdistances=NULL;    /* File to specify the names of the tree/distance files */
  distanceList *distmats;        /* Linked list structure that holds the distance matrix
				    and alignment for each tree (among other things) */
  /* These read flags are used to check for errors to ensure
     that the user doesn't use the incorrect options together */
  int readnewickdistances=0;     
  int readnexusdistances=0;
  int readnexusalignments=0;
  int readphylipalignments=0; 
  int readfastaalignments=0;
  int readnexusboth=0;
  /* The following two flags are used to tell the read function that
     gets the distances and alignments from the correct files, what the
     format of the distance and alignment files is */
  int distanceformat=-1;
  int alignmentformat=-1;
  struct timeval tvals[2];
  struct stat sb;


  /* This command gets the program name from the first command line argument, without
     the directory specification */
  if ((strrchr (argv[0], '/')) == NULL)
    progname = argv[0];
  else
    progname = "null";
  
  /* While there are more options to get */
   while (1)
    {
      /* Get the next option. */
      c = myGetOpt(&commandcount,argc, argv, "ht:p:d:n:bf:o:");
      // c = myGetOpt(&commandcount,argc, argv, "ht:p:d:n:b:f:");
      /* Make sure the option is valid, if not output help, otherwise
	 exit while loop since finished reading command line args */
      if (c == -1){
	if(commandcount == 0)
	  helpFull();
	break;
      }
      
      switch (c)
	{
	  
	  /* Print help */
	case 'h':
	  fprintf(stdout, "Help\n");
	  helpFull ();
	  exit(1);
	  break;
	
	  /* option to specify that each file in list of files contains both alignment 
	     and distance matrix in one nexus file */
	case 'b':

	  fprintf(stderr, "Error: This option is currently not usable due to bugs in the code with more complex nexus files.  If you wish to use this option with VERY SIMPLE nexus files (like Buttercups.nex in the bin directory) then please uncomment the code in the switch statement under case b, and uncomment the second myGetOpt command.\n Exiting program...\n");
	 //  listofdistances=fopen(myoptarg, "r");
// 	  /* make sure that file opened properly */
// 	  if(listofdistances==NULL){
// 	    fprintf(stderr, "ERROR: could not open file %s, exiting program...\n",
// 		    myoptarg);
// 	    exit(1);
// 	  }
// 	  readnexusboth=1;
	  
	  break;
	
  case 'o':
    if (access(myoptarg, 0) != 0){
      fprintf(stderr, "ERROR: directory %s not accessible\n", myoptarg);
      exit(1);
    }
    if(stat(myoptarg, &sb)==0 && S_ISDIR(sb.st_mode)){
        dirName = myoptarg;
    }
    break;

	  /* option to specify that each file in list of alignments is in nexus format */
	case 'n':
	  listofalignsfile = fopen(myoptarg, "r");
	  /* make sure that file opened properly, otherwise exit */ 
	  if(listofalignsfile==NULL){
	    fprintf(stderr, "ERROR: could not open file %s, exiting program...\n",
		    myoptarg);
	    exit(1);
	  }
	  //fprintf(stderr, "alignment file name: %s",myoptarg);
	  readnexusalignments=1;
	  alignmentformat=1;
	  break;

	  /* option to specify that each file in list of alignments is in phylip format
	     both interleaved and sequential work */
	case 'p':
	  listofalignsfile = fopen(myoptarg, "r");
	  /* make sure that file opened properly, otherwise exit */ 
	  if(listofalignsfile==NULL){
	    fprintf(stderr, "ERROR: could not open file %s, exiting program...\n",
		    myoptarg);
	    exit(1);
	  }
	  readphylipalignments=1;
	  alignmentformat=0;
	  break;
	
	  /* option to specify that each file in list of trees/distance matrices is in
	     nexus format */
	case 'd':
	  listofdistances=fopen(myoptarg, "r");
	  /* make sure that file opened properly, otherwise exit */ 
	  if(listofdistances==NULL){
	    fprintf(stderr, "ERROR: could not open file %s, exiting program...\n",
		    myoptarg);
	    exit(1);
	  }
    fprintf(stderr, "Distance file read\n");
	  readnexusdistances=1;
	  distanceformat=1;
	  break;

	  /* option to specify that each file in list of trees/distance matrices is in
	     newick format */
	case 't':
	  listofdistances = fopen(myoptarg, "r");
	  /* make sure that file opened properly, otherwise exit */ 
	  if(listofdistances==NULL){
	    fprintf(stderr, "ERROR: could not open file %s, exiting program...\n",
		    myoptarg);
	    exit(1);
	  }
	  readnewickdistances=1;
	  distanceformat=0;
	  break;

	case 'f':
	  listofalignsfile = fopen(myoptarg, "r");
	  /* make sure that file opened properly, otherwise exit */ 
	  if(listofalignsfile==NULL){
	    fprintf(stderr, "ERROR: could not open file %s, exiting program...\n",
		    myoptarg);
	    exit(1);
	  }
	  readfastaalignments=1;
	  alignmentformat=2;
	  break;
	  
	default:
	  fprintf(stdout, "Help\n");
	  helpFull ();
	  exit(1);
	  break;
	}
    }

   /* if trees/distance matrices have not been specified print error message
      and exit program, or if both trees and distance matrices have 
      been specified then print an error message and exit program */
   if((readnewickdistances!=1)&&(readnexusdistances!=1)){
    fprintf(stderr, "ERROR: Must specify either a list of trees (in NEWICK format), or a list of distance matrices (in NEXUS format) to analyze. Exiting program...\n");
    exit(1);
    
  }
  else if((readnewickdistances==1)&&(readnexusdistances==1)){
    fprintf(stderr, "ERROR: Must specify either a list of trees (in NEWICK format), or a list of distance matrices (in NEXUS format) to analyze, not both. Exiting program...\n");
    exit(1);
  }
   
   /* If both types of alignment format have been specified then print message and exit 
      program */
   if((readnexusalignments==1)&&(readphylipalignments==1)){
     fprintf(stderr, "ERROR: Must specify either a list of alignments in PHYLIP format, or a list of alignments in NEXUS format - not BOTH.  If neither is specified the program with default to PHYLIP formatted alignments with the same name as the tree/distance files with a .phy extension.  Exiting program...\n");
     exit(1);
   }

   if((readnexusalignments==1)&&(readfastaalignments==1)){
     fprintf(stderr, "ERROR: Must specify either a list of alignments in FASTA format, or a list of alignments in NEXUS format - not BOTH.  If neither is specified the program with default to PHYLIP formatted alignments with the same name as the tree/distance files with a .phy extension.  Exiting program...\n");
     exit(1);
   }
   if((readfastaalignments==1)&&(readphylipalignments==1)){
     fprintf(stderr, "ERROR: Must specify either a list of alignments in PHYLIP format, or a list of alignments in FASTA format - not BOTH.  If neither is specified the program with default to PHYLIP formatted alignments with the same name as the tree/distance files with a .phy extension.  Exiting program...\n");
     exit(1);
   }
   /* If using the option to have both distances and alignment in one nexus file, then
      make sure that no other option was specified since only one list of files should
      be given */
   if(readnexusboth==1){
     if((readnexusalignments==1)||(readphylipalignments==1)||(readnewickdistances!=1)||(readnexusdistances!=1)){
       fprintf(stderr, "ERROR: When using the -b option no other files should be specified since both alignments and distances should be specified in the same file in NEXUS format. Exiting program...\n");
       exit(1);
     }
     readnexusalignments=1;
     readnexusdistances=1;
     distanceformat=2;
   }
   
   /* If there is a list of trees/distance matrices to read then get the input from file
      and calculate the tree rates */
   if((readnewickdistances==1)||(readnexusdistances==1)){
     //(void)time(&t_beg);
     gettimeofday(&(tvals[0]),NULL);
     distmats=readInputFromFiles(listofdistances, listofalignsfile, &numtrees, &alltaxa, &numtaxa,distanceformat, alignmentformat);
     fprintf(stderr, "Number of Trees %i\nNumber of taxa %i\n", numtrees,numtaxa);
    //  if(numtaxa>150){
//        fprintf(stderr, "Too many taxa (>150).  Exiting program...\n");
//        exit(1);
//      }
     calculateTreeRates(distmats, numtrees, alltaxa, numtaxa, dirName);
     FreeData(distmats,numtrees,alltaxa,numtaxa);
     gettimeofday(&(tvals[1]),NULL);
     //(void)time(&t_end);     
     //hour = div((int)(t_end-t_beg),3600);
     //min  = div((int)(t_end-t_beg),60  );
     //min.quot -= hour.quot*60;
     //fprintf(stdout,"\n\n. Time used %dh%dm%ds\n", hour.quot,min.quot,(int)(t_end-t_beg)%60);
     //fprintf(stdout, "\n\n. Time used %.15lf\n", difftime(t_end,t_beg));
     //fprintf(stdout, "Time start %d, Time end %ld\n",tvals[0].tv_usec,tvals[1].tv_usec);
     fprintf(stdout, "Time used seconds: %ld microseconds %ld\n",(tvals[1].tv_sec-tvals[0].tv_sec),(tvals[1].tv_usec-tvals[0].tv_usec));
   }
   return;
}


char myGetOpt(int *currargc, int argc, char **argv, char *options){
  
  char retchar;     /* character option */
  char *findchar;   /* used to find retchar in the list of possible options */

 
  /* each option must start with a -.  If random junk is put on the
     command line, then print an error message and exit program. */
 
  
  /* Check to make sure haven't exceeded the number of command line args, if have
     return -1 to signal to calling function that have finished parsing command line
     options */
  if((*currargc)<argc){
    /* Get the first character after - */
    if((strchr(argv[*currargc],'-')!=NULL)&&(strlen(argv[*currargc])==2)){
      retchar=argv[*currargc][1];
      /* Find retchar in the list of possible options */
      findchar=strchr(options,retchar);
      /* Move to the next command line arg */
      (*currargc)++;
      /* if the character is found in string options, then it is a valid option,
	 so parse further to determine if there is an argument to the option.  Otherwise
	 the option is not valid, so return retchar and let main deal with it */
      if(findchar!=NULL){
	/* If the character is followed by a : in the options string, then
	 it means that the next command line arg is to be used by this option.
         If not then just copy a blank string into myoptarg */
	if(findchar[1]==':'){
	  /* Make sure that the next command line arg is actually specified by the user
	     otherwise print an error message and exit */
	  if(((*currargc)<argc)){
	    /* Make sure that the next command line arg is not a new option, and if not
	       copy the string.  Otherwise print an error message and exit */
	    if(strchr(argv[(*currargc)],'-')==NULL){
	      strcpy(myoptarg,argv[(*currargc)]);
	      (*currargc)++;
	    }
	    else{
	      fprintf(stderr, "ERROR: missing argument for -%c option.  Please see -h option for further details. Exiting program\n",retchar);
	      exit(1);
	    }
	  }
	  else{
	    fprintf(stderr, "ERROR: incorrect commandline options, missing argument for -%c option.  Please see -h option for further details. Exiting program...\n",retchar);
	    exit(1);
	  }
	}
	else
	  strcpy(myoptarg,"");
	return retchar;
      }
      else{
	return retchar;
      }
    }
    else{
      fprintf(stderr, "ERROR: incorrect command line option.  Options are specified by -c, where c is a given character. Please see -h option for further details. Exiting program...\n");
      exit(1);
    }
  }
  else
    return -1;
  
}

void FreeData(distanceList *distmats, int numtrees, char **alltaxa, int numtaxa){

  int i,j;
  int currnumspecies;
  distanceList * currnode=distmats,*next;
  
  //fprintf(stderr, "free step 1\n");
  for(i=0;i<numtaxa;i++){
    Free(alltaxa[i],(strlen(alltaxa[i])*sizeof(char)));

  }
  Free(alltaxa,(numtaxa*sizeof(char*)));
  //fprintf(stderr, "free step 2\n");
  for(i=0;i<numtrees;i++){
    currnumspecies=currnode->numspecies;
    //fprintf(stderr, "free variance and distance matrices\n");
    for(j=0;j<currnumspecies;j++){
      // fprintf(stderr, "dist/vars %i\n",j);
      Free(currnode->variances[j],(currnumspecies*sizeof(double)));
      Free(currnode->distances[j],(currnumspecies*sizeof(double)));
      // fprintf(stderr, "alignments %i,strlen %i\n",j,strlen(currnode->speciesnames[j]));
      Free(currnode->speciesnames[j],(strlen(currnode->speciesnames[j])*sizeof(char)));
     //  fprintf(stderr, "free alignment names\n");
//       Free(currnode->alignmentnames[j],(strlen(currnode->alignmentnames[j])*sizeof(char)));
      
    }
    //fprintf(stderr, "final free step\n");
    Free(currnode->variances,(currnumspecies*sizeof(double *)));
    Free(currnode->distances,(currnumspecies*sizeof(double *)));
    Free(currnode->speciesnames,(currnumspecies*sizeof(char *)));
    Free(currnode->genename,(strlen(currnode->genename)*sizeof(char)));
    next=currnode->next;
    //fprintf(stderr, "free distance list\n");
    Free(currnode,sizeof(distanceList));
    currnode=next;
  }
}
