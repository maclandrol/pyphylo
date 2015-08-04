/* Copyright November 2005: Rachel Bevan
   Email: rachel@mcb.mcgill.ca
   NOTE: NOT Done commenting
*/

#ifndef fileio_c
#define fileio_c

#include"fileio.h"
#include<stdio.h>
#include <string>
#include"distR.h"
#include"simple_nexus.h"

#define MAXNAMELEN 200
/* This function reads a nexus file, using code by david bryant simple_nexus.cpp
   The results are parsed according to whether or not the calling function wants
   distances, an alignment or both returned */
   int readNexusFormat(char *infile, int readnexusdistances,double *** distances, int *numspecies, char ***speciesnames,char ***alignment,int readnexusalignment){

  int i,j;                            /* counters */
  /* The following variables are necessary to call the function that reads nexus
     format */
    string s;
    bool success;
    vector<string> taxa_names;
    vector<sequence> seqs;
    data_types seq_type;
    vector<string> char_names;
    bootstrap_params bs_params;
    dbl_array D;
    st_assumpt_params sta_params;
  /* The following variables are to hold the alignment and species names
     which need to be taken out of the vector seqs and taxa_names respectively */
    char **myspeciesnames;
    char **curralignment;
    double **currdistances;

  /* Default assumptions */
    sta_params.power=0;
    sta_params.cutoff = 0.0;
    sta_params.constrained=true;
    sta_params.transform_type = NONE;
    sta_params.equal_rates = true;
    sta_params.shape = 0.0;
    sta_params.pinvar = 0.0;
    sta_params.dist_output = "";
    sta_params.switch_rate = 0.0;
    sta_params.ss_weight = 0.0;

  /* Open up the input file and read nexus format */
    ifstream is(infile);
    success = read_nexus_file(is,s,taxa_names,seqs,seq_type,char_names,D,sta_params,bs_params);
  /* If the file has been successfully read then copy distances, alignment and sequence names
     according to the flags readnexusdistance and readnexusalignment */

    if(success){
    /* Get the number of species and species names regardless of whether returning distances
       alignment or both */
      *numspecies=taxa_names.size();
      (myspeciesnames)=(char **)Calloc((*numspecies),sizeof(char *));
      for(i=0;i<(*numspecies);i++){
        (myspeciesnames)[i]=(char *)Calloc((taxa_names[i].size()+1),sizeof(char));
        strcpy((myspeciesnames)[i],taxa_names[i].data());
      }
      *speciesnames=myspeciesnames;
    /* If the distances are to be returned check that they were read from the file.  If
       not then print error message and exit program.  If so then copy the distances from
       the dbl_array data type to the Array2D data type */
       if(readnexusdistances){
        if (D.empty()) {
         fprintf(stderr, "ERROR: No distance matrix provided in %s. Exiting program...\n",infile);
         exit(1);
       }
       else{    

        currdistances=(double **)Calloc((*numspecies),sizeof(double *));
        for(i=0;i<(*numspecies);i++){
          currdistances[i]=(double *)Calloc((*numspecies),sizeof(double ));
          for(j=0;j<(*numspecies);j++){
           currdistances[i][j]=D[i][j];
         }
       }

       (*distances)=currdistances;
     }

   }
    /* If the alignment is to be returned check that they were read from the file.  If
       not then print error message and exit program.  If so then copy the alignment from
       the vector data type to the char ** data type */
       if(readnexusalignment){
        if(!seqs.empty())
        {
         curralignment=(char **)Calloc((*numspecies),sizeof(char *));
         for(i=0;i<(*numspecies);i++){
	    /* need to get sequence value short int by short int and convert to char */
           curralignment[i]=(char *)Calloc((seqs[i].size()+1),sizeof(char));
           for(j=0;j<(*numspecies);j++){
             curralignment[i][j]=get_state(seq_type, seqs[i][j]);

           }
           curralignment[i][seqs[i].size()+1]='\0';
         }
         (*alignment)=curralignment;
       }

       else{

         fprintf(stderr, "ERROR: Could not read sequences from nexus file. Exiting program...\n");
         exit(1);
       }
     }
     return 0;
   }
   else{
    fprintf(stderr, "ERROR: Cannot read file %s. Exiting program...\n",infile);
    fprintf(stderr, "Error is the following %s.\n",s.c_str());
    exit(1);
    
  }
  
}


/* Function responsible for controling how data is read from input files specified, depending
   upon the flags distanceformat and alignment format.  Data is read into a linked distanceList
   which contains the distance matrix and alignment among other info.  Also the function 
   checks to make sure that the number of species in the alignment and the tree/distance matrix
   match etc. to ensure that the data has no errors (as best possible) before running the main
   function calculateTreerates */
   distanceList *readInputFromFiles(FILE *listofdistances, FILE *listofaligns, int *numberofgenes, char ***alltaxa, int *numtaxa,int distanceformat, int alignmentformat){

  distanceList *distanceMats,*curr,*next;  /* For building the list of data */
  char c;                                  /* To read chars from file */
  int treecountfromfile=0;                 /* To count the number of tree files in tree input
					      file specified by user */
  int aligncountfromfile=0;                /* To count the number of alignment files in
					      alignment input file specified by user */
  int distmatcountfromfile=0;              /* To count the number of distances matrix files
					      in distance matrix input
					      file specified by user */
  fpos_t startlistofaligns,startlistofdistances;    /* Used to count number of tree/distance
						       files and alignment files and make
						       sure same */
  /* The following 3 vars are used to read in the names of the tree file, aligment file
     and distance matrix file respectively, from the user specified list of files */
                   char treename[MAXSTRSIZE],alignname[MAXSTRSIZE],distname[MAXSTRSIZE];
  /* These are the current files from the lists of files specified by the user */
                   FILE *currenttree=NULL;
                   FILE *currentalignment=NULL;
                   FILE *currentdistmat=NULL;

  char **taxa_names;                 /* Taxa names for current tree/distance matrix */
  int numspecies;                    /* Number of species for current tree/distance matrix */
  char **speciesnames;               /* Taxa names for current alignment */
  int alignnumspecies=0;             /* Number of species for current tree/distance matrix */
  char **alignment=NULL;             /* Current alignment read from currentalignment */
  Treenode *mytree=NULL;             /* Used to hold newick tree read from file */
  int startlist=0;                   /* A flag to indicate if distanceList has a node in it or not */
  char *genename;                    /* Name of current gene, put in distanceList */
                   char *tempstr;

  /* The following three read flags are used to indicate whether alignments need to be read from
     separate files or not, and alignments/distances are to be read from nexus format or not.  Necessary
     for control to determine file format being read */
     int readalignments=0;              
     int readnexusdistances=0, readnexusalignments=0;
     double **distances;
     double **variances;
     int i,l;

     if(listofdistances==NULL){
      fprintf(stderr, "ERROR: List of distances (as trees in newick format, or as distances in nexus format) has NOT been opened properly. Exiting program...\n");
      exit(1);
    }
  /* Initialize the number of genes read from file to 0 */
    (*numberofgenes)=0;
  /* Check that a list of alignments has been specified, if so count the
     number of alignments.  This value will be compared to the number of tree/distance matrices
     to ensure the two are equal */
     if(alignmentformat!=-1){
    /* Make sure opening the list of alignment files properly */
      if(listofaligns==NULL){
        fprintf(stderr, "ERROR: List of alignments has NOT been opened properly. Exiting program...\n");
        exit(1);
      }
    /* Count number of alignment files given */
      fgetpos(listofaligns,&startlistofaligns);
      while((c=fgetc(listofaligns))!=EOF){
        ungetc(c, listofaligns);
        fscanf(listofaligns, "%s\n",alignname);
      //fprintf(stderr, "alignname %s\n",alignname);
        aligncountfromfile++;
      }
      fsetpos(listofaligns,&startlistofaligns);
    }

  /* If distance format is 0 then reading tree from newick format */
    if(distanceformat==0){
    /* Compute number of trees specified */
      fgetpos(listofdistances,&startlistofdistances);
      while((c=fgetc(listofdistances))!=EOF){
        ungetc(c, listofdistances);
        fscanf(listofdistances, "%s\n",treename);
        treecountfromfile++;

      }
      fsetpos(listofdistances,&startlistofdistances);
    /* Make sure number of alignments and number of trees agree if a list of alignments has
       been specified.  If a list of alignments hasn't been specified then don't worry because the
       names of the tree files are used to get alignment file names */
       if((alignmentformat!=-1)&&(treecountfromfile!=aligncountfromfile)){
        fprintf(stderr, "ERROR: Number of Alignments does not correspond to number of trees.\nPlease check input files.  Exiting program...\n");
        exit(1);
      }

    /* Calculate default names of alignment files to open, just in case it is needed while opening tree
       files and reading in trees */

      while(((c=fgetc(listofdistances))!=EOF)&&(c!='\n')){
        ungetc(c, listofdistances);
        fscanf(listofdistances, "%s\n",treename);
        currenttree=fopen(treename, "r");
      /* Make sure tree file can be opened properly */
        if(currenttree==NULL){
         fprintf(stderr, "ERROR: Not opening treefile %s properly. Exiting program...\n", treename);
         exit(1);
       }

       /* Read in tree and calculate distances */
       mytree=readPhylipTree(currenttree, &(taxa_names), &(numspecies));
       fclose(currenttree);

       next=newDistanceListnode(numspecies,taxa_names);

       if(startlist==0){
         distanceMats=next;
         curr=distanceMats;
         startlist=1;
         (*numberofgenes)++;
       }
       else{
         curr->next=next;
         curr=curr->next;
         (*numberofgenes)++;
       }

       distances=calculateDistanceMatrixfromTree(mytree,numspecies);
       removeTree(mytree);
       for(l=0;l<numspecies;l++){
         Free(taxa_names[l],(strlen(taxa_names[l])*sizeof(char)));

       }
       Free(taxa_names,(numspecies*sizeof(char *)));
       setDistancesInDistanceList(curr, distances);
       setTotalDistanceInDistanceList(curr);
      /* If the name of the tree file has a . then copy the name upto the first . RACHEL NOTE: CHECK
         that this works if have > 1 . in file name! */
       setGeneNameInDistanceList(curr,treename);



      /* Now need to read in alignment */
     }

     readalignments=1;

   }
   else if((distanceformat==1)||(distanceformat==2)){

    fgetpos(listofdistances,&startlistofdistances);
    while((c=fgetc(listofdistances))!=EOF){
      ungetc(c, listofdistances);
      fscanf(listofdistances, "%s\n",distname);
      distmatcountfromfile++;
      
    }
    fsetpos(listofdistances,&startlistofdistances);

    if((alignmentformat!=-1)&&(distmatcountfromfile!=aligncountfromfile)){
      fprintf(stderr, "Number of Alignments does not correspond to number of distance matrices.\nPlease check input files.  Exiting program...\n");
      exit(1);
    }


    while((c=fgetc(listofdistances))!=EOF){
     ungetc(c, listofdistances);
     fscanf(listofdistances, "%s\n",distname);
     currentdistmat=fopen(distname, "r");
     if(currentdistmat==NULL){
      fprintf(stderr, "ERROR: Not opening distance matrix file %s properly. Exiting program...\n",distname);
      exit(1);
    }
    fclose(currentdistmat);

       /* Here I want to read the nexus files, and create a 
	  new distanceList node */

       /* should fix this function to take in flags of whether or not
	  to read (and return) distances, alignments or both */
    if(distanceformat==2){
      readnexusdistances=1;
      readnexusalignments=1;
      readalignments=0;
    }
    else{
      readnexusdistances=1;
      readalignments=1;
    }
    readNexusFormat(distname,readnexusdistances,&distances,&numspecies, &taxa_names,&alignment,readnexusalignments); 
    next=newDistanceListnode(numspecies,taxa_names);
    if(startlist==0){
      distanceMats=next;
      curr=distanceMats;
      startlist=1;
      (*numberofgenes)++;
    }
    else{
      curr->next=next;
      curr=curr->next;
      (*numberofgenes)++;
    }

    setDistancesInDistanceList(curr, distances);
    setTotalDistanceInDistanceList(curr);
    setGeneNameInDistanceList(curr,distname);

    if(readalignments==0){
	 /* MODALIGNMENT!*/
	 //fprintf(stderr, "read alignment 0, calculating variances\n");
      variances=calculateVarianceMatrixfromAlignment(distances, alignment, taxa_names, taxa_names, numspecies);
      setVariancesInDistanceList(curr,variances);
      for(i=0;i<numspecies;i++)
        free(alignment[i]);
      free(alignment);
	 /* Alignment AND distances have been read from Nexus file so need to fill in
	    rest of distanceList node */

    }


  }
}

  //fprintf(stderr, "time to read alignments!\n");
curr=distanceMats;
int curralign=0;
while(curr!=NULL){
    //fprintf(stderr,"read curralign %i\n",curralign);
  if(readalignments==1){

    if(alignmentformat==-1){
     tempstr=getGeneNameFromDistanceList(curr);
     strcpy(alignname,tempstr);
	//fprintf(stderr, "in if alignname %s\n",alignname);
     if((tempstr=strrchr(alignname,'.'))!=NULL){
       genename=(char*)Calloc((strlen(alignname)+5),sizeof(char));
       strncpy(genename,alignname, (strlen(alignname)-strlen(tempstr)));
     }
     else{
       genename=(char*)Calloc((strlen(alignname)+5),sizeof(char));
       strcpy(genename,alignname);
     }
     strcat(genename, ".phy");
     genename[strlen(genename)]='\0';
	//	fprintf(stderr, "genename %s\n",genename);
     currentalignment=fopen(genename, "r");

     if(currentalignment==NULL){
       fprintf(stdout, "ERROR: Not opening alignment file '%s' properly. Exiting program...\n",genename);
       exit(1);
     }
     Free(curr->genename,(strlen(curr->genename)*sizeof(char)));
     setGeneNameInDistanceList(curr,genename);
     alignment=readPhylipAlignment(currentalignment,&alignnumspecies,&speciesnames);
   }
   else{
     fscanf(listofaligns, "%s\n", alignname);
	//fprintf(stderr, "in else alignname %s\n",alignname);
     alignname[strlen(alignname)]='\0';
     if((alignmentformat==0)||(alignmentformat==2)){
       currentalignment=fopen(alignname, "r");
       if(currentalignment==NULL){
         fprintf(stdout, "ERROR: Not opening alignment file '%s' properly. Exiting program...\n",genename);
         exit(1);
       }

       if(alignmentformat==0){
         alignment=readPhylipAlignment(currentalignment,&alignnumspecies,&speciesnames);
       }
       else if(alignmentformat==2){
         alignment=readFastaAlignment(currentalignment,&alignnumspecies,&speciesnames);
       }
       fclose(currentalignment);
     }
     else if(alignmentformat==1){
       readnexusdistances=0;
       readnexusalignments=1;

       readNexusFormat(genename,readnexusdistances,&distances,&alignnumspecies, &speciesnames,&alignment,readnexusalignments); 

     }
     Free(curr->genename,(strlen(curr->genename)*sizeof(char)));
     setGeneNameInDistanceList(curr,alignname);
   }
   if(compareSpeciesLists(speciesnames,alignnumspecies,getSpeciesNamesFromDistanceList(curr),getNumSpeciesFromDistanceList(curr))==-1){
     fprintf(stderr, "ERROR: Number of species in alignment doesn't match number of species in tree or distance matrix for gene %s\nPlease check input files, exiting program....\n",getGeneNameFromDistanceList(curr));
     exit(1);
   }
      //fprintf(stderr, "read alignment 1 calculating variances\n");
   variances=calculateVarianceMatrixfromAlignment(getDistancesFromDistanceList(curr), alignment, speciesnames, getSpeciesNamesFromDistanceList(curr), getNumSpeciesFromDistanceList(curr));
      //fprintf(stderr,"set variances in curr alignnumspecies %i\n",alignnumspecies);
   setVariancesInDistanceList(curr,variances);
      //fprintf(stderr,"free alignment\n");
   for(i=0;i<alignnumspecies;i++){
      	//fprintf(stderr, "free alignment %i strlen %i %s\n",i,strlen(alignment[i]),alignment[i]);
     Free(alignment[i],(strlen(alignment[i])*sizeof(char)));
     Free(speciesnames[i],(MAXNAMELEN*sizeof(char)));
	//fprintf(stderr, "mem freed\n");

   }

   Free(alignment,(alignnumspecies*sizeof(char *)));
   Free(speciesnames,(alignnumspecies*sizeof(char *)));
      //fprintf(stderr,"done freeing alignment\n");
 }
    // fprintf(stderr,"get next distancelist\n");
 curr=curr->next;
 curralign++;
}
  //exit(1);
  //fprintf(stderr, "done reading alignments\n");
(*alltaxa)=calculateTotalListOfSpeciesFromDistanceList(distanceMats,numtaxa);

  //fprintf(stderr, "Done reading distance matrices\n");
return distanceMats;


}



char **readPhylipAlignment(FILE *currentalignment, int *numspecies, char ***speciesnames){

  char **alignment=NULL;
  int i;
  int length;
  char c,d;
  int speciescount=0;
  int counter=0;
  int placeinstring=0;
  int doneheader=0;
  fpos_t startoffile;

  fgetpos(currentalignment,&startoffile);
  while((c=fgetc(currentalignment))!='\n'){
    if(!(((c>='0')&&(c<='9'))||(c==' ')||(c=='\t')||(c=='\r'))){
      fprintf(stderr, "ERROR: Incorrect char %c in first line of phylip file.\nPlease check format of alignment\nExiting program...\n;", c);
      exit(1);
    }
    
  }

  fsetpos(currentalignment, &startoffile);
  fscanf(currentalignment, "%i\t%i\n", numspecies, &length);
  (*speciesnames)=(char**)Calloc((*numspecies),sizeof(char*));
  alignment=(char**)Calloc((*numspecies),sizeof(char*));

  for(i=0;i<*numspecies;i++){
    (*speciesnames)[i]=(char *)Calloc(MAXNAMELEN,sizeof(char));
    alignment[i]=(char *)Calloc((length+2),sizeof(char));
  }
  
  while((c=fgetc(currentalignment))!=EOF){
    ungetc(c, currentalignment);
    
    if(doneheader==0){
      while(speciescount!=*numspecies){

       fscanf(currentalignment, "%s", (*speciesnames)[speciescount]);
	//fprintf(stderr, "speciesname %s\n",(*speciesnames)[speciescount]);

       if((d=fgetc(currentalignment))==EOF){
         fprintf(stderr, "ERROR: Number of sequences in file does not match specified number of sequences. Please check input alignment\nExiting program...\n");
         exit(1);
       }



       counter=0;
       while(((d=fgetc(currentalignment))!='\n')&&(d!='\r')&&(d!=EOF)){
         if(d!=' '){
           alignment[speciescount][counter]=d;
           counter++;
	    //fprintf(stderr, "char d %c\n",d);
         }
       }
       if((d==EOF)&&(*numspecies!=speciescount)){
         fprintf(stderr, "ERROR: Number of sequences in file does not match specified number of sequences. Please check input alignment\nExiting program...\n");
         exit(1);
       }
       speciescount++;
     }
     placeinstring=counter;
     doneheader=1;
   }
   else{
      // fprintf(stderr, "In else should be done reading!\n");
    speciescount=0;
    while(((d=fgetc(currentalignment))==' '));
    d=fgetc(currentalignment);
    while(speciescount!=*numspecies){
     counter=placeinstring;
     while(((d=fgetc(currentalignment))!='\n')&&(d!='\r')&&(d!=EOF)){
       if(d!=' '){
         alignment[speciescount][counter]=d;
         counter++;
       }
     }
     speciescount++;
   }
   placeinstring=counter;

 }


}

if(speciescount!=*numspecies){
  fprintf(stderr, "Number of sequences in file does not match specified number of sequences. Please check input alignment\nExiting program...\n");
  exit(1);

}

for(i=0;i<*numspecies;i++){
  alignment[i][placeinstring]='\0';
    //fprintf(stderr,"alignment ptr %i\n",alignment[i]);
    // fprintf(stderr, "Strlen alignment %i\n",strlen(alignment[i]));
    // fprintf(stderr, "Speciesname: %s\nAlignment: %s\n",(*speciesnames)[i],alignment[i]);

}


return alignment;
};

char **readFastaAlignment(FILE *currentalignment, int *numspecies, char ***speciesnames){

  char **alignment;
  int speciescount=0;
  fpos_t startpos;
  char *curralignname;
  int i,strcount=0;
  char c;
  
  fgetpos(currentalignment,&startpos);
  while((c=fgetc(currentalignment))!=EOF){
    if(c=='>')
      speciescount++;
  }
  alignment=(char **)Calloc(speciescount,sizeof(char *));
  (*speciesnames)=(char **)Calloc(speciescount,sizeof(char*));
  for(i=0;i<speciescount;i++){
    (*speciesnames)[i]=(char *)Calloc(MAXNAMELEN,sizeof(char));
  }
  fsetpos(currentalignment,&startpos);
  while((c=fgetc(currentalignment))!=EOF){
    ungetc(c,currentalignment);
    alignment[strcount]=getStringFromFile(currentalignment,&curralignname);
    strcpy((*speciesnames)[strcount],curralignname);
    strcount++;
  }
  
  (*numspecies)=speciescount;
  return alignment;
  
}

char *getStringFromFile (FILE *datafile, char **title)
{
  char *first_str;
  unsigned int count;
  int k;  
  char c;                  /* c is used to read from file*/
  fpos_t startpos;         /* Used to keep track of the start position of
  			      various strings in the file */

  /* If the file doesn't exist return from the function */
  if (datafile == NULL)
    return NULL;
  
  
  
  
  
  
  
  /* Initialize c*/
  c = '\0';
  
  /* Scan through file until > character is reached.  The file is
     in FASTA format */
  while (c != '>')
  {
    if (c == EOF)
     return NULL;
   c = getc (datafile);
 }

 Fgetpos(datafile, &startpos);


  /* Calculate the length of the title, and read the title from the
     file after allocating space for it */
 c = '\0';
 count = 0;
 while ((c != '\n') && (c != EOF))
 {
  c = getc (datafile);
  count++;
}

*title = (char *) Calloc (count,sizeof (char ));

Fsetpos (datafile, &startpos);


c = '\0';
count = 0;
while ((c = getc(datafile)) != '\n') 
{
  if (c == EOF)
   break;

 (*title)[count] = c;
 count++;
}

  //(*title)[count] = '\0';

  //fprintf(stdout, "title in read function:%s. count %i\n", *title, count);

for(k=count-1;k>0;k--){
    //fprintf(stdout, "k %i char %c ", k, (*title)[k]);
  if((*title)[k]==' ')
    count--;
  else
    k=0;
}

(*title)[count] = '\0';

  //fprintf(stdout, "title in read function:%s. count %i\n", *title, count);
  /* Set start position of DNA string to be read*/
Fgetpos (datafile, &startpos);
c='\0';
count = 0;
  /* Calculate length of string, ignoring new-line characters
     since they intersperse the string in FASTA format */
while((c!='>') && (c!=EOF))
{
  if((c=='\n')||(c== ' '))
   c=getc(datafile);
 if(c!='>'){
   c = getc(datafile);
   count++;
 }
}

  /* Allocate space for the string */
first_str=(char *)Calloc(count,sizeof (char));

c = '\0';
count = 0;
Fsetpos (datafile, &startpos);

  /* Set the position of the file pointer back to the beginning
     of the string and read in the string ignoring new-line
     characters */
     while ((c = getc(datafile)) != '>' )
     {
      if ((c == '\n')||(c==' '))
       c=fgetc(datafile);
      /* Not sure if below commented region is needed when reading! */
      /*  if(c=='>') */
      /* 	{ */
      /* 	  fprintf(stdout,"have returned > char in while loop\n"); */
      /* 	  ungetc(c, datafile); */
      /* 	  break; */
      /* 	} */
     if ((c==EOF)||(c=='>'))
       break;

     first_str[count] = c;
     if(c!='\n')
       count++;
   }
  /* Terminate the string appropriately */
   first_str[count] = '\0';

   if(c=='>')
   {
    ungetc(c, datafile);
  }
  
  /* Convert uppercase characters to lowercase*/
  /*for(count=0;count<strlen(temp);count++)
    {
    if (temp[count] == 'A')
    temp[count]='a';
    else if (temp[count] == 'C')
    temp[count]='c';
    else if (temp[count] == 'T')
    temp[count]='t';
    else if (temp[count] == 'G')
    temp[count]='g';
    }
  */
  /* return temp*/
    fprintf(stdout, "done getting string\n");
    return first_str;
  }


  void printNexusFile(const Array2D<double> &solutions,const Array2D<int> &indexconversion, char **alltaxa, int numtaxa,const Array1D<int> &removedflags,double coeff, char *dirname){

    int i,j;
    FILE *outfile;
    int index;
    int unknowns=0;
    int knowns=0;
    FILE *missingdistcount;
    std::string dirName = std::string(dirname);
    std::string distanceFile = "distances.nexus";
    std::string distanceShort = "distances.short";
    std::string distanceCount = "distance.counts";
    if(dirname != NULL){
      distanceFile = dirName + "/" + distanceFile;
      distanceShort = dirName + "/" + distanceShort;
      distanceCount = dirName + "/" + distanceCount;
    }

    outfile=fopen(distanceFile.c_str(), "w");

    fprintf(stdout, "Distances\n");
    fprintf(outfile, "#NEXUS\n\nBEGIN taxa;\n\tDIMENSIONS ntax=%i;\nTAXLABELS\n", numtaxa);
    for(i=0;i<numtaxa;i++){
      fprintf(outfile, "[%i]\t%s\n", (i+1), alltaxa[i]);
    }
    fprintf(outfile, ";\nEND;\n\nBEGIN distances;\n\tDIMENSIONS ntaxa=%i;\n\tFORMAT\n\t\ttriangle=UPPER\n\t\tdiagonal\n\t\tlabels\n\t;\n\n\tMATRIX\n",numtaxa);
    for(i=0;i<numtaxa;i++){
      fprintf(outfile, "\t%s\t\t 0 ", alltaxa[i]);
      for(j=i+1;j<numtaxa;j++){
        index=indexconversion[i][j];

        if((solutions[index][0]!=0.0)&&(removedflags[index]==0)){
         fprintf(outfile, "%lf ", coeff*solutions[index][0]);
         fprintf(stdout, "name %s name %s %lf\n", alltaxa[i],alltaxa[j],coeff*solutions[index][0]);
         knowns++;
       }
       else{
         fprintf(outfile, "* ");
         fprintf(stdout, "*\n");
         unknowns++;
       }
     }
     fprintf(outfile, "\n");
   }
   fprintf(outfile, "\t;\nEND;\n");
   fclose(outfile);
   outfile=fopen(distanceShort.c_str(),"w");
   fprintf(outfile, "%i\n",numtaxa);
   for(i=0;i<numtaxa;i++){
    fprintf(outfile,"%s ", alltaxa[i]);
    for(j=0;j<numtaxa;j++){
      if(j>i)
       index=indexconversion[i][j];
     else
       index=indexconversion[j][i];
     if(i==j)
       fprintf(outfile, "0.0 ");
     else
       fprintf(outfile, "%lf ", coeff*solutions[index][0]);

   }
   fprintf(outfile,"\n");
 }
 fclose(outfile);
 missingdistcount=fopen(distanceCount.c_str(), "w");
 fprintf(missingdistcount, "Number of missing distances %i\n", unknowns);
 fprintf(missingdistcount, "Number of estimated distances %i\n",knowns);
 fclose(missingdistcount);
 return;
}
void printRatesFile(const Array2D<double> &rates,distanceList *distlist, char *dirname){

  int i;
  FILE *outfile;
  int name;
  int dim1=rates.dim1();
  double normalizingfactor = 0;
  double frac;
  distanceList *curr;
  std::string dirName = std::string(dirname);

  std::string treerates = "tree.rates";
    if(dirname != NULL){
      treerates = dirName + "/" + treerates;
    } 

  outfile=fopen(treerates.c_str(), "w");
  // outfile2=fopen("tree.rates.normalized","w");

  fprintf(stdout, "Rates\n");

  for(i=0;i<dim1;i++){
    normalizingfactor=normalizingfactor+(1/rates[i][0]);
  }
  //fprintf(stderr, "normalizing factor %lf, x0 %lf\n", normalizingfactor, rates[0][0]);
  name=1;
  i=0;
  curr=distlist;
  while(curr!=NULL){
    frac=dim1*((1/rates[i][0]))/normalizingfactor;
    fprintf(stdout, "'%s' %lf\n", getGeneNameFromDistanceList(curr),frac);
    fprintf(outfile, "%s %lf\n", getGeneNameFromDistanceList(curr),frac);
    curr=curr->next;
    i++;
    name++;
  }
  
  fprintf(stdout, "\n\n\n");
  
  
  return;
}
void printMatrixToFile(const Array2D<double> &coeffs_subset){
  FILE *outfile;
  int i,j, dim1, dim2;


  outfile=fopen("A.mat", "w");
  dim1=coeffs_subset.dim1();
  dim2=coeffs_subset.dim2();
  fprintf(stdout, "Output matrix Rows: %i, Cols: %i\n", dim1, dim2);
  for(i=0;i<dim1;i++){
    for(j=0;j<dim2;j++){
      fprintf(outfile, "%lf\t", coeffs_subset[i][j]);
    }
    fprintf(outfile, "\n");
  }
  return;
}


void printVectorToFile(const Array1D<double> &constants_subset){
  FILE *outfile;
  int i,dim1;


  outfile=fopen("B.vec", "w");
  dim1=constants_subset.dim1();
  fprintf(stdout, "Output matrix Rows: %i\n", dim1);
  for(i=0;i<dim1;i++){
    fprintf(outfile, "%lf\n", constants_subset[i]);
  }
  
  return;

}


void printMatrixtoMatlabformat(const Array2D<double> & currmat,FILE *outfile){
  int dim1=currmat.dim1();
  int dim2=currmat.dim2();
  int i,j;

  for(i=0;i<dim1;i++){
    for(j=0;j<dim2;j++){
      fprintf(outfile,"%lf\t",currmat[i][j]);
    }
    fprintf(outfile,"\n");
  }
  
  return;
}





#endif
