/* Copyright November 2005: Rachel Bevan
   Email: rachel@mcb.mcgill.ca
   NOTE: NOT Done commenting
*/



#ifndef distR_cpp
#define distR_cpp

#include"distR.h"
#include"matrixOperations.h"
#include"fileio.h"



/* Need to fix so no alingment needed */

   double LeastSquaresError(const Array2D<double> &X, const Array2D<double> &truedistances, distanceList *mylist,char **alltaxa, int numtaxa, int numtrees){

    double sol=0.0;
    int i;
    distanceList *curr;
    Array2D<double> currarr;
    int id1, id2, id1corr, id2corr,index1,index2;
    char ** treespeciesnames;
    int currnumspecies;
    double distance, variance;
    double diff;
    double rate;
    double **distances;
    double **variances;


    curr=mylist;
    i=0;
    while(curr!=NULL){

      currnumspecies=getNumSpeciesFromDistanceList(curr);

      treespeciesnames=getSpeciesNamesFromDistanceList(curr);


      distances=getDistancesFromDistanceList(curr);
      variances=getVariancesFromDistanceList(curr);


      for(id1=0;id1<currnumspecies;id1++){
        id1corr=CorrectID(treespeciesnames[id1], alltaxa, numtaxa);

        for(id2=id1+1;id2<currnumspecies;id2++){

         id2corr=CorrectID(treespeciesnames[id2], alltaxa, numtaxa);


         if(id1corr<id2corr){
           index1=id1corr;
           index2=id2corr;
         }
         else{
           index1=id2corr;
           index2=id1corr;
         }

         distance=distances[id1][id2];
	//fprintf(stderr, "done with distances %lf\n", distance);
	//fprintf(stderr, "seqindex1 %i seqindex2 %i\n",seqindex1, seqindex2);
         variance=variances[id1][id2];

         if(distance>0){
	  //if(i==1)
	  // rate=1.0;
	  //else
           rate=X[i][0];
           diff=(truedistances[index1][index2]-(distance*rate));
           sol=sol+((diff*diff)/variance);
         }
	//fprintf(stderr, "sol %lf\n", sol);


       }
     }

     i++;
     curr=curr->next;
   }

  //fprintf(stdout, "sol %lf\n", sol);
   return sol;
 }


 int compareSpeciesLists(char **list1, int length1, char**list2, int length2){

  int i,j,found=0;


  for(i=0;i<length1;i++){
    for(j=0;j<length2;j++){
      //fprintf(stderr, "string i %s string j %s\n",list1[i],list2[j]);
      if(strcmp(list1[i],list2[j])==0){
	//fprintf(stderr,"found\n");
       found++;
     }
   }


 }
 
 if(length1!=length2)
  return -1;
if(length1!=found)
  return -1;
else
  return 0;



}

double ** calculateDistanceMatrixfromTree(Treenode *mytree, int numspecies){
  Treenode *nextleaf, *firstleaf;
  int id1,id2;
  double ** pairwisedistances;
  int i,j;

  pairwisedistances=(double **)Calloc(numspecies,sizeof(double*));
  for(i=0;i<numspecies;i++){
    pairwisedistances[i]=(double *)Calloc(numspecies,sizeof(double));
    for(j=0;j<numspecies;j++)
      pairwisedistances[i][j]=0.0;
  }
  
  /* Here I want to calculate pairwise branch distances */
  /* This works by traversing the leaves in the tree postorder,
     and at every leaf, traversing to every leaf remaining in the
     postorder traversal and calculating the distance between the two.
     This way the resulting distance matrix is triangular (probably
     upper in every case), and since this matrix should be symmetric
     this isn't a problem 
  */

     firstleaf=getLeftMostLeafFromTreenode(mytree);

     while(firstleaf!=NULL){
      if(isTreenodeLeaf(firstleaf)){
        id1=getIDFromTreenode(firstleaf);

        nextleaf=nextTreenodePostOrder(firstleaf);
        while((nextleaf!=mytree)&&(nextleaf!=NULL)){

         while(!isTreenodeLeaf(nextleaf)){
           nextleaf=nextTreenodePostOrder(nextleaf);
         }
         if(nextleaf!=NULL){
	  /* Get ID from this leaf so distances can be put in matrix */
           id2=getIDFromTreenode(nextleaf);

	  /* Call function to calculate the pairwise distances */
           pairwisedistances[id1][id2]=calculatePairwiseDistance(firstleaf, nextleaf);
           pairwisedistances[id2][id1]=pairwisedistances[id1][id2];

	  /*Now have next leaf through post order traversal */
           nextleaf=nextTreenodePostOrder(nextleaf);
         }
       }
     }

    /* Get next leaf in tree */
     firstleaf=nextTreenodePostOrder(firstleaf);
   }


   return pairwisedistances;
 }


 double ** calculateVarianceMatrixfromAlignment(double **currdistances, char **myalign, char ** alignspecies, char **distspecies, int numspecies){

  int i,j;
  double **pairwisevars;
  double distance=0.0;
  int seqindex1, seqindex2;
  int strlen1;
  int strlen2;
  
  pairwisevars=(double **)Calloc(numspecies,sizeof(double*));
  for(i=0;i<numspecies;i++){
    pairwisevars[i]=(double *)Calloc(numspecies,sizeof(double));
    for(j=0;j<numspecies;j++)
      pairwisevars[i][j]=0.0;
  }

 //  for(i=0;i<numspecies;i++){
//     fprintf(stderr, "index %i sequence %s\n",i,myalign[i]);
//   }

  for(i=0;i<numspecies;i++){
    seqindex1=CalculateSequenceIndex(distspecies[i],alignspecies,numspecies);
    strlen1=strlen(myalign[seqindex1]);
    for(j=i+1;j<numspecies;j++){

      seqindex2=CalculateSequenceIndex(distspecies[j],alignspecies,numspecies);
      strlen2=strlen(myalign[seqindex2]);
      
      distance=currdistances[i][j];

      // fprintf(stderr, "seqindex1 %i seqindex2 %i\n",seqindex1, seqindex2);
      // fprintf(stderr, "Seq1 %s\nstrlen %i\nseq2 %s\nstrlen %i\n",myalign[seqindex1],strlen1,myalign[seqindex2],strlen2);
      pairwisevars[i][j]=bulmerVarianceCalculation(myalign[seqindex1], myalign[seqindex2],distance);
      pairwisevars[j][i]=pairwisevars[i][j];
    }
  }

  // fprintf(stderr, "time to return\n");
  return pairwisevars;

}



void calculateTreeRates(distanceList *distmats, int numtrees, char **alltaxa, int numtaxa, char *dirName){

  Array2D<double> Dcoeffs(numtrees, numtrees,0.0);
  Array2D<double> Cinvcoeffs((((numtaxa*numtaxa) - numtaxa)/2),1,0.0);
  Array2D<double> BTrans(numtrees, (((numtaxa*numtaxa) - numtaxa)/2), 0.0);
  Array2D<double> B((((numtaxa*numtaxa) - numtaxa)/2),numtrees,0.0);
  Array2D<double> Vconstants((((numtaxa*numtaxa) - numtaxa)/2),1, 0.0);
  Array2D<double> result1;
  Array2D<double> vecresult;
  Array2D<double> matresult;
  Array2D<double> X;
  Array2D<double> Y;
  Array1D<int> removedflags((((numtaxa*numtaxa) - numtaxa)/2),0);
  Array2D<double> distances(numtaxa,numtaxa,0.0);
  Array2D<double> VarCoVarMatrix((((numtaxa*numtaxa) - numtaxa)/2),(((numtaxa*numtaxa) - numtaxa)/2),0.0);
  vector<Array2D<double> > treedistances(numtrees);
  char ** treespeciesnames;
  int currnumspecies;
  int id1=-1, id2=-1, id1corr, id2corr;
  int i,j;
  int coeffsdim=0;
  int index1=0,index2=0;
  double totaldistance = 0.0, distance=0.0;
  int pindex=0;
  Array2D<int> indexconversions(numtaxa,numtaxa,0); 
  int counter=0;
  double variance=0.0;
  double distvarratio=0.0;
  Array2D<double> Uconstants(numtrees,1,0.0);
  Array2D<double> Dcoeffsmod(numtrees, numtrees-1,0.0);
  Array2D<double> BigMatrix(numtrees+(((numtaxa*numtaxa) - numtaxa)/2),numtrees-1+(((numtaxa*numtaxa) - numtaxa)/2),0.0);
  distanceList *curr;
  double sumvariance=0.0;
  double sumdistance=0.0;
  double sumdistance2=0.0;
  double **currdistances;
  double distancesum=0.0;
  double largestdist=0.0;
  double totaldist=0.0;
  int largestindex=0;
  Array2D<double> lambdacoeffs;
  Array2D<double> constraintvals((((numtaxa*numtaxa) - numtaxa)/2),1,0.0);
  Array2D<double> constraintvalstrans(1,(((numtaxa*numtaxa) - numtaxa)/2),0.0);
  double kconst1=0.0;
  double kconst2=0.0;
  double kconst=0.0;
  Array2D<double> ACinvV;
  Array2D<double> ACinvVVtransCinvB;
  Array2D<double> ACinvVVtrans;
  Array2D<double> ACinvVVtransCinv;
  double VtransV;
  Array2D<double> finalScoeffs;
  Array2D<double> V, Vtrans;
  Array2D<double> constraintvalsD((((numtaxa*numtaxa) - numtaxa)/2),1,0.0);
  Array2D<double> constraintvalsDtrans(1,(((numtaxa*numtaxa) - numtaxa)/2),0.0);
  double totaldistvarratio=0.0;
  coeffsdim=((numtaxa*numtaxa)-numtaxa)/2;
  double **currvariances;

  fprintf(stderr, "starting calculateTreeRates\n");
  counter = 0;
  for(i=0;i<numtaxa;i++){
    for(j=i+1;j<numtaxa;j++){
      indexconversions[i][j]=counter;
      counter++;
    }
  }

  curr=distmats;

  i=0;
  

  while(curr!=NULL){
    totaldistance=0.0;
    sumvariance=0.0;
    sumdistance=0.0;
    sumdistance2=0.0;
    distancesum=0.0;
    //fprintf(stderr, "new iteration %i\n",i);
    currnumspecies=getNumSpeciesFromDistanceList(curr);
    
    treespeciesnames=getSpeciesNamesFromDistanceList(curr);
    
    
    currdistances=getDistancesFromDistanceList(curr);
    currvariances=getVariancesFromDistanceList(curr);

    totaldist=getTotalDistanceFromDistanceList(curr);
    if(totaldist>largestdist){
      totaldist=largestdist;
      largestindex=i;
    }
    for(id1=0;id1<currnumspecies;id1++){
      id1corr=CorrectID(treespeciesnames[id1], alltaxa, numtaxa);
      for(id2=id1+1;id2<currnumspecies;id2++){

       id2corr=CorrectID(treespeciesnames[id2], alltaxa, numtaxa);

       if(id1corr<id2corr){
         index1=id1corr;
         index2=id2corr;
       }
       else{
         index1=id2corr;
         index2=id1corr;
       }

       distance=currdistances[id1][id2];
       distancesum=distance+distancesum;
       variance=currvariances[id1][id2];
       if(variance==0.0)
         variance=0.000000000000001;
       fprintf(stderr, "currdistance %lf currvar %lf\n",distance,variance);
       sumvariance=sumvariance+variance;

       sumdistance=sumdistance+distance;
       sumdistance2=sumdistance2+(distance*distance);
       distvarratio=distance/variance;

       totaldistvarratio=totaldistvarratio+distvarratio;

	/* The following two lines were used to make sure
	   that the inverse rates are being calculated properly */

       pindex=indexconversions[index1][index2];
       constraintvalsD[pindex][0]=constraintvalsD[pindex][0]+(1/variance);
       constraintvalsDtrans[0][pindex]=constraintvalsD[pindex][0];
       if(i==0){
         constraintvals[pindex][0]=-1/variance;
         constraintvalstrans[0][pindex]=1/variance;


       }

       totaldistance=totaldistance+(distance*distvarratio);

       BTrans[i][pindex]=-2*distvarratio;

       B[pindex][i]=-2*distvarratio;



	//here I am calculating C
       Cinvcoeffs[pindex][0]=Cinvcoeffs[pindex][0]+(2/variance);

     }
   }



   Dcoeffs[i][i]=2*totaldistance;
   if(i==0)
    kconst1=totaldistance;
  

  i++;

  curr=curr->next;
}
kconst2=totaldistvarratio;

for(i=0;i<Cinvcoeffs.dim1();i++){

  if(Cinvcoeffs[i][0]==0.0){
    removedflags[i]=1;
  }
}




  /* Calculate C-1 since this is what I need to solve for
     X and Y */
for(i=0;i<Cinvcoeffs.dim1();i++){
  if(Cinvcoeffs[i][0]!=0.0)
    Cinvcoeffs[i][0]=1/Cinvcoeffs[i][0];

}




V=constraintvalsD;
  //fprintf(stderr, "V\n");

Vtrans=constraintvalsDtrans;
kconst=kconst2;

  //multiply A*C-1
matrixMultiplyBTransDiag(result1,BTrans,Cinvcoeffs);
  //fprintf(stderr, "ACinv\n");

  //get AC-1B into matresult
matrixMultiply(matresult,result1,B);

  //not used
  //matrixMultiply(vecresult,result1,Vconstants);



  //multipy AC-1 with V
matrixMultiply(ACinvV,result1,V);


  //multiply result with Vtrans
matrixMultiply(ACinvVVtrans,ACinvV,Vtrans);

  //multiply result with C-1
matrixMultiplyBTransDiag(ACinvVVtransCinv,ACinvVVtrans,Cinvcoeffs);

  //multiply result with B
matrixMultiply(ACinvVVtransCinvB,ACinvVVtransCinv,B);
Array2D<double> VtransCinv;
  // multiply VT and C-1
matrixMultiplyBTransDiag(VtransCinv,Vtrans,Cinvcoeffs);
  // get scalar quatinity of prev result multiplied by V
vectorMultiply(&VtransV,VtransCinv,V);

  // matresult has Btran * Cinv *B

  // calculate constant vector
for(i=0;i<ACinvV.dim1();i++){
  for(j=0;j<ACinvV.dim2();j++){

    ACinvV[i][j]=-ACinvV[i][j]*kconst/VtransV;
    
  }
}
  //calculate other coeffects of S
for(i=0;i<ACinvVVtransCinvB.dim1();i++){
  for(j=0;j<ACinvVVtransCinvB.dim2();j++){

    ACinvVVtransCinvB[i][j]=ACinvVVtransCinvB[i][j]/VtransV;

  }
}

//   Array2D<double> matmod3(numtrees,numtrees,0.0);
//   /*settup of my solution here */
//   for(i=0;i<(matmod3.dim1());i++){
//     for(j=0;j<(matmod3.dim2());j++)
//       matmod3[i][j]=Dcoeffs[i][j];
//   }

  //largestindex=0;
 //  matmod3[0][0]=0.0;

//   Array2D<double> vecmod3(numtrees,1,0.0);
//   vecmod3[0][0]=-Dcoeffs[0][0];

  // get first set of S coeffs ready
matrixDifference(matresult,Dcoeffs,matresult);
  // add to second set of s coeffs
matrixAddition(finalScoeffs,matresult, ACinvVVtransCinvB);

QR<double> decomp(finalScoeffs);
  /* solve for rates */
if(decomp.isFullRank()){
  X=decomp.solve(ACinvV);
}
else
  fprintf(stdout, "not performing QR solve\n");


  /* solve for distances*/
double normalizingfactor=0.0;
for(i=0;i<X.dim1();i++){
  normalizingfactor=normalizingfactor+(1/X[i][0]);
}

int dim1 = X.dim1();

 //  fprintf(stdout, "Rates:\n");

//   for(i=0;i<X.dim1();i++){
//     fprintf(stdout, "%.15lf ", (1/X[i][0]));
//   }
//   fprintf(stdout, "\n");

matrixMultiply(result1,B,X);

Array2D<double> VtransCinvBS;

matrixMultiply(VtransCinvBS,VtransCinv,result1);



Array2D<double>VtransCinvBSV(V.dim1(),V.dim2(),0.0);

for(i=0;i<VtransCinvBSV.dim1();i++){
  for(j=0;j<VtransCinvBSV.dim2();j++){
    VtransCinvBSV[i][j]=VtransCinvBS[0][0]*V[i][j]/VtransV;

  }

}


Array2D<double> newvec(V.dim1(),V.dim2(),0.0);
for(i=0;i<V.dim1();i++){
  for(j=0;j<V.dim2();j++){
    newvec[i][j]=kconst*V[i][j]/VtransV;

  }


}
Array2D<double> result2;
matrixAddition(result2,newvec,VtransCinvBSV);

matrixDifference(result1,result2,result1);


matrixMultiplyDiagA(Y,Cinvcoeffs,result1);
 //  fprintf(stdout, "Distances:\n");
//    for(i=0;i<Y.dim1();i++){
//      for(j=0;j<Y.dim2();j++){
//        fprintf(stdout, "%.15lf ", Y[i][j]);
//      }
//      fprintf(stdout, "\n");
//    }
//    fprintf(stdout, "\n");






printRatesFile(X,distmats, dirName);

double coeff = normalizingfactor/X.dim1();

printNexusFile(Y,indexconversions,alltaxa,numtaxa,removedflags,coeff, dirName);

  //createDistanceMatrix(Y,indexconversions,numtaxa,removedflags,distances);
  /* use matrix distances to pass to leastsquareserror function */
  //double lse= LeastSquaresError(X,Y, distmats,alltaxa, numtaxa,numtrees);
  //fprintf(stderr, "lse %lf\n", lse);
return;
}






int CalculateSequenceIndex( char *taxaname, char **alignmentnames,int numspecies){

  int i;
  int found=-1;

  for(i=0;i<numspecies;i++){
    if(strcmp(taxaname,alignmentnames[i])==0)
      found=i;
  }
  return found;


};

double bulmerVarianceCalculation(char * seq1, char *seq2,double distance){

  int nogaps=0,i;
  int length=strlen(seq1);
  int length2=strlen(seq2);
  double val;
  double b=0.93;

  //fprintf(stderr, "strlen seq1 %i strlen seq2 %i\n",strlen(seq1),strlen(seq2));

  if(length2!=length){
    fprintf(stderr, "Problem with alignment.  Trying to calculate variance between 2 sequences of different length.  Exiting program...\n");
    exit(1);
  }
  else{
    for(i=0;i<length;i++){
      if((seq1[i]!='-')&&(seq2[i]!='-'))
       nogaps++;
   }
    //fprintf(stdout, "nogaps %i\n", nogaps);
    /* does val need to be negated or not? */
    // nogaps=2000;
   val=(exp(2*distance/b)*b*(1-exp(-distance/b))*(1-b+b*exp(-distance/b)))/nogaps;
    //fprintf(stdout, "val %lf distance %lf\n", val,distance);

   return val;
 }
}

int CorrectID( char * taxaname, char **alltaxa, int numtaxa){

 int i, index=-1;


   //fprintf(stdout, "taxaname %s\n",taxaname); 
 for(i=0;i<numtaxa;i++){
     // fprintf(stdout, "alltaxa[i] %s\n", alltaxa[i]);
   if(strcmp(taxaname, alltaxa[i])==0)
     index=i;
 }



 return index;

}

void createDistanceMatrix(const Array2D<double> &solutions,const Array2D<int> &indexconversions,int numtaxa,const Array1D<int> &removedflags, Array2D<double> &distances ){

  int i,j,index;

  for(i=0;i<numtaxa;i++){
    for(j=i+1;j<numtaxa;j++){
      index=indexconversions[i][j];
      if((solutions[index][0]!=0.0)&&(removedflags[index]==0))
       distances[i][j]=solutions[index][0];
     else
       distances[i][j]=INF;
   }
 }



 return;
};



#endif
