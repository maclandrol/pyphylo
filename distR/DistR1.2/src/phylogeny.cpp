/***********************************************
author: Rachel Bevan 
copyright: Rachel Bevan November 2005.

date (creation): 08/01/04
version: 1.0

The "bare bones" phylogeny. Apart from structure in tree we store
an int (corresponding to taxa#) for each leaf and a double
for the edge length.  Also taxa names, and alignment are stored
in the 'root' of the tree.

we have some basic utilities for this structure.

**************************************************/
#include "phylogeny.h"
#include "wrapper.h"
#include <ctype.h>

Treenode * newTreenode(){
  Treenode *t;
  t=(Treenode*)Calloc(1,sizeof(Treenode));
  t->parent = NULL;
  t->leftsib = NULL;
  t->rightsib = NULL;
  t->leftchild = NULL;
  t->rightchild = NULL;
  t->id = -1;
  t->length =-1;
  t->genename=NULL;
  t->alignment=NULL;
  t->speciesnames=NULL;
  t->numspecies=-1;
  t->isroot=0;
  return t;

}

Treenode * nextTreenodePreOrder(Treenode * t){

  Treenode *curr;
    
  if(t!=NULL){
    if (!isTreenodeLeaf(t))
      curr = getLeftchildFromTreenode(t);
    else {
      curr = t;
      while ((curr!=NULL) && (getRightsibFromTreenode(curr) == NULL))
	curr = getParentFromTreenode(curr);
      if (curr!=NULL)
	curr=getRightsibFromTreenode(curr);
    }
    return curr;
  }
  else
    return NULL;


}

/* NOTE: For this function to work it is imperative to send the leftmost
   node as the root for the first call! */
Treenode * nextTreenodePostOrder(Treenode *t){
  Treenode *curr;
  
  if(t!=NULL){
    if (getRightsibFromTreenode(t)!=NULL){
      curr = getLeftMostLeafFromTreenode(getRightsibFromTreenode(t));
    }
    else
      curr = getParentFromTreenode(t);
    
    return curr;
    
  }
  else
    return NULL;
}


Treenode * getLeftMostLeafFromTreenode(Treenode *t){
  Treenode *curr;
  Treenode *next;

  if(t!=NULL){
    curr = t;
    while((next=getLeftchildFromTreenode(curr))!=NULL)
      curr=next;
    return curr;
  }
  else
    return NULL;
}

Treenode * getRightMostLeafFromTreenode(Treenode *t){
  Treenode *curr;
  Treenode *next;

  if(t!=NULL){
    curr = t;
    while((next=getRightchildFromTreenode(curr))!=NULL)
      curr=next;
    return curr;
  }
  else
    return NULL;
}



int insertTreenodeAsLeftChild(Treenode *t,Treenode *child){
  Treenode *lchild,*rchild,*rsib,*rightmost;

  if((t!=NULL)&&(child!=NULL)){
    lchild = getLeftchildFromTreenode(t);
    rchild = getRightchildFromTreenode(t);
    setParentInTreenode(child, t);
    rightmost=child;
    while((rsib=getRightsibFromTreenode(child))!=NULL){
      setParentInTreenode(rsib,t);
      rightmost=rsib;
    }
    if(lchild!=NULL){
      setRightsibInTreenode(rightmost,lchild);
      setLeftsibInTreenode(lchild,rightmost);
    }
    setLeftchildInTreenode(t,child);
    if (rchild == NULL)
      setRightchildInTreenode(t,child);
    return 0;
  }
  else 
    return -1;
}

int insertTreenodeAsRightSibling(Treenode *t, Treenode *child, Treenode *leftsib){

  Treenode *rightsib;

  if((t!=NULL)&&(child!=NULL)&&(leftsib!=NULL)){
    rightsib=getRightsibFromTreenode(leftsib);
    setParentInTreenode(child, t);
    if(getParentFromTreenode(leftsib)!=t)
      fprintf(stderr, "Error: Trying to insert child next to a non-child");
    setLeftsibInTreenode(child,leftsib);
    
    setRightsibInTreenode(leftsib,child);
    
    if(rightsib!=NULL){
      setRightsibInTreenode(child,rightsib);
      setLeftsibInTreenode(rightsib,child);
    }
    else
      setRightchildInTreenode(t,child);
    return 0;
  }
  else
    return -1;
}

int insertTreenode(Treenode *t, Treenode *child){
  Treenode *leftsib;

  if(t!=NULL){
    if(isTreenodeLeaf(t))
      insertTreenodeAsLeftChild(t,child);
    else{
      leftsib = getRightchildFromTreenode(t);
      insertTreenodeAsRightSibling(t, child, leftsib);
    }
    
    return 0;
  }
  else
    return -1;
}


/* Returns parent of deleted node */

Treenode * removeTreenodeFromTree(Treenode *t){
  Treenode *leftsib, *rightsib, *parent, *Prightchild, *Pleftchild, *child1,*child2;

  if(t!=NULL){
    leftsib=getLeftsibFromTreenode(t);
    rightsib=getRightsibFromTreenode(t);
    parent=getParentFromTreenode(t);
    
    if(leftsib!=NULL)
      setRightsibInTreenode(leftsib,rightsib);
    if(rightsib!=NULL){
      setLeftsibInTreenode(rightsib,leftsib);
    }
    if(parent!=NULL){
      Prightchild=getRightchildFromTreenode(parent);
      Pleftchild=getLeftchildFromTreenode(parent);
      if(Pleftchild==t)
	setLeftchildInTreenode(parent,rightsib);
      else if(Prightchild==t)
	setRightchildInTreenode(parent,leftsib);
    }
    
   
    
    child1=getLeftchildFromTreenode(t);
    
    if(leftsib!=NULL){
      child2=getRightsibFromTreenode(child1); /* Should work even if child1 NULL
						 since code checks for NULL pointers */
      while(child1!=NULL){
	insertTreenodeAsRightSibling(parent, child1, leftsib); /*D has child1 and leftsib
								 reversed in his code but I 
								 think that it is wrong */
	leftsib=child1;
	child1=child2;
	child2=getRightsibFromTreenode(child1);
	
	
      }
    }
    else if(leftsib==NULL){
      if(child1!=NULL)
	insertTreenodeAsLeftChild(parent,child1);
    }
   
    free(t);
    t=NULL;
  }
  return parent;
};  

void removeTree(Treenode *root){
  Treenode *leftmostchild, *curr,*next;
  
  if(root!=NULL){
    leftmostchild = getLeftMostLeafFromTreenode(root);
    curr = leftmostchild;
    while(curr!=root){
      /* Need to check nextTreenodePostOrder */
      //fprintf(stdout, "Branch Length %lf, node ID %i\n", getLengthFromTreenode(curr),curr->id);
      next = nextTreenodePostOrder(curr);
      Free(curr,sizeof(Treenode));
      curr = next;
      /* Dangling pointers?*/
    }
    Free(root,sizeof(Treenode));
    root=NULL;
  }
  return;
  
}

/********************
clone_tree

Allocates memory and returns the root of a copy of T
*********************/

Treenode * cloneTree(Treenode * root) {
    Treenode* newT = newTreenode();
    Treenode* p;
    Treenode*child;
    
    

    newT->length = root->length;
    newT->id = root->id;
    newT->numspecies=root->numspecies;
    newT->speciesnames=root->speciesnames;
    newT->genename=root->genename;

    p=getLeftchildFromTreenode(root);
    while(p!=NULL){
      child=cloneTree(p);
      insertTreenode(newT, child);
      p=getRightsibFromTreenode(p);
      
    }
    return newT;
}


Treenode * readPhylipTree(FILE *newickformat, char ***taxa_names, int *numspecies){
  
  char next_tok,c,prevc;
  int id,namelength,currchar,counter,size=0;
  Treenode *treeptr, *root=NULL, *currptr;
  double weight;
  fpos_t start;
  int prevstate=0;
  int state = 0;  /*state 1 after (
		    state 2 after ,
		    state 3 after )
		    state 4 while reading in data */
  
  
  /* count the number of species */
  

  Fgetpos(newickformat,&start);
  
  (*numspecies)=0;
  prevc='a';
  while(((c=fgetc(newickformat))!=';')&&(c!=EOF)){
    //if((c==':')&&(prevc!=')'))
    if(c==',')
      (*numspecies)++;
    prevc=c;
  }
  
  if(*numspecies==0){
    fprintf(stderr, "ERROR:  Problem with newick format.  Please check input file.  Exiting program...\n");
    exit(1);
  }
  //fprintf(stderr, "commacount %i,numspecies %i\n",commacount, *numspecies);
  (*numspecies)++;
  
  *taxa_names = (char**)Calloc((*numspecies),sizeof(char *));
  id=0;
  Fsetpos(newickformat,&start);
  
  
  
  while((next_tok=fgetc(newickformat))!=';'){
   
    counter = 0;
    if(next_tok=='('){
     
      if(root==NULL){
	
	root = newTreenode();
	currptr=root;
	size++;
      }
      else{
	
	size++;
	prevstate=state;
	treeptr=newTreenode();
	insertTreenode(currptr,treeptr);
	currptr=getRightchildFromTreenode(currptr);

	if(currptr==NULL)
	  fprintf(stderr, "Error - child not properly inserted into tree\n");
      }
      state =1;
      
    }
   
    else if(next_tok==','){
      
      prevstate=state;
      state=2;
      /* Not sure if should move up level here*/
    }
    else if(next_tok==')'){
      
      if(root == NULL)
	fprintf(stderr, "Error - treefile has incorrect format!\n");
      else{
	prevstate=state;
	if((prevstate==4)||(prevstate==3)){
	  next_tok=fgetc(newickformat);
	  while(isdigit(next_tok))
	    next_tok=fgetc(newickformat);
	  if(next_tok==':'){
	    weight=0.0;
	    fscanf(newickformat,"%lf",&weight);
	    setLengthInTreenode(currptr,weight);
	    
	  }
	  else
	    ungetc(next_tok,newickformat);
	  currptr=getParentFromTreenode(currptr);
	}
	else{
	  fprintf(stderr, "in code didn't think would use check!\n");
	  size++;
	  treeptr=newTreenode();
	  insertTreenode(currptr,treeptr);
	  currptr=getRightchildFromTreenode(currptr);
	}
      }
      state=3;
      
    }
    else if((next_tok=='\n')||(next_tok==' '))
      continue;
    
    else{
      
      prevstate=state;
      state=4;
      ungetc(next_tok,newickformat);
      Fgetpos(newickformat,&start);
      namelength=0;
      while(((next_tok=fgetc(newickformat))!=':')&&(next_tok!=',')&&(next_tok!=')')){
	namelength++;
      }
      namelength++;
      (*taxa_names)[id]=(char*)Calloc(namelength,sizeof(char));
      
      Fsetpos(newickformat,&start);
      currchar=0;
      while(((next_tok=fgetc(newickformat))!=':')&&(next_tok!=',')&&(next_tok!=')')){
	(*taxa_names)[id][currchar]=next_tok;
	//fprintf(stderr,"%c\n",next_tok);
	currchar++;
      }
      (*taxa_names)[id][currchar]='\0';
      if((next_tok==',')||(next_tok==')'))
	ungetc(next_tok,newickformat);
      weight=0.0;
      if(next_tok==':'){
	
	fscanf(newickformat,"%lf",&weight);
      }
      
      
      if((prevstate==1)||(prevstate==2)){
	size++;
	treeptr=newTreenode();
	setIDInTreenode(treeptr,id);
	setLengthInTreenode(treeptr,weight);
	
	insertTreenode(currptr,treeptr);
      }
      else if(prevstate==3){
	fprintf(stderr, "Error - with treefile, cannot create leaf after a leaf\n");
	currptr=getParentFromTreenode(currptr); /*not sure if this is necessary
						  should check */
      }
      id++;
      
    }
    
  }
   
  root->isroot=1;
  return root;
  

}

Treenode * makeTreeUnrooted(Treenode *root){
  int children;
  Treenode *newrightchild;
  Treenode *leftchild;
  double totallen;
  
  
  
  if(root!=NULL){
    
    if(root->isroot==1){
      children=calculateNumberOfChildren(root);
      
      if(children==2){
	//fprintf(stdout, "pointers - root %i root->lchild %i root->rchild %i root->lchild->parent %i root->lchild->lsibling %i root->lchild->rsibling %i root->rchild->parent %i root->rchild->lsibling %i\n,root->rchild->rsibling %i\n", root, root->leftchild, root->rightchild, root->leftchild->parent, root->leftchild->leftsib, root->leftchild->rightsib,root->rightchild->parent,root->rightchild->leftsib,root->rightchild->rightsib);
	//printTreePostOrder(root,taxa_names,root->numspecies);
	//fprintf(stderr, "have 2 children\n");
	
	
	newrightchild=getRightchildFromTreenode(root);
	leftchild=getLeftchildFromTreenode(root);
	//fprintf(stderr, "leftchild %i rightchild %i\n",leftchild, newrightchild);
	
	
	
	totallen=getLengthFromTreenode(newrightchild)+getLengthFromTreenode(leftchild);
	if(!isTreenodeLeaf(leftchild)){
	  setLengthInTreenode(newrightchild,totallen);
	  insertTreenode(leftchild,newrightchild);
	  newrightchild->leftsib=NULL;
	  
	  leftchild->parent=NULL;
	  leftchild->rightsib=NULL;
	  leftchild->isroot=1;
	  leftchild->length=0.0;
	  leftchild->numspecies=root->numspecies;
	  free(root);
	  root=leftchild;
	}
	else{
	  setLengthInTreenode(leftchild,totallen);
	  
	  leftchild->rightsib=NULL;
	  
	  leftchild->leftsib=NULL;
	  //fprintf(stderr, "insert leftchild\n");
	  insertTreenodeAsLeftChild(newrightchild,leftchild);
	  //fprintf(stderr, "insert worked\n");
	  
	  
	  
	  
	  
	  newrightchild->parent=NULL;
	  
	  newrightchild->leftsib=NULL;
	  
	  newrightchild->rightsib=NULL;
	  
	  newrightchild->isroot=1;
	  
	  newrightchild->length=0.0;
	  newrightchild->numspecies=root->numspecies;
	  free(root);
	  root=newrightchild;
	  //fprintf(stderr, "leave else\n");
	}
	//fprintf(stdout, "pointers - root %i root->lchild %i root->rchild %i root->lchild->parent %i root->lchild->lsibling %i root->lchild->rsibling %i root->rchild->parent %i root->rchild->lsibling %i\n,root->rchild->rsibling %i root->parent %i\n", root, root->leftchild, root->rightchild, root->leftchild->parent, root->leftchild->leftsib, root->leftchild->rightsib,root->rightchild->parent,root->rightchild->leftsib,root->rightchild->rightsib,root->parent);
	//printTree(leftchild,taxa_names,leftchild->numspecies);
	//fprintf(stderr, "done printtree postorder\n");
      }
      //fprintf(stderr, "done\n");
      return root;
    }
    else
      return NULL;
  }
  else
    return NULL;
  
} 

void writePhylipTree(FILE *outfilenewick, char **taxa_names, Treenode *root){

  if(root!=NULL){
    if(isTreenodeLeaf(root)) {
      fprintf(outfilenewick, "%s", taxa_names[getIDFromTreenode(root)]);
    }
    else {
      fprintf(outfilenewick,"(");
      Treenode *p=getLeftchildFromTreenode(root);
      while(p!=NULL){
      	if (p!=getLeftchildFromTreenode(root)) 
	  fprintf(outfilenewick, ",");
	writePhylipTree(outfilenewick,taxa_names,p);
	p=getRightsibFromTreenode(p);
      }
      fprintf(outfilenewick,")");
    }
    if ((root->isroot==0)) //comment if you don't want to print lengths 
      fprintf(outfilenewick, ":%.15lf",getLengthFromTreenode(root));
    else
      fprintf(outfilenewick, ";");
  }
  else
    fprintf(stderr, "NULL NODE!\n");
}

int calculateNumberOfChildren(Treenode *root){
  int count=0;
  Treenode *child;
  
  if(root!=NULL){
    child=getLeftchildFromTreenode(root);
    while(child!=NULL){
      count++;
      child=getRightsibFromTreenode(child);
    }
    return count;
  }
  else
    return -1;

}


Treenode * findMostRecentCommonAncestor(Treenode *t1, Treenode *t2){

  Treenode *searchnode,*commonparent;
  
  if((t1!=NULL)&&(t2!=NULL)){

    searchnode=t2;
    commonparent=getParentFromTreenode(t1);
    while((searchnode!=commonparent)){
      if(searchnode!=NULL){
	searchnode=getParentFromTreenode(searchnode);
      }
      else
	{
	  searchnode=t2;
	  commonparent=getParentFromTreenode(commonparent);
	  if(commonparent==NULL){
	    fprintf(stderr, "ERROR:  No common ancestor for these two leaves, problem with tree. Exiting program...\n");
	    exit(1);
	  }
	}
      
    }
    return commonparent;
  }

  else
    return NULL;

}

Treenode * fixSubtree(Treenode *root){

  Treenode *leftchild;
  Treenode *rightchild;
  Treenode *curr;
  Treenode *retnode;
  double totallen;
  int next = 0;

  
  if(root!=NULL){
    curr=root;
    while(curr!=NULL){
      leftchild=getLeftchildFromTreenode(curr);
      rightchild=getRightchildFromTreenode(curr);
      if((leftchild!=NULL)&&(rightchild!=NULL)){
	
	if(leftchild==rightchild){
	
	  if(curr!=root){
	
	    /* need to check if curr is root since root length is I think
	       -1, but don't want to actually add these values */
	
	    totallen=getLengthFromTreenode(curr)+getLengthFromTreenode(leftchild);
	
	    setLengthInTreenode(leftchild,totallen);
	    
	    curr=removeTreenodeFromTree(curr);
	
	    
	  }
	  else{
	
	    setRootInTreenode(leftchild);
	    setParentInTreenode(leftchild,NULL);
	    free(curr);
	    curr=leftchild;
	    next = 1;
	  }
	}
      }
      if(next==0){
	
	curr=nextTreenodePreOrder(curr);
      }
      else
	next = 0;
    }
    retnode=makeTreeUnrooted(root);
    return root;
  }
  else
    return NULL;
  
  
}


double calculatePairwiseDistance(Treenode *t1, Treenode *t2){

  Treenode *commonparent;
  double distance1=0.0, distance2 = 0.0;

  if((t1!=NULL)&&(t2!=NULL)){
    commonparent = findMostRecentCommonAncestor(t1, t2);
    while(t1!=commonparent){
      distance1=getLengthFromTreenode(t1)+distance1;
      t1=getParentFromTreenode(t1);
    }
    while(t2!=commonparent){
      distance2=getLengthFromTreenode(t2)+distance2;
      t2=getParentFromTreenode(t2);
    }
    return(distance1+distance2);
  }
  else
    return -1.0;

}


int printTreePostOrder(Treenode *root, char **taxa_names, int numspecies){

  Treenode *curr;
  int id;

  
  if(root!=NULL){
    fprintf(stdout, "POSTORDER traversal of tree\n");
    curr=getLeftMostLeafFromTreenode(root);
    if(numspecies!=0)
      fprintf(stdout, "Number of species %i\n", numspecies);
    while(curr!=NULL){
      /* Uncomment next line if you wish to look at the pointer values
	 as well as the data in the tree */
      /*fprintf(stderr, "currptr %i\n", curr);*/
      id=getIDFromTreenode(curr);
      fprintf(stdout, "Branch Length %lf, ID %i", getLengthFromTreenode(curr),id);
      if(id!=-1)
	fprintf(stdout, " %s", taxa_names[id]);
      fprintf(stdout, "\n");
      curr=nextTreenodePostOrder(curr);
    }
    return 0;
  }
  else
    return -1;
  
  
}

int printTreePreOrder(Treenode *root, char **taxa_names, int numspecies){

  Treenode *curr;
  int id;

  if(root!=NULL){
    curr=root;
    fprintf(stdout, "PREORDER traversal of tree\n");
    if(numspecies!=0)
      fprintf(stdout, "Number of species %i\n", numspecies);
    while(curr!=NULL){
      id=getIDFromTreenode(curr);
      fprintf(stdout, "Branch Length %lf, ID %i", getLengthFromTreenode(curr),id);
      if(id!=-1)
	fprintf(stdout, " %s", taxa_names[id]);
      fprintf(stdout, "\n");

      curr=nextTreenodePreOrder(curr);
    }

    return 0;
  }
  else
    return -1;
}

Treenode * getParentFromTreenode(Treenode *t){

  if(t!=NULL)
    return t->parent;
  else
    return NULL;


}
int setParentInTreenode(Treenode *t, Treenode *p){

  if(t!=NULL){
    t->parent=p;
    return 0;
  }
  else
    return -1;


}
Treenode * getLeftsibFromTreenode(Treenode *t){

  if(t!=NULL)
    return t->leftsib;
  else
    return NULL;


}
int setLeftsibInTreenode(Treenode *t, Treenode *lsib){

  if(t!=NULL){
    t->leftsib=lsib;
    return 0;
  }
  else
    return -1;

}
Treenode * getRightsibFromTreenode(Treenode *t){
  
  if(t!=NULL)
     return t->rightsib;
  else
    return NULL;

}
int setRightsibInTreenode(Treenode *t, Treenode *rsib){

  if(t!=NULL){
    t->rightsib=rsib;
    return 0;
  }
  else
    return -1;

}
Treenode * getLeftchildFromTreenode(Treenode *t){
  
  if(t!=NULL)
     return t->leftchild;
  else
    return NULL;

}
int setLeftchildInTreenode(Treenode *t, Treenode *lchild){

  if(t!=NULL){
    t->leftchild=lchild;
    return 0;
  }
  else
    return -1;

}
Treenode * getRightchildFromTreenode(Treenode *t){
  
  if(t!=NULL)
     return t->rightchild;
  else
    return NULL;

}
int setRightchildInTreenode(Treenode *t, Treenode *rchild){

  if(t!=NULL){
    t->rightchild=rchild;
    return 0;
  }
  else
    return -1;

}
int getIDFromTreenode(Treenode *t){
  
  if(t!=NULL)
     return t->id;
  else
    return -1;

}
int setIDInTreenode(Treenode *t, int ID){

  if(t!=NULL){
    t->id=ID;
    return 0;
  }
  else
    return -1;

}
double getLengthFromTreenode(Treenode *t){
  
  if(t!=NULL)
     return t->length;
  else
    return -1;

}
int setLengthInTreenode(Treenode *t, double length){

  if(t!=NULL){
    t->length=length;
    //fprintf(stdout, "LENGTH SET HERE %lf\n", t->length);
    return 0;
  }
  else
    return -1;

}
int isTreenodeLeaf(Treenode *t){
  if(t!=NULL){
    if(getLeftchildFromTreenode(t)!=NULL)
      return 0;
    else
      return 1;
  }
  else
    return -1;

}

int setGeneNameInRoot(Treenode *root, char *genename){

  if(root!=NULL){
    root->genename=genename;
    return 0;
  }
  else
    return -1;

}
char * getGeneNameFromRoot(Treenode *root){

  if(root!=NULL)
    return root->genename;
  else
    return NULL;

};

int setAlignmentInRoot(Treenode *root, char **alignment){

  if(root!=NULL){
    root->alignment=alignment;
    return 0;
  }
  else
    return -1;

};
char ** getAlignmentFromRoot(Treenode *root){

  if(root!=NULL)
    return root->alignment;
  else
    return NULL;

};

int setSpeciesNamesInRoot(Treenode *root, char **speciesnames){

  if(root!=NULL){
    root->speciesnames=speciesnames;
    return 0;
  }
  else
    return -1;

};
char ** getSpeciesNamesFromRoot(Treenode *root){

  if(root!=NULL)
    return root->speciesnames;
  else
    return NULL;

};

extern int setNumSpeciesInRoot(Treenode *root, int numspecies){

  if(root!=NULL){
    root->numspecies=numspecies;
    return 0;
  }
  else
    return -1;
};
extern int getNumSpeciesFromRoot(Treenode *root){

  if(root!=NULL){
    return root->numspecies;
  }
  else
    return -1;
};

extern int setRootInTreenode(Treenode *root){

  if(root!=NULL){
    root->isroot=1;
    return 0;
  }
  else
    return -1;

}

extern int getRootFromTreenode(Treenode *root){

  if(root!=NULL){
    return root->isroot;
  }
  else
    return -1;
}


