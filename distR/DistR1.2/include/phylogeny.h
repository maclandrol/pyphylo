/* Copyright November 2005: Rachel Bevan
   Email: rachel@mcb.mcgill.ca
   NOTE: NOT Done commenting
*/

#ifndef phylogeny_h
#define phylogeny_h

#include<stdio.h>
#include<stddef.h>
#include<stdlib.h>

typedef struct mynode{

  struct mynode *parent;     /* Parent of current node */
  struct mynode *leftsib;    /* Left sibling of current node */
  struct mynode *rightsib;   /* Right sibling of current node */
  struct mynode *leftchild;  /* Left child of current node */
  struct mynode *rightchild; /* Right child of current node */
  int id;             /* Corresonds to an array of names
			 of the taxa - array indexed by
			 this id value */
  double length;      /*Length of the branch leading
			from the parent to the node */
  /* Here can include some struct that contains info
     dependent upon the algorithm being run */
  char *genename;
  char **alignment;
  char **speciesnames;  /* Corresponding to species in alignment
			   index of species won't necessarily match
			   ID */
  int numspecies;
  int isroot;

} Treenode;

extern Treenode * newTreenode();

extern Treenode * nextTreenodePreOrder(Treenode * t);
extern Treenode * nextTreenodePostOrder(Treenode *t);
extern Treenode * getLeftMostLeafFromTreenode(Treenode *t);
extern Treenode * getRightMostLeafFromTreenode(Treenode *t);
extern int insertTreenodeAsLeftChild(Treenode *t,Treenode *child);
extern int insertTreenodeAsRightSibling(Treenode *t, Treenode *child, Treenode *leftsib);
extern int insertTreenode(Treenode *t, Treenode *child);
extern Treenode * removeTreenodeFromTree(Treenode *t);  
extern void removeTree(Treenode *root);
extern Treenode * cloneTree(Treenode* root);
extern Treenode * readPhylipTree(FILE *newickformat, char ***taxa_names, int *numspecies);
extern Treenode * makeTreeUnrooted(Treenode *root);
extern void writePhylipTree(FILE *outfilenewick, char **taxa_names, Treenode *root); 
extern int calculateNumberOfChildren(Treenode *root);
extern Treenode * findMostRecentCommonAncestor(Treenode *t1, Treenode *t2);
extern Treenode * fixSubtree(Treenode *root);
extern double calculatePairwiseDistance(Treenode *t1, Treenode *t2);
extern int printTreePostOrder(Treenode *root, char **taxa_names, int numspecies);
extern int printTreePreOrder(Treenode *root, char **taxa_names, int numspecies);
extern Treenode * getParentFromTreenode(Treenode *t);
extern int setParentInTreenode(Treenode *t, Treenode *p);
extern Treenode * getLeftsibFromTreenode(Treenode *t);
extern int setLeftsibInTreenode(Treenode *t, Treenode *lsib);
extern Treenode * getRightsibFromTreenode(Treenode *t);
extern int setRightsibInTreenode(Treenode *t, Treenode *rsib);
extern Treenode * getLeftchildFromTreenode(Treenode *t);
extern int setLeftchildInTreenode(Treenode *t, Treenode *lchild);
extern Treenode * getRightchildFromTreenode(Treenode *t);
extern int setRightchildInTreenode(Treenode *t, Treenode *rchild);
extern int getIDFromTreenode(Treenode *t);
extern int setIDInTreenode(Treenode *t, int ID);
extern double getLengthFromTreenode(Treenode *t);
extern int setLengthInTreenode(Treenode *t, double length);
extern int isTreenodeLeaf(Treenode *t);
extern int setGeneNameInRoot(Treenode *root, char *genename);
extern char * getGeneNameFromRoot(Treenode *root);
extern int setAlignmentInRoot(Treenode *root, char **alignment);
extern char ** getAlignmentFromRoot(Treenode *root);
extern int setSpeciesNamesInRoot(Treenode *root, char **speciesnames);
extern char ** getSpeciesNamesFromRoot(Treenode *root);
extern int setNumSpeciesInRoot(Treenode *root, int numspecies);
extern int getNumSpeciesFromRoot(Treenode *root);
extern int setRootInTreenode(Treenode *root);
extern int getRootFromTreenode(Treenode *root);
#endif
