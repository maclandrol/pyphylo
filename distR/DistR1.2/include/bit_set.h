/*******************************************************
filename: 	bit_set.h
author: 	David Bryant (bryant@mcb.mcgill.ca)
created: 	June 2001
copyright:	David Bryant 2001

A simple and efficient implementation of a set of non-negative integers. The
set can be a subset of 0,1,...,max_cap-1. The size of max_cap can
be extended using the reserve command.

For efficiency, there is no dynamic resizing of the set, and all comparisons between sets assume 
that the two sets have the same capacity max_cap. To perform runtime checking of this, uncomment
the following line:
*******************/
//#define BITSET_BOUNDS
/*****************

IMPLEMENTATION DETAILS
The elements are stored in a vector of words. Each word is an unsigned long, and
each bit in each word corresponds to a different element.

The smallest elements correspond to the rightmost bits in the first word.
In the final word the bits corresponding to elements exceeding max_cap are 
kept false - this is important for many set comparisons.

Some of the ideas here (and the names of procedures) were based on the d_set structure in LEDA, but I got
sick of debugging the LEDA class and emails sent to LEDA were ignored so 
designed my own.

*******************************************************/

#ifndef BIT_SET_H
#define BIT_SET_H

#include"global.h"

typedef unsigned long  word; 

/* Constants defined for quick access and comparison of words */
const word EMPTY = (word)0; 
const word FULL = ~EMPTY;
const word FIRST_BIT = (word)1; /* Smallest bit */
const int NUM_BITS = sizeof(word)*8;
const word LAST_BIT = ((word)1)<<(NUM_BITS-1); /* Largest bit */


class bit_set {
public:
  bit_set();
  bit_set(int cap); /* Create set with max_cap equal to cap */
  bit_set(const bit_set& bset);
  ~bit_set() {clear();} /* Clear set and free memory used */

  void copy(const bit_set& bset); 
  void clear(); /* Empties set and frees all memory */


  /* Standard operators */
  bit_set& operator=(const bit_set& bset); /* assign: error if bset has different capacity */
  int operator==(const bit_set& bset) const; /* test equality: error if bset has different capacity */

  /* access */
  int min() const; /* Minimum element: runs in time O(num_words + NUM_BITS). Error if empty. */
  int max() const; /* Maximum element: runs in time O(num_words + NUM_BITS). Error if empty. */	
  int member(int x) const; /* Test whether x is in the set */
  int size() const {return num_elements;} 
  inline int empty() const {return (num_elements==0);}
  int capacity() const {return max_cap;}

  /* Comparison */
  bool disjoint(const bit_set& bset) const; /* Checks if this set is disjoint with bset */
  bool operator<=(const bit_set& bset) const
    {return ((*this<bset)||(*this==bset));} /* lexical ordering */
  bool operator<(const bit_set& bset) const; /* lexical ordering: largest element has highest significance */
  
  /* manipulation */
  void insert(int x);
  void del(int x);
  void erase(); /* Makes the set empty, but does not affect capacity */
  void reserve(int cap); /* Sets the capacity of the set and allocates memory.  */
  void flip(); /* Computes the complement */

  /* Operators */
  bit_set operator!() const 
    { bit_set tmp(*this); tmp.flip(); return tmp;} 
  int operator!=(const bit_set& bset) const {return !((*this)==bset);}
  bit_set  operator+(const bit_set& bset) const /* union */
    { bit_set tmp(*this); tmp.join(bset); return tmp;}
  bit_set& operator+=(const bit_set& bset)
    { this->join(bset); return *this; }
  bit_set  operator&(const bit_set& bset) const /* intersection */
    { bit_set tmp(*this); tmp.intersect(bset); return tmp;}
  bit_set& operator&=(const bit_set& bset)
    { this->intersect(bset); return *this; }
  bit_set& operator-=(const bit_set& bset) /* difference */
    { this->diff(bset); return *this; }
  bit_set  operator-(const bit_set& bset) const
    { bit_set tmp(*this); tmp.diff(bset); return tmp; }
  
    
  
private:
  int max_cap;
  int num_words;
  int num_elements;
  vector<word> array;

  int manual_count(); /*recalculate number of elements */
  void intersect(const bit_set& bset);
  void join(const bit_set& bset);
  void diff(const bit_set& bset);

  friend ostream& operator<<(ostream& os, const bit_set& bset);
};

ostream& operator<<(ostream& os, const bit_set& bset);

#endif

