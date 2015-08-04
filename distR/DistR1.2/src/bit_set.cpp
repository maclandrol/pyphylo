/*******************************************************
filename: 	bit_set.cpp
author: 	David Bryant (bryant@math.mcgill.ca)
created: 	June 2001
copyright:	David Bryant 2001 

********************************************************/

#include "bit_set.h"

/* Default constructor: creates a set with zero capacity */
bit_set::bit_set() {
  max_cap=0;
  num_words=0;
  num_elements=0;
  array.clear();
 }

/* Preferred constructor - creates a set with given capacity */
bit_set::bit_set(int cap) {  
# ifdef BITSET_BOUNDS
  if (cap<0) error_message("Illegal capacity passed to bit_set");
# endif
  max_cap=0;
  num_words=0;
  num_elements=0;
  array.clear();
  reserve(cap);
}

/* Constructor that makes an exact copy of another bit_set */ 
bit_set::bit_set(const bit_set& bset) {
  array.clear();
  copy(bset);
}

/* Copy: this procedure CAN change the capacity of the set */
void bit_set::copy(const bit_set& bset) {
  max_cap=bset.max_cap;
  num_words = bset.num_words;
  num_elements = bset.num_elements;
  array.resize(num_words);
  //TO-DO check that array was allocated OK
  std::copy(bset.array.begin(),bset.array.end(),array.begin());
}

/* Clear everything in the set and set capacity to zero. If you want to just remove the elements, use erase */
void bit_set::clear() {
  max_cap=0;
  num_words=0;
  num_elements=0;
  array.clear();
}



/* Assignment */
bit_set& bit_set::operator=(const bit_set& bset) {
  if (this!=&bset) {
    copy(bset);
  }
  return *this;
}

/* Equality test 

The equality test (and the other comparisons) assume that in the last word of the array, those bits
corresponding to indices greater than capacity are set to zero. 
*/
int bit_set::operator==(const bit_set& bset) const 
{
#ifdef BITSET_BOUNDS 
  if (max_cap!=bset.max_cap)
    error_message("Testing equality between two sets with different universes");
#endif
  vector<word>::const_iterator p,q;
  if (this==&bset)  return true;
  if (size()!=bset.size()) return false;
  for (p=array.begin(),q=(bset.array).begin();p!=array.end(); p++,q++) {
    if ((*p)!=(*q))
      return false;
  }
  return true; /* Assumes that last word has extra bits set to zero */
}

/* minimum of a set 

We search for the first non-zero word then look for the smallest non-zero bit 
*/

int bit_set::min() const {

  int i=0,j;
  word theword;
  vector<word>::const_iterator p=array.begin();

  if (empty())
    error_message("Programming stuff up: Trying to compute minimum of an empty set");

  while(*p==EMPTY) {/* Guaranteed not to overflow since num_elements > 0 */
	p++;
	i++;
  }

  theword = *p;
  for (j=0; (theword & FIRST_BIT)==EMPTY; theword>>=1, j++) {} /* find the rightmost bit flagged */  
  return i*(NUM_BITS) + j;
}

/* minimum of a set 

We search for the last non-zero word then look for the largest non-zero bit 
*/

int bit_set::max() const {

  int i,j;
  word theword;
  vector<word>::const_iterator p=array.end();

  if (num_elements==0)
    error_message("Trying to compute maximum of an empty set");

  p--; /* If array.end()=array.begin() then num_elements=0 */
  i=num_words-1;
  
  while(*p==EMPTY) { /* Guaranteed not to overflow since num_elements > 0 */
	p--;
	i--;
  }

  theword = *p;
  for (j=NUM_BITS-1; (theword & LAST_BIT)==EMPTY; theword<<=1, j--) {} /* find the leftmost bit flagged */

  return i*(NUM_BITS) + j;
}

/* Test whether x is a member of the set */
int bit_set::member(int x) const {
  int i,j;

  if ((x<0)||(x>=max_cap))
	return false;

  if (empty())
	return false;
  
  i=x / NUM_BITS;
  j=x % NUM_BITS;
  return (array[i]&(1<<j));
}


/* Determine whether or not this set intersects with bset */
bool bit_set::disjoint(const bit_set& bset) const {

vector<word>::const_iterator p,q;

#ifdef BITSET_BOUNDS
 if (max_cap!=bset.max_cap)
   error_message("Computing disjoint between sets with different capacities");
#endif
 for (p=array.begin(), q=bset.array.begin();p!=array.end(); p++,q++)
   if ( ((*p)&(*q))!=EMPTY) 
     return false;
 return true;
}
		
/* Lexical ordering

Returns true if bset comes after *this in the lexical ordering.

The largest element has highest significance, so
it is as if we are turning the set into numbers 

n(X) = \sum_{i \in X} 2^i

and comparing n(X) and n(Y)

*/

bool bit_set::operator<(const bit_set& bset) const {

#ifdef BITSET_BOUNDS
  if (max_cap!=bset.max_cap)
    error_message("Doing a lexical comparison between sets with different capacities");
#endif
  
  vector<word>::const_iterator p,q; // Maybe better to use reverse iterators?    
  if (num_words == 0)
    return false; /* Two empty sets with 0 capacity */
  p=array.end();
  q=bset.array.end();
  do {
      p--;
      q--;
      if ((*p) < (*q))
	return true;
      else if ((*p) > (*q))
	return false;
  } while (p!=array.begin());
  return false;
}


/* Insert an element into the set. Will return an error if the element>=max_cap */

void  bit_set::insert(int x) {
  int i,j;

#ifdef BITSET_BOUNDS
  if (x>=max_cap)
    error_message("Trying to insert an element that is out of bounds");
#endif

  //To Do: check whether we save anything by performing membership test explicitly (save recomputing i,j)
  if (!member(x)) {
	i=x / NUM_BITS;
	j=x % NUM_BITS;
	array[i]|=(1<<j);
	num_elements++;
  }
}

/* Remove an element from the set */
void  bit_set::del(int x) {
  int i,j;
  
  //To Do: check whether we save anything by performing membership test explicitly (save recomputing i,j)
  if (member(x)) {
	i=x / NUM_BITS;
	j=x % NUM_BITS;
	array[i]&=~(1<<j);
	num_elements--;
  }
}

/* Empty the set, but don't change the capacity */
void  bit_set::erase() {
  num_elements=0;
  fill(array.begin(),array.end(),EMPTY);
}

/* Modify the capacity of the set, freeing or allocating memory as necc */
void bit_set::reserve(int cap) {
  vector<word>::iterator p;

# ifdef BITSET_BOUNDS
  if (cap<0)
	error_message("Trying to create set with negative capacity");
# endif

  int old_cap = max_cap;
  if (max_cap==cap)
    return;
  max_cap = cap;
  int old_words=num_words;
  num_words = ((cap-1)/NUM_BITS)+1;
  array.resize(num_words);
  //To do: check that memory is available
  /* Initialise any new memory */
  if (num_words>old_words) {
    p = array.begin();
    advance(p,old_words);
    fill(p,array.end(),EMPTY);
    
    /* Number of elements will be unchanged */
  }
  else if (max_cap<old_cap) {
    /* We have potentially lost some elements */
    num_elements = manual_count();
  }
}




/*
Complement of the set 

Note:- we set to zero the bits of array[num_words-1] that do not correspond to
elements in the universe.

*/
void bit_set::flip() {
  int num_end_bits;
  word lastmask;
  vector<word>::iterator p;

  if (num_words==0)
    return;

  for(p=array.begin();p!=array.end();p++)
    (*p) = ~(*p);
  
  num_end_bits=(max_cap)%NUM_BITS; /* Number of bits in the last word that correspond to elements of universe*/
  if (num_end_bits!=0) {
    lastmask =  ~(FULL<< num_end_bits); /* Last mask has a 1 for all elements<capacity, and a 0 for the others*/
    p--;
    (*p) &= (lastmask);
  }
  num_elements = max_cap - num_elements;

}

/*recalculate number of elements */
int bit_set::manual_count(){  
  int count=0;
  vector<word>::iterator p;
  
  int j;


  word theword;
  for(p=array.begin();p!=array.end();p++) {
	theword = *p;
	if (theword!=EMPTY)
	  for(j=0;j<NUM_BITS;j++, theword>>=1) {
		if (theword&FIRST_BIT != EMPTY)
		  count++;
	  }
  }
  return count;
}

/* Make this set equal to the intersection of itself and bset */
void bit_set::intersect(const bit_set& bset) {
  
  vector<word>::iterator p;
  vector<word>::const_iterator q;
  
# ifdef BITSET_BOUNDS
  if (bset.capacity()!=max_cap)
    error_message("Trying to interset two sets with different universes");
# endif

  for(p=array.begin(),q=bset.array.begin();p!=array.end();p++,q++)
    (*p)&=(*q);

  num_elements = manual_count();

}

/* Make this set equal to the union of itself and bset */
void bit_set::join(const bit_set& bset) {

  vector<word>::iterator p;
  vector<word>::const_iterator q;
  
# ifdef BITSET_BOUNDS
  if (bset.capacity()!=max_cap)
    error_message("Trying to join two sets with different universes");
# endif

  for(p=array.begin(),q=bset.array.begin();p!=array.end();p++,q++)
    (*p)|=(*q);
  num_elements = manual_count();
}

/* Remove from this set all elements in bset */

void bit_set::diff(const bit_set& bset) {

  vector<word>::iterator p;
  vector<word>::const_iterator q;

# ifdef BITSET_BOUNDS
  if (bset.capacity()!=max_cap)
    error_message("Trying to subtract two sets with different universes");
# endif

  for(p=array.begin(),q=bset.array.begin();p!=array.end();p++,q++)
    (*p)&=~(*q);
  num_elements = manual_count();
}


/* Output the set - used mainly for debugging */
ostream& operator<<(ostream& os, const bit_set& bset) {
  int firstelem=true;
  int i;

  os<<"{";  
  for (i=0;i<bset.capacity();i++) {
	if (bset.member(i)) {
	  if (!firstelem)
		os<<",";
	  os<<i;
	  firstelem=false;
	}
  }
  os<<"}";
  return os;
}



	
  
													  


  

