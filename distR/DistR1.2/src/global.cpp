/*******************************************************
filename: 	global.cpp
author: 	David Bryant (bryant@math.mcgill.ca)
created: 	12 Aug 1999
copyright:	David Bryant 1999

Utility routine(s)

error_message(" ") prints error message and aborts.

********************************************************/

#include"global.h"


void error_message(const string& message) {
  cerr<<endl;
  cerr<<"ERROR: "<<message<<endl;
  cerr<<"\t Please contact me (David Bryant) and tell me what happened: bryant@mcb.mcgill.ca"<<endl;
  exit(1);
}


