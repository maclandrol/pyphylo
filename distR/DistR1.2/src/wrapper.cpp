#ifndef wrapper_c
#define wrapper_c

#include <errno.h>
#include "wrapper.h"

/* Copyright November 2005: Rachel Bevan
   Email: rachel@mcb.mcgill.ca
   NOTE: NOT Done commenting
*/
/* This file contains wrapper functions for some of the c commands which
   might return NULL.  They make it easier to signal this to the user */


static size_t allocated=0;
static int files_open=0;

size_t alloc_size()
{
  return allocated;
}

FILE *
Fopen (char *path,  char *mode)
{
  FILE *temp = fopen (path, mode);
  
  if (temp == NULL)
    {
      fprintf(stderr, "fopen: cannot open file %s",path);
      exit (1);
    }
  else{
    files_open++;
    return temp;
  }
}

char Fclose(FILE *handle)
{
  char f;
  
  f=fclose(handle);
  
  if(f == EOF)
    {
      fprintf(stderr,"Could not close file handle\n");
      exit(1);
    }
  files_open--;
  return f;
}

int
Fgetpos (FILE *stream, fpos_t *pos)
{
  int temp = fgetpos (stream, pos);
  
  if (temp != 0)
    {
      fprintf(stderr,"Error in fgetpos\n");
      return 1;
    }
  else
    return 0;
}

int
Fsetpos (FILE *stream, fpos_t *pos)
{
  int temp = fsetpos (stream, pos);
  
  if (temp != 0)
    {
      fprintf(stderr, "Error in fsetpos\n");
      return 1;
    }
  else
    return 0;
}

void *
Malloc (size_t size)
{
  void *temp = malloc (size);
    
  if (temp == NULL)
    {
      fprintf(stderr,"malloc:unable to allocate memory");
      exit (1);
    }
  else{
    allocated=allocated+size; 
    return temp;

  }
}

void *
Calloc (size_t nmemb,size_t size)
{
  void *temp = calloc (nmemb,size);
    
  if (temp == NULL)
    {
      fprintf(stderr,"calloc:unable to allocate memory");
      exit (1);
    }
  else{
    allocated=allocated+nmemb*size; 
    return temp;
  }
}

void *
Realloc (void *ptr, size_t old_size,size_t new_size)
{
  void *temp = realloc (ptr, new_size);
  
  if ((temp == NULL))
    {
      fprintf(stderr,"realloc:unable to allocate memory");
      exit (1);
    }
  else{
    allocated=allocated + new_size-old_size;
    return temp;
  }
}

void Free(void *ptr, size_t size)
{
  if(ptr!=NULL){
    free(ptr);
    allocated=allocated-size; 
  }
  else{
    fprintf(stderr, "attempting to free NULL pointer!\n");

  }
  return;
}

#endif







