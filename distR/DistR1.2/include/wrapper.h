/* Copyright November 2005: Rachel Bevan
   Email: rachel@mcb.mcgill.ca
   NOTE: NOT Done commenting
*/

/*
 *
 * This module is just a bunch of wrappers for standard library
 * functions that can set errno and/or cause major malfunctions.  It's
 * tiresome to check for errors from them as well.  As of now, the
 * default behaviour is to exit when it fails.  I think this should
 * hold for files and memory allocation, because if they fail, there
 * is no way for the program to finish.  File opening should be done
 * before any calculations happen to prevent things like doing all the
 * calculations then finding out an output file won't open, thus
 * wasting all the CPU time.
 *
 * Fgetpos and Fsetpos do not call exit.
 */

#ifndef SW_WRAPPERS_MODULE
#define SW_WRAPPERS_MODULE

#include <stdio.h>
#include <stdlib.h>
extern size_t alloc_size();
extern FILE *Fopen (char *path, char *mode);
extern char Fclose(FILE *handle);
extern int Fgetpos (FILE *stream, fpos_t *pos);
extern int Fsetpos (FILE *stream, fpos_t *pos);
extern void *Calloc (size_t nmemb,size_t size);
extern void *Malloc (size_t size);
extern void *Realloc (void *ptr, size_t old_size,size_t new_size);
extern void Free(void *ptr, size_t size);
#endif /* SW_WRAPPERS_MODULE */

