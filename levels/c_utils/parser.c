/*
 *  This file is part of libc_utils.
 *
 *  libc_utils is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libc_utils is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libc_utils; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libc_utils is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Parameter file parser
 *
 *  Copyright (C) 2008, 2009 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parser.h"
#include "c_utils.h"

typedef struct
  {
  size_t nparams;
  size_t reserved;
  char **key, **value;
  } params;

static void strtrim (char *string)
  {
  size_t len=strlen(string);
  size_t lw=0;
  while ((len>0) && ((string[len-1]==' ')
                  || (string[len-1]=='\t')
                  || (string[len-1]=='\n')))
    --len;
  string[len]=0;
  while ((string[lw]==' ') || (string[lw]=='\t') || (string[lw]=='\n'))
    ++lw;
  if (lw>0)
    memmove (string,string+lw,len+1-lw);
  }

static void push_params (params *par, const char *key, const char *value)
  {
  if (par->reserved<=par->nparams)
    {
    par->reserved*=2;
    if (par->reserved == 0) par->reserved=1;
    REALLOC(par->key,char *,par->reserved);
    REALLOC(par->value,char *,par->reserved);
    }
  par->key[par->nparams] = RALLOC(char,strlen(key)+1);
  strcpy (par->key[par->nparams],key);
  par->value[par->nparams] = RALLOC(char,strlen(value)+1);
  strcpy (par->value[par->nparams],value);
  ++par->nparams;
  }

void *openParams (const char *file)
  {
  FILE *fd;
  params *par;
  int cnt=0;
  char line[2048];

  fd = fopen(file,"r");
  if (!fd) return NULL;

  par = malloc (sizeof(params));
  par->nparams=0;
  par->reserved=0;
  par->key=NULL;
  par->value=NULL;

  while (fgets (line, 2000, fd)!=NULL)
    {
    size_t len;
    strtrim(line);
    ++cnt;
    len=strlen(line);

    if (len>0)
      {
      if (line[0]!='#')
        {
        char key[2048], value[2048];
        size_t eqpos = 0;
        while ((eqpos<len) && (line[eqpos]!='=')) ++eqpos;
        if ((eqpos==0)||(eqpos==len))
          {
          fprintf(stderr, "unrecognizable syntax at line %d:\n%s\n",cnt,line);
          fclose(fd); return par;
          }
        strncpy (key,line,eqpos);
        key[eqpos]=0;
        strcpy (value,line+eqpos+1);
        strtrim(key);
        strtrim(value);
        push_params(par,key,value);
        }
      }
    }
  fclose(fd);
  return par;
  }

void closeParams (void *vpar)
  {
  params *par=vpar;
  if (par!=NULL)
    {
    size_t m;
    for (m=0; m<par->nparams; ++m)
      {
      DEALLOC(par->key[m]);
      DEALLOC(par->value[m]);
      }
    DEALLOC(par->key);
    DEALLOC(par->value);
    DEALLOC(par);
    }
  }

int locateParam (void *vpar, const char *key)
  {
  params *par=vpar;
  size_t m;
  for (m=0; m<par->nparams; ++m)
    if (strcmp(key,par->key[m])==0) return m;
  return -1;
  }

int getVarInt (void *vpar, const char *key, int *result)
  {
  params *par=vpar;
  int idx = locateParam (par, key);
  char *endptr;
  if (idx<0) return 1;                     /* key not found */
  if (*(par->value[idx]) == 0) return 1;   /* empty value */
  *result = strtol(par->value[idx],&endptr,10);
  if ((*endptr)!=0) return 1;              /* did not parse entire string */
  return 0;
  }

int getVarLongLong (void *vpar, const char *key, long long *result)
  {
  params *par=vpar;
  int idx = locateParam (par, key);
  char *endptr;
  if (idx<0) return 1;                     /* key not found */
  if (*(par->value[idx]) == 0) return 1;   /* empty value */
  *result = strtoll(par->value[idx],&endptr,10);
  if ((*endptr)!=0) return 1;              /* did not parse entire string */
  return 0;
  }

int getVarDouble (void *vpar, const char *key, double *result)
  {
  params *par=vpar;
  int idx = locateParam (par, key);
  char *endptr;
  if (idx<0) return 1;                     /* key not found */
  if (*(par->value[idx]) == 0) return 1;   /* empty value */
  *result = strtod(par->value[idx],&endptr);
  if ((*endptr)!=0) return 1;              /* did not parse entire string */
  return 0;
  }

int getVarString(void* vpar, const char* key, char** result)
  {
  return getVarStringMem(vpar, key, result);
  }

int getVarStringMem (void *vpar, const char *key, char **result)
  {
  params *par=vpar;
  int idx = locateParam (par, key);
  if (idx<0) return 1;                     /* key not found */
  *result = malloc (strlen(par->value[idx])+1);
  if (!(*result)) return 1;                /* out of memory */
  strcpy (*result,par->value[idx]);
  return 0;
  }

int getVarStringProt (void *vpar, const char *key, char *result, size_t capacity)
  {
  params *par=vpar;
  int idx = locateParam (par, key);
  if (idx<0) return 1;                     /* key not found */
  if (strlen(par->value[idx])>capacity-1) return 1;   /* result does not fit */
  strcpy (result,par->value[idx]);
  return 0;
  }

int getVarStringUnProt (void *vpar, const char *key, char *result)
  {
  params *par=vpar;
  int idx = locateParam (par, key);
  if (idx<0) return 1;                     /* key not found */
  strcpy (result,par->value[idx]);
  return 0;
  }
