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

#ifndef PLANCK_PARSER_H
#define PLANCK_PARSER_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

void *openParams (const char *file);
void closeParams (void *vpar);
int locateParam (void *vpar, const char *key);
int getVarInt (void *vpar, const char *key, int *result);
int getVarLongLong (void *vpar, const char *key, long long *result);
int getVarDouble (void *vpar, const char *key, double *result);
int getVarString (void *vpar, const char *key, char **result);
int getVarStringMem (void *vpar, const char *key, char **result);
int getVarStringProt (void *vpar, const char *key, char *result, size_t capacity);
int getVarStringUnProt (void *vpar, const char *key, char *result);

#ifdef __cplusplus
}
#endif

#endif
