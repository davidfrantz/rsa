/**+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

This file is part of FORCE - Framework for Operational Radiometric 
Correction for Environmental monitoring.

Copyright (C) 2013-2022 David Frantz

FORCE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

FORCE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with FORCE.  If not, see <http://www.gnu.org/licenses/>.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/

/**+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Handle tables (csv-styled)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/


#ifndef TABLE_CL_H
#define TABLE_CL_H

#include <stdio.h>   // core input and output functions
#include <stdlib.h>  // standard general utilities library
#include <stdbool.h> // boolean data type

#include "const.h"
#include "alloc.h"
#include "string.h"
#include "stats.h"


#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  int nrow;
  int ncol;
  bool has_row_names;
  bool has_col_names;
  char **row_names;
  char **col_names;
  double **data;
  bool *row_mask;
  bool *col_mask;
  int n_active_cols;
  int n_active_rows;
  double *mean;
  double *sd;
  double *min;
  double *max;
  double *sum;
} table_t;

table_t read_table(char *fname, bool has_row_names, bool has_col_names);
table_t allocate_table(int nrow, int ncol, bool has_row_names, bool has_col_names);
int find_table_col(table_t *table, const char *name);
int find_table_row(table_t *table, const char *name);
void init_table(table_t *table);
void print_table(table_t *table, bool truncate, bool skip_rows);
void write_table(table_t *table, char *fname, const char *separator, bool skip_rows);
void free_table(table_t *table);

#ifdef __cplusplus
}
#endif

#endif


