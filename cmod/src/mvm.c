/*
#dasp, the Durham Adaptive optics Simulation Platform.
#Copyright (C) 2004-2016 Alastair Basden and Durham University.

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as
#published by the Free Software Foundation, either version 3 of the
#License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

#You should have received a copy of the GNU Affero General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
//
// Matrix - vector multiplication.
//
// Author: Urban Bitenc
// Date:   June 2012
//
// This function is called by function "dot" in cmod/src/utils.c
//
#include "mvm.h"

void mvm_float( mvm_data_t* args )
{
  int n    = args->n;
  int m    = args->m;
  float* M = args->M;  // the only difference from mvm_double is that
  float* V = args->V;  //    these three variables are float and you add 6 terms per for loop iteration
  float* R = args->R;

  int i,j,cnt; // 'for'-loop indeces

  for(i=0, cnt=0; i<m; i++)
    {
      float tmp = 0;

      // Add 6 elements a time - this works fastest on MacBook (see UB's log book 1, p. 108, 2012Jun13)
      for(j=0; j<n-5;)
      	{
      	  tmp += M[cnt]*V[j] + M[cnt+1]*V[j+1] + M[cnt+2]*V[j+2] + M[cnt+3]*V[j+3] + M[cnt+4]*V[j+4] + M[cnt+5]*V[j+5];
      	  cnt += 6;
      	  j   += 6;
      	}

      // If N is not a multiple of X, add the last elements:
      for(; j<n;)
      	{
      	  tmp += M[cnt++]*V[j++];
      	}

      // Save the result into the result vector:
      R[i] = tmp;
    }
}

void mvm_double( mvm_data_t* args )
{
  int n    = args->n;
  int m    = args->m;
  double* M = args->M;  // the only difference from mvm_float is that
  double* V = args->V;  //    these three variables are double and you add 7 terms per for loop iteration
  double* R = args->R;

  int i,j,cnt; // 'for'-loop indeces

  for(i=0, cnt=0; i<m; i++)
    {
      float tmp = 0;

      // Add 7 elements a time - this works fastest on MacBook (see UB's log book 1, p. 146, 2012Jul09)
      for(j=0; j<n-6;)
      	{
      	  tmp += M[cnt]*V[j] + M[cnt+1]*V[j+1] + M[cnt+2]*V[j+2] + M[cnt+3]*V[j+3] + M[cnt+4]*V[j+4] + M[cnt+5]*V[j+5] + M[cnt+6]*V[j+6];
      	  cnt += 7;
      	  j   += 7;
      	}

      // If N is not a multiple of X, add the last elements:
      for(; j<n;)
      	{
      	  tmp += M[cnt++]*V[j++];
      	}

      // Save the result into the result vector:
      R[i] = tmp;
    }
}
