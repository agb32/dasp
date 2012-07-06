//
// Matrix - vector multiplication.
//
// Author: Urban Bitenc
// Date:   June 2012
//
// This function is called by function "dot" in cmod/src/utils.c
//
#include "mvm.h"

void mvm( mvm_data_t* args )
{
  int n    = args->n;
  int m    = args->m;
  float* M = args->M;
  float* V = args->V;
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
