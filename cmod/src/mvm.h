// 
// Performs matrix-vector multiplication:
// 
//
// Author: Urban Bitenc
// Date:   June 2012
//
// This function is called by function "dot" in cmod/src/utils.c
//

// data structure necessary to pass the arguments to a thread:
typedef struct {
  int n;     // input vector size
  int m;     // result vector size
  float* M;  // input matrix
  float* V;  // input vector
  float* R;  // result vector
} mvm_data_t;

// function to perform the multiplication:
void mvm( mvm_data_t* arguments );
