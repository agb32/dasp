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
  void* M;  // input matrix
  void* V;  // input vector
  void* R;  // result vector
} mvm_data_t;

// function to perform the multiplication:
void mvm_float( mvm_data_t* arguments );  // for float input
void mvm_double( mvm_data_t* arguments ); // for double input
