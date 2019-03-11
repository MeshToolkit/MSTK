/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>
#include "MSTK_private.h"

/* Adapted from the Wikipedia page on quicksort */

/* These are specialized quicksort functions that will sort an integer
   array A according to an auxilliary integer key. It is useful when
   you cannot use the C qsort function because you don't want to form
   structures embedding the auxilliary information in each array entry
   or use global information. */


int mstk_quicksort_partition(int A[], int key[], int i_beg, int i_end, 
                             int i_pivot);


/*

  A     (in)  --- Array to be sorted
  key   (in)  --- Key values for array entries (so one looks at key[A[i]]
                  No check if the key array actually is valid and contains
                  values for all values of A[i]. Also, clearly A must be
                  an array of unsigned ints. *** CAN BE NULL ***
  i_beg (in)  --- beginning index in A (usually 0 at the topmost level)
  i_end (in)  --- ending index in A (inclusive) (usually N-1 at the topmost level)

*/

void mstk_quicksort(int A[], int key[], int i_beg, int i_end) {
  int n, tmp;

  n = i_end-i_beg+1;

  if (n <=1) return; 

  if (n == 2) {
    if (key) {
      if (key[A[i_beg]] > key[A[i_end]]) {
        tmp = A[i_beg];
        A[i_beg] = A[i_end];
        A[i_end] = tmp;
      }
    }
    else {
      if (A[i_beg] > A[i_end]) {
        tmp = A[i_beg];
        A[i_beg] = A[i_end];
        A[i_end] = tmp;
      }
    }
  }
  else {
    int p, nlo, nhi;

    p = i_beg + (int) (i_end-i_beg+1)/2;
    p = mstk_quicksort_partition(A, key, i_beg, i_end, p);
    
    nlo = p-i_beg+1;
    nhi = i_end-p;

    if (nlo < nhi) { /* sort the small sub-array first */      
      mstk_quicksort(A,key,i_beg,p-1);
      mstk_quicksort(A,key,p+1,i_end);
    }
    else {
      mstk_quicksort(A,key,p+1,i_end);
      mstk_quicksort(A,key,i_beg,p-1);
    }

  }
}

/* Partition an array to be sorted about a pivot point with values less 
   than key[A[i_pivot]] being to the left of the pivot and values greater
   than key[A[i_pivot]] being to the right

  A     (in)  --- Array to be sorted
  key   (in)  --- Key values for array entries (so one looks at key[A[i]]
                  No check if the key array actually is valid and contains
                  values for all values of A[i]. Also, clearly A must be
                  an array of unsigned ints. *** CAN BE NULL ***
  i_beg (in)  --- beginning index in A 
  i_end (in)  --- ending index in A (inclusive) 
  i_pivot (int) --- pivot index
*/

int mstk_quicksort_partition(int A[], int key[], int i_beg, int i_end, 
                        int i_pivot) {
  int val_pivot;
  int tmp, i, i_store;

  val_pivot = key ?  key[A[i_pivot]] : A[i_pivot];
 
  /* Move pivot value to rightmost point and out of the way */

  tmp = A[i_end];
  A[i_end] = A[i_pivot];
  A[i_pivot] = tmp;

  /* iteratively move elements smaller than
     pivot value to the left part of the array */

  i_store = i_beg;
  for (i = i_beg; i <= i_end-1; i++) {
    if (key) {
      if (key[A[i]] < val_pivot) {
        tmp = A[i];
        A[i] = A[i_store];
        A[i_store] = tmp;
        i_store++;
      }
    }
    else {
      if (A[i] < val_pivot) {
        tmp = A[i];
        A[i] = A[i_store];
        A[i_store] = tmp;
        i_store++;
      }
    }
  }

  /* Move the pivot back into position */

  tmp = A[i_store];
  A[i_store] = A[i_end];
  A[i_end] = tmp;

  return i_store; /* new position of the pivot */
}
