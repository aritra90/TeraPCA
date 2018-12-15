#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "structures.h"
#include "mkl.h"
#include "mkl_lapacke.h"

//=========================
// Define min, max routines
//=========================
#define min(a,b) (a<=b?a:b)
#define max(a,b) (a>=b?a:b)
//=========================

void subspaceIteration(double *MAT, double *RHS2, struct logistics *logg);

void BlockSubspaceIter(std::ifstream& in, double *RHS2, logistics *logg);

void benchmarking(std::ifstream& in, double *RHS2, logistics *logg);
