#include <stdio.h>              // Needed for printf()
#include <stdlib.h>             // Needed for exit() and ato*()
#include <math.h>               // Needed for sqrt() and log()

//----- Defines -------------------------------------------------------------
#define PI         3.14159265   // The value of pi

double norm2(double mean, double std_dev);  // Returns a normal rv
double rand_val(int seed);                 // Jain's RNG
