/*
==========================================================================
TeraPCA project.
Last update: 06 / 28 / 2018, USA.

 * Copyright (C) 2017-2018
 *  * All rights reserved.

==========================================================================

Authors:   Vassilis Kalantzis, kalan019@umn.edu
 * 	   Aritra Bose, bose6@purdue.edu
 * 	   Eugenia Kontopoulou, ekontopo@purdue.edu
==========================================================================
  * This program is free software; you can redistribute it and/or modify
  * it under the terms of the GNU General Public License as published by
  * the Free Software Foundation; either version 3 of the License, or
  * (at your option) any later version.
==========================================================================
*/

//==============================================================
// Header files
//==============================================================
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include "structures.h"
#include "utilities.h"
#include "methods.h"
#include "gaussian.h"
#include "mkl.h"
#include "omp.h"
#include "mkl_lapacke.h"
#include "gennorm.h"
#include "io.h"
//==============================================================

//==============================================================
//                 Define min, max routines
//==============================================================
#define max(a,b) (a>=b?a:b)
//==============================================================

//==============================================================
//                        Main driver
//==============================================================
int main(int argc, char **argv){

   //===========================================================
   // Declaration of variables
   //===========================================================
   double tt1, tt2;
   int ii, jj, kk;
   struct logistics logg;

   char fname[1024];
   char pname[1024];
   double fone = 1.0, fzero = 0.0;

   //===========================================================
   // Initialize logg structure
   //===========================================================
   initialize_structure(&logg);
   //===========================================================

   //============================================================================
   // This part is taken from EVSL -- it allows us to give the arguments in any 
   // particular order -- it requires a pre-defined value for each input argument
   //============================================================================
   int flg = findarg("help", NA, NULL, argc, argv);
   if (flg) {
     printf("\nUsage: ./TeraPCA.exe -bfile /path/to/matrix/ [char*] -nsv (default is 10) [int] -nrhs (default 2*nsv) [int] -rfetched [int] or -memory(in GB, default is 2) [int] -power [int] -print [int] -filewrite [int] -toll [double] -blockPower_maxiter [int] -blockPower_conv_crit [int]\n\n");
     printf("bfile: Filename containing binary genotypes along with the family and marker files \n");
     printf("nsv: # of leading left singular vectors sought\n");
     printf("rfetched: # of matrix rows per block \n");
     printf("memory: Allowed amount of memory (in GB)\n");
     printf("print: print statistics (1, 2) or not (0). Choice 2 prints approximate singular values and relative errors (the latter only if rfetched == rows)\n");
     printf("filewrite: write (1) or not (0)\n");
     printf("prefix: prefix for the name of the files\n");
     printf("toll: stopping tolerance for blockPower SVD\n");
     printf("blockPower_maxiter: max # of iterations in blockPower SVD \n");
     printf("blockPower_conv_crit: 0-> trace, 1-> individual \n\n");
     printf("power: applies Subspace Iteration to (AA')^power \n\n");
     printf("trueSVD: compute true SVD of A or not (applies only when the dataset is fully loaded in RAM) \n\n");
     return 0;
   }
   flg = findarg("about", NA, NULL, argc, argv);
   if (flg) {
     printf("\nAbout: Name of the library, version, last update (?), developers (?), funding (?)\n\n");
     return 0;
   }
   int pflag = findarg("prefix",STR, pname,argc,argv);
   findarg("bfile",STR, fname, argc, argv);
   findarg("nrhs", INT, &logg.NRHS, argc, argv);
   findarg("nsv", INT, &logg.NSV, argc, argv);
   findarg("memory", DOUBLE, &logg.mem, argc, argv);
   findarg("rfetched", INT, &logg.rows_fetched, argc, argv);
   findarg("print", INT, &logg.PRINT_INFO, argc, argv);
   findarg("filewrite", INT, &logg.filewrite, argc, argv);
   findarg("toll", DOUBLE, &logg.toll, argc, argv);
   findarg("blockPower_maxiter", INT, &logg.blockPower_maxiter, argc, argv);
   findarg("blockPower_conv_crit", INT, &logg.blockPower_conv_crit, argc, argv);
   findarg("power", INT, &logg.power, argc, argv);
   findarg("trueSVD", INT, &logg.trueSVD, argc, argv);
   //============================================================================

   std::string bfile(fname);
   if(pflag)
     logg.prefixname = pname; // assign the prefix name to string in logg
   logg.filename = bfile.c_str();
   logg.pure_name = ExtractFileName(logg.filename); // This stores the name of the file without the extension

   //===========================================================
   // Read BED files
   //===========================================================
   std::string strb(".bed");
   std::string bedfile = bfile+".bed";

   std::ifstream bedin(bedfile.c_str(), std::ios::in | std::ios::binary);
   bedin.seekg(3, std::ifstream::beg);  //The first 3 bytes are magic bytes in PLINK BED format
   if(!bedin){
       std::string err = std::string("[Data::read_bed] Error reading file ")
	 + bedfile;
       throw std::runtime_error(err);
   }
   //===========================================================

   //===========================================================
   // Read FAM and BIM files
   //===========================================================
   string famfile = bfile + ".fam";
   string bimfile = bfile + ".bim";
   logg.show_timestamp = 1;
   string strf(".fam");
   string strbim(".bim");
   //===========================================================

   //===========================================================
   // Read Fam file
   //===========================================================
   unsigned int line_num = 0;
   int nrows = -1;
   ifstream famin(famfile.c_str(), ios::in);
   cout <<  endl << timestamp(&logg) << "Reading .fam file: " << famfile << endl;
   if(!famin){
     string err = string("Error reading file '")
       + famfile + "': " + strerror(errno);
     cout <<  endl << err <<  endl;
     throw  runtime_error(err);
   }
   vector< string> famlines;

   while(famin){
     string line;
     getline(famin, line);
     if(!famin.eof() && (nrows == -1 || line_num < nrows)){
       if(line_num >= 0)
	 famlines.push_back(line);
       line_num++;
     }
   }
   GetFamInfo(famlines, &logg);
   logg.M = line_num;  //Number of Individuals
   famin.close();
   //===========================================================

   //===========================================================
   // Read Bim file
   //===========================================================
   line_num = 0;
   nrows = -1;
   ifstream bimin(bimfile.c_str(),  ios::in);
   cout <<  endl << timestamp(&logg) << "Opened .bim file: " << bimfile <<  endl;
   if(!bimin){
     string err =  string("Error reading file '")
       + bimfile + "': " + strerror(errno);
     throw  runtime_error(err);
   }
   vector< string> bimlines;
   while(bimin){
     string line1;
     getline(bimin, line1);
     if(!bimin.eof() && (nrows == -1 || line_num < nrows)){
       if(line_num >= 0)
	 bimlines.push_back(line1);
       line_num++;
     }
   }
   GetBimInfo(bimlines, &logg);
   logg.N = line_num;
   bimin.close();
   //===========================================================

   //===========================================================
   // Determine size of RAM
   //===========================================================
   int     ram_KB;
   double  ram_GB;
   ram_KB = GetRamInKB();
   std::cout << "Size of RAM (in KB): " << ram_KB << std::endl;
   ram_GB = (double) ram_KB;
   ram_GB = ram_GB / 1000000;
   logg.ram_KB = ram_KB;
   logg.ram_GB = ram_GB;
   //===========================================================

   //===========================================================
   // Extract number of OMP/MKL threads
   //===========================================================
   if (getenv("OMP_NUM_THREADS")) {
     logg.threads = atoi(getenv("OMP_NUM_THREADS"));
   }
   //===========================================================

   //===========================================================
   // Compute the rows to be fetched, dependent on the user defined RAM memory used
   //===========================================================
   double blksize;
   if (logg.rows_fetched <= 0){
	//Set aside a memory buffer from the allocated memory to accommodate the transient computations of the eigenvectors
    	double membuff = (3*(logg.N*8.0)) + (logg.M*logg.NSV*8.0) + (2*(logg.N*logg.NRHS)) + (3*logg.M*8.0) + 2048*100000;
	//Compute the workable memory
        double workmem = (logg.mem*1000000000) - membuff;
        if (workmem > 0)
           blksize = workmem/(8.0*logg.M);
        else
	   blksize = (membuff + (logg.M*8.0))/10000000;

   	logg.rows_fetched = (int)blksize;
  	if (logg.rows_fetched <= 0)
		logg.rows_fetched = logg.NSV;
   }
   // check if there is enough space to store the entire matrix
   if (logg.rows_fetched >= logg.N) {
     logg.rows_fetched = logg.N;
   }
   std::cout << std::endl << "Number of rows per block: " << logg.rows_fetched << std::endl;
   //=========================================================

   //=========================================================
   // Check # of rows/columns of matrix A
   //=========================================================
   if (logg.M <= 0 || logg.N <= 0) {
     printf("M and/or N were either zero or negative. Aborting...\n");
     exit(1);
   }
   //=========================================================

   //=========================================================
   // # of singular pairs sought
   //=========================================================
   if (logg.NSV <= 0) {
     printf("NSV was either zero or negative. Aborting...\n");
     exit(1);
   }
   if ( logg.NSV > min(logg.M,logg.N) ) {
     logg.NSV = min(logg.M,logg.N);
     printf("The value of NSV given was larger than min(M,N). Adjusting to NSV=min(M,N)...\n");
   }
   //==========================================================

   //==========================================================
   // # of columns in the initial subspace
   // =========================================================
   if ( logg.NRHS <= 0 ) {
     printf("NRHS was either zero or negative. Adjusting to NRHS=min(2xNSV,min(M,N))...\n");
     logg.NRHS = min(2*logg.NSV,min(logg.M,logg.N));
   }
   if ( logg.NRHS > min(logg.M,logg.N) ) {
     logg.NRHS = min(logg.M,logg.N);
     printf("The value of NRHS given was larger than min(M,N). Adjusting to NRHS=min(M,N)...\n");
   }
   if ( logg.NRHS < logg.NSV ) {
     printf("NRHS can not be smaller than NSV. Adjusting to NRHS=min(2xNSV,min(M,N))...\n");
     logg.NRHS = min(2*logg.NSV,min(logg.M,logg.N));
   }
   //==========================================================

   //==========================================================
   // If rows_fetched == rows --> load entire matrix
   //==========================================================
   double *MAT;
   long int offset = 0;
   if (logg.rows_fetched == logg.N) {
     uint64_t malloc_size = (uint64_t) logg.N*logg.M*sizeof(double);

     MAT = (double*) malloc(malloc_size);

     if (MAT==NULL) {
       printf("MAT malloc failed\n");
       exit(1);
     }

     tt1 = dsecnd();
     cout <<  endl << timestamp(&logg) << "Reading .bed file: " << bedfile << endl;
     Read_Bed(bedin,MAT,&logg);

     tt2 = dsecnd() - tt1;
     logg.TIME_2_LOAD_MATRIX = tt2;
   }
   //==========================================================

   //==========================================================
   // Fill the RHS matrix with normal random numbers
   //==========================================================
   double  mean, std_dev, norm_rv;
   mean    = 0;
   std_dev = 1;
   double *RHS = (double*)malloc(logg.M*logg.NRHS*sizeof(double));

   tt1 = dsecnd();
   for ( jj = 0; jj < logg.M; jj++ ) {
     for ( kk = 0; kk < logg.NRHS; kk++ ) {
       rand_val(jj*logg.M+kk);
       norm_rv = norm2(mean, std_dev);
       RHS[jj*logg.NRHS+kk] = (double) norm_rv;
     }
   }
   tt2 = dsecnd() - tt1;
   logg.TIME_2_GENERATE_RHS = tt2;
   //==========================================================

   //===========================================================
   // Compute leading left singular vectors
   //===========================================================
   if (logg.rows_fetched == logg.N) {
     subspaceIteration(MAT, RHS, &logg);
   } else {
     cout <<  endl << timestamp(&logg) << "Reading .bed file: " << bedfile << " by blocks." << endl;
     BlockSubspaceIter(bedin, RHS, &logg);
   }
   free(RHS);
   bedin.close();
   //==========================================================

   //==========================================================
   // Print approximate singular
   //==========================================================
   if (logg.PRINT_INFO > 1) {
     for (jj = 0; jj < logg.NSV; jj++) {
       printf("Approx. sing. value %d: %02.13f\n", jj, logg.sing_values[jj]);
     }
   }
   //===========================================================

   //===========================================================
   // Compute true SVD (if the entire matrix can fit in RAM)
   //===========================================================
   int min_dim = min(logg.M,logg.N);
   double* TRUE_SING_VALUES;
   double* TRUE_LEFT_SING_VECS;
   double* TRUE_RIGHT_SING_VECS;
   double* TRUE_superb;
   double* sing_vecs_relerror;
   double* sing_vals_relerror;
   double* copy_singular_vectors1;
   double* copy_singular_vectors2;
   double* UhatU;
   double* CosineValues;

   if ( logg.rows_fetched == logg.N && logg.trueSVD == 1 ){
     //===========================================================
     // Declarations
     //===========================================================
     double fone = 1.0;
     sing_vals_relerror = new double[logg.NSV];
     mkl_dimatcopy('R', 'T', logg.N, logg.M, fone, MAT, logg.M, logg.N);
     TRUE_SING_VALUES = (double*) malloc(min_dim*sizeof(double));
     TRUE_LEFT_SING_VECS = (double*) malloc(logg.M*min_dim*sizeof(double));
     TRUE_RIGHT_SING_VECS = (double*) malloc(logg.N*min_dim*sizeof(double));
     TRUE_superb = (double*) malloc((min_dim-1)*sizeof(double));
     sing_vecs_relerror = (double*) malloc(logg.M*logg.NSV*sizeof(double));
     copy_singular_vectors1 = (double*) malloc(logg.M*logg.NSV*sizeof(double));
     copy_singular_vectors2 = (double*) malloc(logg.M*logg.NSV*sizeof(double));
     UhatU = (double*) malloc(logg.NSV*logg.NSV*sizeof(double));
     //===========================================================

     //===========================================================
     // Compute SVD
     //===========================================================
     tt1 = dsecnd();
     int info_svd_lapacke = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'S', 'S', logg.M, logg.N, MAT, logg.N, TRUE_SING_VALUES, TRUE_LEFT_SING_VECS, min_dim, TRUE_RIGHT_SING_VECS, logg.N, TRUE_superb );
     tt2 = dsecnd()-tt1;
     logg.TIME_2_TRUE_SVD = tt2;
     if (logg.PRINT_INFO > 1)
       printf("True SVD was computed by DGESVD in: %02.13f seconds, code returned: %d \n", logg.TIME_2_TRUE_SVD, info_svd_lapacke);
     //===========================================================

     //===========================================================
     // Compute the (componentwise) relative error of
     // the singular values and vectors
     //===========================================================
     for ( ii = 0; ii < logg.NSV; ii++ ) {
       sing_vals_relerror[ii] = fabs(TRUE_SING_VALUES[ii]-logg.sing_values[ii])/TRUE_SING_VALUES[ii];
       if (logg.PRINT_INFO > 1)
	 printf("Rel. error of approx singval %d: %02.13f. True singval: %02.13f\n", ii, sing_vals_relerror[ii], TRUE_SING_VALUES[ii]);
     }
     for ( ii = 0; ii < logg.M; ii++ ) {
       for ( jj = 0; jj < logg.NSV; jj++ ) {
	 sing_vecs_relerror[ii*logg.NSV+jj] = fabs(fabs(TRUE_LEFT_SING_VECS[ii*logg.M+jj])-fabs(logg.left_sing_vecs[ii*logg.NSV+jj]))/fabs(TRUE_LEFT_SING_VECS[ii*logg.M+jj]);
	 if(logg.PRINT_INFO > 1)
	   printf("Rel. error of approx singvec %d: %02.13f\n", ii, sing_vecs_relerror[ii*logg.NSV+jj]);
       }
     }
     //===========================================================

     //===========================================================
     // Write True Left Singular Vectors and Singular Values into file
     //===========================================================
     if (logg.filewrite == 1) {
       string tempname1,tempname2;
       if (logg.prefixname.empty()){
	 tempname1 = ConstructFilename(logg,"realLeftSingularVectors");
	 tempname2 = ConstructFilename(logg,"realSingularValues");
       }
       else{
	 tempname1 = logg.prefixname + "_realLeftsingularVectors.txt";
	 tempname2 = logg.prefixname + "_realSingularValues.txt";
       }
       FILE *fwrite_realleftsingvecs = fopen(tempname1.c_str(), "w");
       if (fwrite_realleftsingvecs==NULL) {
	 printf("Unable to write to file. Aborting...");
	 exit(1);
       }
       for(ii = 0; ii < logg.M; ii++ ) {
	 for(jj = 0; jj < logg.NSV; jj++ )
	   fprintf(fwrite_realleftsingvecs, "% 2.13f   ", TRUE_LEFT_SING_VECS[ii*logg.M+jj]);
	 fprintf(fwrite_realleftsingvecs, "\n");
       }
       fclose(fwrite_realleftsingvecs);
       FILE *fwrite_realsingularvalues = fopen(tempname2.c_str(), "w");
       if (fwrite_realsingularvalues==NULL) {
	 printf("Unable to write to file. Aborting...");
	 exit(1);
       }
       for(jj = 0; jj < logg.NSV; jj++ )
	 fprintf(fwrite_realsingularvalues, "% 2.13f\n", TRUE_SING_VALUES[jj]);
       fclose(fwrite_realsingularvalues);
     }
     //===========================================================

     //===========================================================
     // Write singular vectors accuracy in file
     //===========================================================
     if (logg.filewrite == 1) {
       string tempname;
       if (logg.prefixname.empty())
	 tempname = ConstructFilename(logg,"singvecs_accuracy");
       else
	 tempname = logg.prefixname + "_singvecs_accuracy.txt";
       FILE *fwrite_errors = fopen(tempname.c_str(), "w");
       if (fwrite_errors==NULL) {
	 printf("Unable to write to file. Aborting...");
	 exit(1);
       }
       for(ii = 0; ii < logg.M*logg.NSV; ii++){
	 fprintf(fwrite_errors, "%2.13lf\n", sing_vecs_relerror[ii]);
       }
       fclose(fwrite_errors);
     }
     //===========================================================

     //===========================================================
     //
     //===========================================================
     for ( ii = 0; ii < logg.M; ii++ ) {
       for ( jj = 0; jj < logg.NSV; jj++ ) {
	 copy_singular_vectors1[ii*logg.NSV+jj] = logg.left_sing_vecs[ii*logg.NSV+jj];
	 copy_singular_vectors2[ii*logg.NSV+jj] = TRUE_LEFT_SING_VECS[ii*logg.M+jj];
       }
     }
     //===========================================================

     //===========================================================
     // Compute Uhat^T x U
     //===========================================================
     cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, logg.NSV, logg.NSV, logg.M, fone, copy_singular_vectors1, logg.NSV, copy_singular_vectors2, logg.NSV, fzero, UhatU, logg.NSV);
     //===========================================================

     //===========================================================
     // Compute Frobenius norm
     //===========================================================
     double frob_norm = 0;
     for ( ii = 0; ii < logg.NSV; ii++ ){
       for ( jj = 0; jj < logg.NSV; jj++ ){
	 if ( ii == jj )
	   frob_norm = frob_norm + (fabs(UhatU[ii*logg.NSV+jj])-1.0)*(fabs(UhatU[ii*logg.NSV+jj])-1.0);
	 else
	   frob_norm = frob_norm + UhatU[ii*logg.NSV+jj]*UhatU[ii*logg.NSV+jj];
       }
     }
     //===========================================================

     //===========================================================
     // Angle
     //===========================================================
     logg.frob_norm_angle = sqrt(frob_norm);
     //===========================================================

     //===========================================================
     // Compute cosine error per principal direction
     //===========================================================
     CosineValues = (double*) malloc(logg.NSV*sizeof(double));
     logg.cos_values.resize(logg.NSV);
     computeCosineError(copy_singular_vectors1, copy_singular_vectors2, logg.M, logg.NSV, CosineValues, &logg.cos_error);
     if (logg.PRINT_INFO > 1){
       printf("Cosine of true and approximate leading singular vectors\n");
       printf(" - - -- - - -- - - -- - - -- - - -- - - -- - - -- - -- - -\n");
     }
     for (ii = 0; ii < logg.NSV; ii++){
       if (logg.PRINT_INFO > 1)
	 printf("CosineValue(%d): %lf\n", ii, CosineValues[ii]);
       logg.cos_values[ii] = CosineValues[ii];
     }
     //===========================================================

     //===========================================================
     // Store cosine values in file
     //===========================================================
     if (logg.filewrite == 1){
       string tempname;
       if (logg.prefixname.empty())
         tempname = ConstructFilename(logg,"cosineValues");
       else
         tempname = logg.prefixname + "_cosineValues.txt";
       FILE *fwrite_cosinevalues = fopen(tempname.c_str(), "w");
       if (fwrite_cosinevalues==NULL){
	 printf("Unable to write to file. Aborting...");
	 exit(1);
       }
       for(int ii = 0; ii < logg.NSV; ii++)
	 fprintf(fwrite_cosinevalues, "% 2.13lf\n", logg.cos_values[ii]);
       fclose(fwrite_cosinevalues);
      }
      //===========================================================

      //===========================================================
      // Deallocation
      //===========================================================
      free(TRUE_SING_VALUES);
      free(TRUE_LEFT_SING_VECS);
      free(TRUE_RIGHT_SING_VECS);
      free(TRUE_superb);
      free(sing_vecs_relerror);
      free(sing_vals_relerror);
      free(copy_singular_vectors1);
      free(copy_singular_vectors2);
      free(UhatU);
      free(CosineValues);
      //===========================================================
   }
   //===========================================================

   //===========================================================
   // Store date and time of current simulation
   //===========================================================
   time_t rawtime;
   time ( &rawtime );
   logg.timeinfo = localtime ( &rawtime );
   //===========================================================

   //===========================================================
   // Print program statistics
   //===========================================================
   if (logg.PRINT_INFO > 0)
      print_statistics(logg);
   //===========================================================

   //===========================================================
   if (logg.rows_fetched == logg.N)
     free(MAT);
   //===========================================================

   return 0;
}
