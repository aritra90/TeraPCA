#include "methods.h"
#include "utilities.h"

void subspaceIteration(double *MAT, double *RHS2, logistics *logg) {
  //===========================================================
  int M = logg->M, N = logg->N, NRHS = logg->NRHS;
  int max_iter = logg->blockPower_maxiter;
  int powers = logg->power;
  int ii, jj, converged = 0, kk, ii2, jj2;
  double tt1, tt2;
  double fone = 1.0, fzero = 0.0;
  double w[NRHS];
  double* B    = (double*) malloc(NRHS*NRHS*sizeof(double));
  double* B2   = (double*) malloc(NRHS*NRHS*sizeof(double));
  double* RHS  = (double*) malloc(NRHS*N*sizeof(double));
  double* ARHS = (double*) malloc(NRHS*M*sizeof(double));
  double* SING_VALUES     = (double*) malloc(logg->NSV*sizeof(double));
  double* SING_VALUES_OLD = (double*) malloc(logg->NSV*sizeof(double));
  double* tau = (double*) malloc(NRHS*sizeof(double));
  int info, info_sgeqrf_lapacke, info_sorgqr_lapacke;
  double tr1[NRHS],tr2[NRHS];
  logg->delta_iter.resize(max_iter);
  logg->sing_values.resize(logg->NSV);
  logg->left_sing_vecs.resize(M*logg->NSV);
  //============================================================

  for (ii = 0; ii < max_iter; ii++) {
    //============================================================
    // Multiply A' by RHS2 from the right
    //============================================================
    tt1 = dsecnd();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, NRHS, M, fone, MAT, M, RHS2, NRHS, fzero, RHS, NRHS);
    tt2 = dsecnd() - tt1;
    logg->TIME_2_MM = logg->TIME_2_MM + tt2;
    logg->TIME_2_MM_A = logg->TIME_2_MM_A + tt2;
    //============================================================

    //============================================================           
    // Multiply A by RHS from the right
    //============================================================
    tt1 = dsecnd();
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, NRHS, N, fone, MAT, M, RHS, NRHS, fzero, ARHS, NRHS);
    tt2 = dsecnd() - tt1;
    logg->TIME_2_MM = logg->TIME_2_MM + tt2;
    logg->TIME_2_MM_A_TRANSPOSED = logg->TIME_2_MM_A_TRANSPOSED + tt2;
    //============================================================
    
    //============================================================
    //                       (AA')^{powers-1}
    //============================================================
    for ( kk = 0; kk < powers-1; kk++) {
      tt1 = dsecnd();
      info_sgeqrf_lapacke = LAPACKE_dgeqrf( LAPACK_ROW_MAJOR, M, NRHS, ARHS, NRHS, tau );
      info_sorgqr_lapacke = LAPACKE_dorgqr( LAPACK_ROW_MAJOR, M, NRHS, NRHS, ARHS, NRHS, tau );
      tt2 = dsecnd() - tt1;
      logg->TIME_2_GS = logg->TIME_2_GS + tt2;
      
      memcpy(RHS2,ARHS,M*NRHS*sizeof(double));
      
      tt1 = dsecnd();
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, NRHS, M, fone, MAT, M, RHS2, NRHS, fzero, RHS, NRHS);
      tt2 = dsecnd() - tt1;
      logg->TIME_2_MM = logg->TIME_2_MM + tt2;
      logg->TIME_2_MM_A = logg->TIME_2_MM_A + tt2;

      tt1 = dsecnd();
      cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, NRHS, N, fone, MAT, M, RHS, NRHS, fzero, ARHS, NRHS);
      tt2 = dsecnd() - tt1;
      logg->TIME_2_MM = logg->TIME_2_MM + tt2;
      logg->TIME_2_MM_A_TRANSPOSED = logg->TIME_2_MM_A_TRANSPOSED + tt2;
    }
    //============================================================

    //============================================================                               
    // Perform orthogonalization                                                                           
    //============================================================
    tt1 = dsecnd();
    info_sgeqrf_lapacke = LAPACKE_dgeqrf( LAPACK_ROW_MAJOR, M, NRHS, ARHS, NRHS, tau );
    info_sorgqr_lapacke = LAPACKE_dorgqr( LAPACK_ROW_MAJOR, M, NRHS, NRHS, ARHS, NRHS, tau );
    tt2 = dsecnd() - tt1;
    logg->TIME_2_GS = logg->TIME_2_GS + tt2;
    //============================================================
    
    //============================================================
    // Multiply A' by RHS2 from the right
    //============================================================
    tt1 = dsecnd();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, NRHS, M, fone, MAT, M, ARHS, NRHS, fzero, RHS, NRHS);
    tt2 = dsecnd() - tt1;
    logg->TIME_2_MM = logg->TIME_2_MM + tt2;
    logg->TIME_2_MM_A = logg->TIME_2_MM_A + tt2;
    //============================================================

    //============================================================
    // Now compute RHS'RHS
    //============================================================
    tt1 = dsecnd();
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, NRHS, NRHS, N, fone, RHS, NRHS, RHS, NRHS, fzero, B, NRHS);
    tt2 = dsecnd() - tt1;
    logg->TIME_2_MM = logg->TIME_2_MM + tt2;
    //============================================================

    //============================================================
    // Solve projected eigenvalue problem
    //============================================================
    tt1 = dsecnd();
    for (ii2 = 0; ii2 < NRHS; ii2++) {
      for (jj2 = 0; jj2 < NRHS; jj2++) {
        if ( jj2 < ii2 ) {
	  B[ii2*NRHS + jj2] = 0.0;
        }
      }
    }
 
    info = LAPACKE_dsyev( LAPACK_ROW_MAJOR, 'V', 'U', NRHS, B, NRHS, w );
    tt2 = dsecnd() - tt1;
    logg->TIME_2_PROJECTED_SVD = logg->TIME_2_PROJECTED_SVD + tt2;
    //===========================================================

    //===========================================================
    // Copy back to RHS2 -- be careful, singular vectors 
    // come from smallest to largest
    //===========================================================
    memcpy(B2,B,sizeof(double)*NRHS*NRHS);
    for (ii2 = 0; ii2 < NRHS; ii2++) {
      for (jj2 = 0; jj2 < NRHS; jj2++) {
        B[ii2*NRHS+jj2] = B2[ii2*NRHS+(NRHS-1-jj2)];
      }
    }
    //===========================================================

    //===========================================================
    // Now multiply RHS2 by the eigenvectors from the right
    //===========================================================
    tt1 = dsecnd();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, NRHS, NRHS, fone, ARHS, NRHS, B, NRHS, fzero, RHS2, NRHS);
    tt2 = dsecnd() - tt1;
    logg->TIME_2_MM = logg->TIME_2_MM + tt2;
    //==========================================================

    //==========================================================
    // Copy singular values and monitor trace
    //==========================================================
    for (jj2 = 0; jj2 < logg->NSV; jj2++ ) {
      SING_VALUES[jj2] = sqrt(w[NRHS-1-jj2]);
      logg->sing_values[jj2] = SING_VALUES[jj2];
      logg->delta_iter[ii] = logg->delta_iter[ii] + logg->sing_values[jj2];
      if (logg-> PRINT_INFO >1)
	printf("At iteration:%d -->sing.val: %d is %02.13f\n", ii, jj2, SING_VALUES[jj2]);
    }
    //==========================================================

    //==========================================================
    //          Check convergence
    //==========================================================
    if ( logg->blockPower_conv_crit == 0) { // checking convergence based on the partial sum of singular values
      if (ii==0) {
        if (logg->PRINT_INFO > 1) 
	  printf("Partial sum at iteration: %d --> %02.13f\n", ii, logg->delta_iter[ii]);
      } else {
        logg->blockPower_trace_error = fabs(logg->delta_iter[ii-1]-logg->delta_iter[ii])/logg->delta_iter[ii];
        if (logg->PRINT_INFO > 1) 
	  printf("Partial sum at iteration: %d --> %02.13f. Rel. error: %02.13f\n", ii, logg->delta_iter[ii], logg->blockPower_trace_error);
        if (logg->blockPower_trace_error <= logg->toll) {
	  break;
        }
      }
    } else { // checking convergence based on the relative convergence of each individual singular value
      if (ii > 0) {
        for (jj = 0; jj < logg->NSV; jj++ ) {
          if (fabs((SING_VALUES[jj]-SING_VALUES_OLD[jj])/SING_VALUES[jj]) <= logg->toll ) {
            converged++;
          }
        }
	if (logg-> PRINT_INFO >1)
	  printf("At iteration: %d, %d singular values converged\n", ii, converged);
        if (converged == logg->NSV) {
          break;
        } 
      }
    }

    // prepare logistics for next iteration
    converged = 0;
    for (jj = 0; jj < logg->NSV; jj++ ) {
      SING_VALUES_OLD[jj] = SING_VALUES[jj];
    }
    //========================================================== 
  }
  //============================================================

  if ( ii < max_iter ) {
    logg->blockPower_total_its = ii+1;
  } else {
    logg->blockPower_total_its = ii;
  }
  
  for(ii = 0; ii < M; ii++ ) {
    for(jj = 0; jj < NRHS; jj++ ) {
      if ( jj < logg->NSV ) {
	logg->left_sing_vecs[ii*logg->NSV+jj] = RHS2[ii*NRHS+jj];
      }
    }
  }

  //========================================================== 
  // write approx. singular values, left vectors to file
  //========================================================== 
   printf("\n%d\n", logg->filewrite); 
    if (logg->filewrite == 1) {
    	string tempname;
    if (logg->prefixname.empty())
      tempname = ConstructFilename(*logg,"singularValues");
    else
      tempname = logg->prefixname + "_singularValues.txt";
    
    FILE *fwrite_singvalues = fopen(tempname.c_str(), "w");
    if (fwrite_singvalues==NULL) {
      printf("Unable to write to file. Aborting...");
      exit(1);
    }  
    
    std::vector<double> singularvals = logg->sing_values;
    std::transform(singularvals.begin(), singularvals.end(), singularvals.begin(), computeSquare);
    std::vector<string> individs = logg->indiv_ids;

    fprintf(fwrite_singvalues, "EIGENVALUES\n\n"); 

    for(ii = 0; ii < logg->NSV; ii++ )
      fprintf(fwrite_singvalues, "%2.13lf\n", singularvals[ii]);  
    fclose(fwrite_singvalues);

    if (logg->prefixname.empty())
      tempname = ConstructFilename(*logg,"singularVectors");
    else
      tempname = logg->prefixname + "_singularVectors.txt";
    FILE *fwrite_singvecs = fopen(tempname.c_str(), "w");
    if (fwrite_singvecs==NULL) {
      printf("Unable to write to file. Aborting...");
      exit(1);
    }
	fprintf(fwrite_singvecs, "FID");
	for(jj = 0; jj < logg->NSV; jj++)
		fprintf(fwrite_singvecs, "   PC%d",jj);
	fprintf(fwrite_singvecs,"\n");
    for(ii = 0; ii < M; ii++ ){
		fprintf(fwrite_singvecs, "%34s   ", individs[ii].c_str());	
      for(jj = 0; jj < logg->NSV; jj++ ) 
	fprintf(fwrite_singvecs, "% 2.13f   ", RHS2[ii*NRHS+jj]);
      fprintf(fwrite_singvecs, "\n");
    }
    fclose(fwrite_singvecs);
  }
  //==========================================================

  //==========================================================
  // Finalize program and deallocate resources                                                                             
  //==========================================================                                                                   
  free(SING_VALUES);                                                                                                    
  free(SING_VALUES_OLD);                                                                                                    
  free(tau);
  //==========================================================
  
}

void BlockSubspaceIter(std::ifstream& infile, double *RHS2, logistics *logg) {

  //==========================================================                            
  int    M = logg->M, N = logg->N, NRHS = logg->NRHS;
  int    max_iter = logg->blockPower_maxiter, min_dim = min(N,NRHS), powers = logg->power;
  double tt1, tt2;
  int    ione = 1,   converged = 0, ii, jj, kk, ii2, jj2;
  double fone = 1.0, minusfone = -1.0, fzero = 0.0;
  double* ARHS            = (double*) malloc(M*NRHS*sizeof(double));
  double* SING_VALUES     = (double*) malloc(min_dim*sizeof(double));
  double* SING_VALUES_OLD = (double*) malloc(min_dim*sizeof(double));
  double* tau             = (double*) malloc(NRHS*sizeof(double));
  int info_svd_lapacke, info_sgeqrf_lapacke, info_sorgqr_lapacke;
  logg->sing_values.resize(logg->NSV);
  logg->delta_iter.resize(max_iter);
  logg->left_sing_vecs.resize(M*logg->NSV);
  //================================================================                                          

  //================================================================
  int rows_fetched      = logg->rows_fetched;
  int loops             = N / rows_fetched;
  int remaining_rows    = N - rows_fetched*loops;
  unsigned int ik;
  int      colss        = logg->M;
  double* RHS           = (double*) malloc(max(remaining_rows,rows_fetched)*NRHS*sizeof(double));
  double* B2            = (double*) malloc(NRHS*NRHS*sizeof(double));
  double* B2_duplicate  = (double*) malloc(NRHS*NRHS*sizeof(double));
  double tr1[NRHS],tr2[NRHS];
  unsigned int np       = (unsigned long long)ceil((double)M/PACK_DENSITY);  //size of the packed data, in bytes, per SNP
  unsigned int actual_block_size=0, startval; 
  double *LOC_MAT;
  unsigned char *decbin  = (unsigned char*)malloc(np*PACK_DENSITY*sizeof(unsigned char*));
  unsigned char *readbin = (unsigned char*)malloc(np*sizeof(unsigned char*));
  double *norm_tmp = (double*)malloc(M*sizeof(double)); 
  //================================================================
  if (remaining_rows > 0) {
    loops = loops + 1;
  } else {
    remaining_rows = rows_fetched;
  }

  vector<int>start(loops);
  vector<int>stop(loops); 

  for(ik = 0 ; ik < loops ; ik++){
    start[ik]= ik * rows_fetched;
    stop[ik] = start[ik] + rows_fetched - 1;
    stop[ik] = stop[ik] >= N ? N - 1 : stop[ik];
  }

  LOC_MAT = (double*)malloc(max(remaining_rows,rows_fetched)*M*sizeof(double));
  //================================================================
  double *norm_precomp = (double*)malloc((4*N)*sizeof(double));
  memset(norm_precomp,0,(4*N)*sizeof(double));
  bool *seen_snp = new bool[N]();
  //================================================================
  
  //================================================================
  for (ii = 0; ii < max_iter; ii++) {
    //====================================
    // Multiply AA' by RHS2 from the right
    //====================================
    infile.seekg(3, std::ifstream::beg);
    memset(ARHS, 0, M*NRHS*sizeof(double));
    for (jj = 0; jj < loops; jj++) {
      // prepare to fetch data from memory
      // pack more bytes than required because it needs to read all of the byte
      if(jj<loops-1){
        actual_block_size = stop[jj] - start[jj] + 1;
      }else{
	actual_block_size = remaining_rows;
      }
      startval = start[jj];
      infile.seekg(3 + np * startval);

      if ( jj == 0 || jj < loops-1 ) {
	// Load chunk of SNPs       
        tt1 = dsecnd();
	Read_Bed_Blocks(infile, np, actual_block_size, LOC_MAT, startval, logg, decbin, readbin, norm_tmp, seen_snp, norm_precomp);
        tt2 = dsecnd() - tt1;
        logg->TIME_2_LOAD_MATRIX = logg->TIME_2_LOAD_MATRIX + tt2;
	
        // Multiply with A'
	tt1 = dsecnd();
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rows_fetched, NRHS, M, fone, LOC_MAT, M, RHS2, NRHS, fzero, RHS, NRHS);	
	tt2 = dsecnd() - tt1;
	logg->TIME_2_MM   = logg->TIME_2_MM + tt2;
	logg->TIME_2_MM_A = logg->TIME_2_MM_A + tt2;
	
        // Multiply with A
	tt1 = dsecnd();
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, NRHS, rows_fetched, fone, LOC_MAT, M, RHS, NRHS, fone, ARHS, NRHS);
	tt2 = dsecnd() - tt1;
	logg->TIME_2_MM = logg->TIME_2_MM + tt2;
	logg->TIME_2_MM_A_TRANSPOSED = logg->TIME_2_MM_A_TRANSPOSED + tt2;} 
      else{
        // Load chunk of SNPs
      	tt1 = dsecnd();
      	Read_Bed_Blocks(infile, np, actual_block_size, LOC_MAT, startval, logg, decbin, readbin, norm_tmp, seen_snp, norm_precomp);
        tt2 = dsecnd() - tt1;
      	logg->TIME_2_LOAD_MATRIX = logg->TIME_2_LOAD_MATRIX + tt2;
      	
        // Multiply with A'
      	tt1 = dsecnd();
      	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, remaining_rows, NRHS, M, fone, LOC_MAT, M, RHS2, NRHS, fzero, RHS, NRHS);
      	tt2 = dsecnd() - tt1;
      	logg->TIME_2_MM = logg->TIME_2_MM + tt2;
      	logg->TIME_2_MM_A = logg->TIME_2_MM_A + tt2;
      	
        // Multiply with A
      	tt1 = dsecnd();
      	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, NRHS, remaining_rows, fone, LOC_MAT, M, RHS, NRHS, fone, ARHS, NRHS);
      	tt2 = dsecnd() - tt1;
      	logg->TIME_2_MM = logg->TIME_2_MM + tt2;
      	logg->TIME_2_MM_A_TRANSPOSED = logg->TIME_2_MM_A_TRANSPOSED + tt2;
      }
    }

    //======================================================================
    // Trick to reduce the number of times that A is fetched from the memory
    //======================================================================
    for ( kk = 0; kk < powers-1; kk++) {
	//==========================
	// Perform orthogonalization
	//==========================
	tt1 = dsecnd();
	info_sgeqrf_lapacke = LAPACKE_dgeqrf( LAPACK_ROW_MAJOR, M, NRHS, ARHS, NRHS, tau );
	info_sorgqr_lapacke = LAPACKE_dorgqr( LAPACK_ROW_MAJOR, M, NRHS, NRHS, ARHS, NRHS, tau );
	tt2 = dsecnd() - tt1;
	logg->TIME_2_GS = logg->TIME_2_GS + tt2;
	//
	memcpy(RHS2,ARHS,M*NRHS*sizeof(double));
        memset(ARHS, 0,  M*NRHS*sizeof(double));
	//==========================
	
	//====================================
	// Multiply AA' by RHS2 from the right
	//====================================
	infile.seekg(3, std::ifstream::beg);
	for (jj = 0; jj < loops; jj++) {
	  //pack more bytes than required because it needs to read all of the byte
          if(jj<loops-1){
            actual_block_size = stop[jj] - start[jj] + 1;
          }else{
            actual_block_size = remaining_rows;
          }
	  startval = start[jj];
	  infile.seekg(3 + np * startval);
	  if ( jj == 0 || jj < loops-1 ) {
	     
	    // Load chunk of SNPs
	    tt1 = dsecnd();
	    Read_Bed_Blocks(infile, np, actual_block_size, LOC_MAT, startval, logg, decbin, readbin, norm_tmp, seen_snp, norm_precomp);
	    tt2 = dsecnd() - tt1;
	    logg->TIME_2_LOAD_MATRIX = logg->TIME_2_LOAD_MATRIX + tt2;
	    
	    // Multiply with A'
	    tt1 = dsecnd();
	    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rows_fetched, NRHS, M, fone, LOC_MAT, M, RHS2, NRHS, fzero, RHS, NRHS);
	    tt2 = dsecnd() - tt1;
	    logg->TIME_2_MM = logg->TIME_2_MM + tt2;
	    logg->TIME_2_MM_A = logg->TIME_2_MM_A + tt2;
	 
	    // Multiply with A
	    tt1 = dsecnd();
	    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, NRHS, rows_fetched, fone, LOC_MAT, M, RHS, NRHS, fone, ARHS, NRHS);
	    tt2 = dsecnd() - tt1;
	    logg->TIME_2_MM = logg->TIME_2_MM + tt2;
	    logg->TIME_2_MM_A_TRANSPOSED = logg->TIME_2_MM_A_TRANSPOSED + tt2;
	  } else {
	 
	    // Load chunk of SNPs
	    tt1 = dsecnd();
	    Read_Bed_Blocks(infile, np, actual_block_size, LOC_MAT, startval, logg, decbin, readbin, norm_tmp, seen_snp, norm_precomp);
	    tt2 = dsecnd() - tt1;
	    logg-> TIME_2_LOAD_MATRIX = logg->TIME_2_LOAD_MATRIX + tt2;
	    
	    // Multiply with A'
	    tt1 = dsecnd();
	    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, remaining_rows, NRHS, M, fone, LOC_MAT, M, RHS2, NRHS, fzero, RHS, NRHS);
	    tt2 = dsecnd() - tt1;
	    logg->TIME_2_MM = logg->TIME_2_MM + tt2;
	    logg->TIME_2_MM_A = logg->TIME_2_MM_A + tt2;

	    // Multiply with A
    	    tt1 = dsecnd();
    	    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, NRHS, remaining_rows, fone, LOC_MAT, M, RHS, NRHS, fone, ARHS, NRHS);
    	    tt2 = dsecnd() - tt1;
    	    logg->TIME_2_MM = logg->TIME_2_MM + tt2;
    	    logg->TIME_2_MM_A_TRANSPOSED = logg->TIME_2_MM_A_TRANSPOSED + tt2;
	  }
	}
    }
    //======================================================================

    //==========================
    // Perform orthogonalization
    //==========================
    tt1 = dsecnd();
    info_sgeqrf_lapacke = LAPACKE_dgeqrf( LAPACK_ROW_MAJOR, M, NRHS, ARHS, NRHS, tau );
    info_sorgqr_lapacke = LAPACKE_dorgqr( LAPACK_ROW_MAJOR, M, NRHS, NRHS, ARHS, NRHS, tau );
    tt2 = dsecnd() - tt1;
    logg->TIME_2_GS = logg->TIME_2_GS + tt2;
    //==========================

    //====================================
    // Multiply AA' by RHS2 from the right
    //====================================
    infile.seekg(3, std::ifstream::beg);
    memcpy(RHS2,ARHS,M*NRHS*sizeof(double));
    memset(ARHS, 0, M*NRHS*sizeof(double));
    for (jj = 0; jj < loops; jj++) {
        // prepare to fetch data from memory
        // pack more bytes than required because it needs to read all of the byte
        if(jj<loops-1){
          actual_block_size = stop[jj] - start[jj] + 1;
        }else{
          actual_block_size = remaining_rows;
        }
        startval = start[jj];
        infile.seekg(3 + np * startval);

        if ( jj == 0 || jj < loops-1 ) {
          // Load chunk of SNPs       
          tt1 = dsecnd();
  	  Read_Bed_Blocks(infile, np, actual_block_size, LOC_MAT, startval, logg, decbin, readbin, norm_tmp, seen_snp, norm_precomp);
          tt2 = dsecnd() - tt1;
          logg->TIME_2_LOAD_MATRIX = logg->TIME_2_LOAD_MATRIX + tt2;
	
          // Multiply with A'
	  tt1 = dsecnd();
	  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rows_fetched, NRHS, M, fone, LOC_MAT, M, RHS2, NRHS, fzero, RHS, NRHS);	
	  tt2 = dsecnd() - tt1;
	  logg->TIME_2_MM   = logg->TIME_2_MM + tt2;
	  logg->TIME_2_MM_A = logg->TIME_2_MM_A + tt2;
	
          // Multiply with A
	  tt1 = dsecnd();
	  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, NRHS, rows_fetched, fone, LOC_MAT, M, RHS, NRHS, fone, ARHS, NRHS);
	  tt2 = dsecnd() - tt1;
	  logg->TIME_2_MM = logg->TIME_2_MM + tt2;
	  logg->TIME_2_MM_A_TRANSPOSED = logg->TIME_2_MM_A_TRANSPOSED + tt2;} 
        else{
          // Load chunk of SNPs
      	  tt1 = dsecnd();
      	  Read_Bed_Blocks(infile, np, actual_block_size, LOC_MAT, startval, logg, decbin, readbin, norm_tmp, seen_snp, norm_precomp);
          tt2 = dsecnd() - tt1;
      	  logg->TIME_2_LOAD_MATRIX = logg->TIME_2_LOAD_MATRIX + tt2;
      	
          // Multiply with A'
      	  tt1 = dsecnd();
      	  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, remaining_rows, NRHS, M, fone, LOC_MAT, M, RHS2, NRHS, fzero, RHS, NRHS);
      	  tt2 = dsecnd() - tt1;
      	  logg->TIME_2_MM = logg->TIME_2_MM + tt2;
      	  logg->TIME_2_MM_A = logg->TIME_2_MM_A + tt2;
      	
          // Multiply with A
      	  tt1 = dsecnd();
      	  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, NRHS, remaining_rows, fone, LOC_MAT, M, RHS, NRHS, fone, ARHS, NRHS);
      	  tt2 = dsecnd() - tt1;
      	  logg->TIME_2_MM = logg->TIME_2_MM + tt2;
      	  logg->TIME_2_MM_A_TRANSPOSED = logg->TIME_2_MM_A_TRANSPOSED + tt2;
        }
    }
    //====================================

    //===================================================
    //                Now compute B'B
    //===================================================
    tt1 = dsecnd();
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, NRHS, NRHS, M, fone, RHS2, NRHS, ARHS, NRHS, fzero, B2, NRHS);
    tt2 = dsecnd() - tt1;
    logg->TIME_2_MM = logg->TIME_2_MM + tt2;
    //===================================================

    //=======================================
    // Solve the projected eigenvalue problem
    //=======================================
    tt1 = dsecnd();
    for (ii2 = 0; ii2 < NRHS; ii2++) {
      for (jj2 = 0; jj2 < NRHS; jj2++) {
        if ( jj2 < ii2 ) {
          B2[ii2*NRHS + jj2] = 0.0;
        }
      }
    }
    //
    double w[NRHS];
    int info = LAPACKE_dsyev( LAPACK_ROW_MAJOR, 'V', 'U', NRHS, B2, NRHS, w );
    tt2 = dsecnd() - tt1;
    logg->TIME_2_PROJECTED_SVD = logg->TIME_2_PROJECTED_SVD + tt2;
    //=======================================

    //==================================================================================
    // Be careful, singular vectors come from smallest to largest -- reverse their order
    //==================================================================================
    memcpy(B2_duplicate,B2,NRHS*NRHS*sizeof(double));
    for (ii2 = 0; ii2 < NRHS; ii2++) {
        for (jj2 = 0; jj2 < NRHS; jj2++) {
          B2[ii2*NRHS+jj2] = B2_duplicate[ii2*NRHS+(NRHS-1-jj2)];
        }
    }
    //==================================================================================
      
    //=====================================================
    // Now multiply RHS2 by the eigenvectors from the right
    //=====================================================
    memcpy(ARHS,RHS2,M*NRHS*sizeof(double));
    tt1 = dsecnd();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, NRHS, NRHS, fone, ARHS, NRHS, B2, NRHS, fzero, RHS2, NRHS);
    tt2 = dsecnd() - tt1;
    logg->TIME_2_MM = logg->TIME_2_MM + tt2;
    //=====================================================
      
    //=======================================
    // Copy singular values and monitor trace
    //=======================================
    for (jj = 0; jj < logg->NSV; jj++ ) {
	SING_VALUES[jj] = sqrt(w[NRHS-1-jj]);
	logg->sing_values[jj] = SING_VALUES[jj];
	logg->delta_iter[ii] = logg->delta_iter[ii] + logg->sing_values[jj];
	if (logg->PRINT_INFO >1){
	  printf("At iteration:%d -->sing.val: %d is %02.13f\n", ii, jj, SING_VALUES[jj]);}
    }
    //=======================================
      
    //=======================================
    //          Check convergence
    //=======================================
    if (logg->blockPower_conv_crit == 0) { // if checking convergence based on the sum
	if (ii==0) {
	  if (logg->PRINT_INFO > 1) {
	    printf("Partial sum at iteration: %d --> %02.13f\n", ii, logg->delta_iter[ii]);
	  }
	} else {
	  logg->blockPower_trace_error = fabs(logg->delta_iter[ii-1]-logg->delta_iter[ii])/logg->delta_iter[ii];
	  if (logg->PRINT_INFO > 1) {
	    printf("Partial sum at iteration: %d --> %02.13f. Rel. error: %02.13f\n", ii, logg->delta_iter[ii], logg->blockPower_trace_error);
	  }
	  if (logg->blockPower_trace_error <= logg->toll) {
	    break;
	  }
	}
      } else { // if checking convergence of each individual singular value
	if (ii > 0) {
	  for (jj = 0; jj < logg->NSV; jj++ ) {
	    if ( fabs((SING_VALUES[jj]-SING_VALUES_OLD[jj])/SING_VALUES[jj]) <= logg->toll ) {
	      converged++;
	    }
	  }
	  if (logg->PRINT_INFO >1){
	    printf("At iteration: %d, %d sing vals converged\n", ii, converged);}
	  if (converged == logg->NSV) {
	    break;
	  } 
	}
      }
      
      //=====================================
      // prepare logistics for next iteration
      //=====================================
      converged = 0;
      for (jj = 0; jj < logg->NSV; jj++ ) {
	SING_VALUES_OLD[jj] = SING_VALUES[jj];
      }
      //=====================================
  }
  //=======================================

  //=========================================================
  // If convergence achieved, increase iteration counter by 1
  //=========================================================
  if ( ii < max_iter ) {
    logg->blockPower_total_its = ii+1;
  } else {
    logg->blockPower_total_its = ii;
  }
  //=========================================================
    
  //=========================================================
  // Copy approx left singular vectors to RHS2
  //=========================================================
  for (ii = 0; ii < M; ii++ ) {
    for (jj = 0; jj < logg->NSV; jj++ ) {
      logg->left_sing_vecs[ii*logg->NSV+jj] = RHS2[ii*NRHS+jj];
    }
  }
  //=========================================================

  //============================================================================
  // If filewrite==1, write approximate singular values, left vectors to file
  //============================================================================
  if (logg->filewrite == 1) {
    string tempname;
    if (logg->prefixname.empty())
      tempname = ConstructFilename(*logg,"singularValues");
    else
      tempname = logg->prefixname + "_singularValues.txt";
    FILE *fwrite_singvalues = fopen(tempname.c_str(), "a");
    if (fwrite_singvalues==NULL) {
      printf("Unable to write to file. Aborting...");
      exit(1);
    }
    
    std::vector<double> singularvals = logg->sing_values;
    std::transform(singularvals.begin(), singularvals.end(), singularvals.begin(), computeSquare); 
    std::vector<string> individ = logg->indiv_ids; 
    
    fprintf(fwrite_singvalues, "EIGENVALUES\n\n"); 
    for(ii = 0; ii < logg->NSV; ii++ ){ 
      fprintf(fwrite_singvalues, "%2.13lf\n", singularvals[ii]);  
    }
    fclose(fwrite_singvalues);
      
    if (logg->prefixname.empty())
	tempname = ConstructFilename(*logg,"singularVectors");
    else
	tempname = logg->prefixname + "_singularVectors.txt";
    FILE *fwrite_singvecs = fopen(tempname.c_str(), "a");
    if (fwrite_singvecs==NULL) {
	printf("Unable to write to file. Aborting...");
	exit(1);
    }
	fprintf(fwrite_singvecs, "FID");
        for(jj = 0; jj < logg->NSV; jj++)
                fprintf(fwrite_singvecs, "   PC%d",jj);
        fprintf(fwrite_singvecs,"\n");

    for (ii = 0; ii < M; ii++ ) {
		fprintf(fwrite_singvecs, "%36s   ",individ[ii].c_str());
      for (jj = 0; jj < logg->NSV; jj++ ) {
	  fprintf(fwrite_singvecs, "% 2.13f   ", RHS2[ii*NRHS+jj]);
      }
	fprintf(fwrite_singvecs, "\n");
    }
    fclose(fwrite_singvecs);
  }   
  //============================================================================
   
  //==========================================                                 
  // Finalize program and deallocate resources 
  //==========================================                   
 
  free(LOC_MAT);
  free(ARHS);
  free(RHS); 
  free(SING_VALUES); 
  free(SING_VALUES_OLD); 
  free(B2);
  free(B2_duplicate);
  free(tau);
  free(readbin);
  free(decbin);
  free(norm_tmp);
  free(norm_precomp);
  delete[] seen_snp;
  start.clear();
  stop.clear(); 
  //==========================================   
}
