#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mkl.h"
#include "omp.h"
#include "structures.h"
#include "mkl_lapacke.h"
#include <string>
#include <sstream>

using namespace std;

string toString(int value) {
  ostringstream s;
  s << value;
  return s.str();
}

void decode_plink(unsigned char * __restrict__ out,
   const unsigned char * __restrict__ in,
   const unsigned int n);
void standardize(double *tmp, unsigned char *nnX, int n, int m, double &p, double &q);

void decode_plink_precomp(unsigned char * __restrict__ out,
   const unsigned char * __restrict__ in,
   const unsigned int n);
 
 
string timestamp(struct logistics *logg);
// The following routine was copied by the following link:
// http://stackoverflow.com/questions/349889/how-do-you-determine-the-amount-of-linux-system-ram-in-c
int GetRamInKB(void) {
  FILE *meminfo = fopen("/proc/meminfo", "r");

  char line[256];
  while(fgets(line, sizeof(line), meminfo)) {
    int ram;
    if(sscanf(line, "MemTotal: %d kB", &ram) == 1) {
      fclose(meminfo);
      return ram;
    }
  }

  fclose(meminfo);
  return -1;
}


// partialy copied from https://stackoverflow.com/questions/14143801/return-filename-from-path            
std::string ExtractFileName(const std::string& fullPath){
  const size_t lastSlashIndex = fullPath.find_last_of("/\\");
  const size_t lastDotIndex   = fullPath.substr(lastSlashIndex + 1).find_last_of(".");
  return fullPath.substr(lastSlashIndex+1,lastDotIndex);
}


string ConstructFilename(struct logistics logg, string fileType){
  string year  = toString(1900 + logg.timeinfo->tm_year);
  string month;
  if (logg.timeinfo->tm_mon + 1 < 9) 
    month = "0" + toString(logg.timeinfo->tm_mon + 1);
  else 
    month =  toString(logg.timeinfo->tm_mon + 1);
  string day   = toString(logg.timeinfo->tm_mday); 
  string atime = toString(logg.timeinfo->tm_hour) + "-" + toString(logg.timeinfo->tm_min) + "-" + toString(logg.timeinfo->tm_sec);
  string cname = logg.pure_name + "_" + year + month + day + "_" + atime + "_" + fileType + ".txt";
  return cname;
}

void print_statistics(struct logistics logg) {

  printf("===================================================================\n");
  printf("TeraPCA library, Version: 1.0.\n");
  printf( "Current local time and date: %s", asctime (logg.timeinfo) );
  printf("===================================================================\n");
  printf("# of threads exploited: %d\n",logg.threads);
  printf("RAM size in KBs: %d\n",logg.ram_KB);
  printf("RAM size in GBs: %lf\n",logg.ram_GB);
  printf("# of matrix rows: %d\n", logg.M);
  printf("# of matrix columns: %d\n", logg.N);
  printf("# of singular pairs sought: %d\n", logg.NSV);
  printf("# of right-hand sides used: %d\n", logg.NRHS);
  printf("Value of power: %d\n", logg.power);
  printf("# of rows fetched from the disk per block: %d\n", logg.rows_fetched);
  if (logg.rows_fetched < logg.N) {
    printf("Total number of times which the data matrix was fetched from the memory: %d\n", logg.blockPower_total_its*(logg.power+1));
    printf("Average amount of time elapsed per matrix fetching: %lf\n", logg.TIME_2_LOAD_MATRIX / (logg.blockPower_total_its*(logg.power+1)));
  }
  printf("===================================================================\n");
  printf("The following times are listed in seconds!\n");
  printf("- - - - -\n");
  printf("Time to generate the right-hand sides matrix: %02.13f\n",logg.TIME_2_GENERATE_RHS);
  printf("Time to load the data matrix: %02.13f\n",logg.TIME_2_LOAD_MATRIX);
  printf("Time to perform the MM products (overall): %02.13f\n",logg.TIME_2_MM);
  printf("Time to perform the MM products (with A): %02.13f\n",logg.TIME_2_MM_A);
  printf("Time to perform the MM products (with A^T): %02.13f\n",logg.TIME_2_MM_A_TRANSPOSED);
  printf("Time to perform the orthonormalization: %02.13f\n",logg.TIME_2_GS);
  printf("Time to solve the projection eigenvalue problem: %02.13f\n",logg.TIME_2_PROJECTED_SVD);
  printf("- - - - -\n");
  printf("Total wall-clock time elapsed: %02.13f\n",logg.TIME_2_PROJECTED_SVD + logg.TIME_2_GENERATE_RHS + logg.TIME_2_LOAD_MATRIX + logg.TIME_2_MM + logg.TIME_2_GS + logg.TIME_2_OTHER);
  printf("Total wall-clock time elapsed (without including the amount of time spent on loading the matrix): %02.13f\n",logg.TIME_2_PROJECTED_SVD + logg.TIME_2_GENERATE_RHS + logg.TIME_2_MM + logg.TIME_2_GS + logg.TIME_2_OTHER);
  printf("===================================================================\n");
  printf("The following times are listed in hours!\n");
  printf("- - - - -\n");
  printf("Time to generate the right-hand sides matrix: %02.13f\n",logg.TIME_2_GENERATE_RHS/3600);
  printf("Time to load the data matrix: %02.13f\n",logg.TIME_2_LOAD_MATRIX/3600);
  printf("Time to perform the MM products (overall): %02.13f\n",logg.TIME_2_MM/3600);
  printf("Time to perform the MM products (with A): %02.13f\n",logg.TIME_2_MM_A/3600);
  printf("Time to perform the MM products (with A^T): %02.13f\n",logg.TIME_2_MM_A_TRANSPOSED/3600);
  printf("Time to perform the orthonormalization: %02.13f\n",logg.TIME_2_GS/3600);
  printf("Time to solve the projection eigenvalue problem: %02.13f\n",logg.TIME_2_PROJECTED_SVD/3600);
  printf("- - - - -\n");
  printf("Total wall-clock time elapsed: %02.13f\n",(logg.TIME_2_PROJECTED_SVD + logg.TIME_2_GENERATE_RHS + logg.TIME_2_LOAD_MATRIX + logg.TIME_2_MM + logg.TIME_2_GS + logg.TIME_2_OTHER)/3600);
  printf("Total wall-clock time elapsed (without including the time spent on loading the matrix): %02.13f\n",(logg.TIME_2_PROJECTED_SVD + logg.TIME_2_GENERATE_RHS + logg.TIME_2_MM + logg.TIME_2_GS + logg.TIME_2_OTHER)/3600);
 
  if (logg.rows_fetched == logg.N && logg.trueSVD == 1) {
    printf("===================================================================\n");
    printf("Total wall-clock time to compute the full ('econ') SVD: %02.13f\n",logg.TIME_2_TRUE_SVD);
    printf("||Uhat^TU - I||_F: %02.13f\n",logg.frob_norm_angle);
    printf("||Uhat^TU - I||_2: %02.13f\n",logg.cos_error);
  }
  printf("===================================================================\n");

}


void initialize_structure(struct logistics *logg) {
  
  logg->filename = "noprefix";
  logg->prefixname = "noprefix";
  logg->pure_name = "noprefix";

  logg->M = 0; 
  logg->N = 0; 
  logg->NSV = 10; 
  logg->NRHS = 20; 
  logg->threads = 1; 
  logg->max_threads = 1;
  logg->rows_fetched = 0;
  logg->PRINT_INFO = 1;
  logg->filewrite = 0;
  logg->power = 1;
  logg->mem = 0;

  logg->blockPower_total_its = 0;
  logg->blockPower_trace_error = 0.0;
  logg->blockPower_maxiter = 100;
  logg->blockPower_conv_crit = 1;
  logg->toll = 1e-3;

  logg->TIME_2_GS = 0; 
  logg->TIME_2_LOAD_MATRIX=0;
  logg->TIME_2_MM = 0; 
  logg->TIME_2_MM_A_TRANSPOSED = 0;
  logg->TIME_2_MM_A = 0; 
  logg->TIME_2_PROJECTED_SVD = 0;
  logg->TIME_2_GENERATE_RHS = 0;
  logg->TIME_2_OTHER = 0;
  
  logg->frob_norm_angle = 0;
  logg->cos_error = 0;
  logg->trueSVD = 0;

}

void  computeCosineError(double *MatA, double *MatB, int M, int NSV, double *CosineValues, double *CosineError){
  /* 
     MatA         : M x NSV -- computed NSV right singular vectors
     MatB         : M x NSV -- real dominant NSV right  singular vectors
     M            : Number of rows of the matrices                                                                                                                                                                 
     NSV          : Number of columns of the matrices
     CosineValues : Computing the quantity diag(MatA'MatB)
     CosineError  : Computing the quantity ||MatA'MatB-I||_2
  */

  int ii, jj;
  double *SingularValues;
  double *MatC, *U, *V, *superb;

  MatC = (double*) malloc(NSV*NSV*sizeof(double));
  SingularValues   = (double*) malloc(NSV*sizeof(double));
  superb = (double*) malloc(NSV*sizeof(double));

  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,NSV,NSV,M,1,MatA,NSV,MatB,NSV,0,MatC,NSV);

  for ( ii = 0; ii < M; ii++ ){
    for ( jj = 0; jj < NSV; jj++){
      if ( ii == jj){
	CosineValues[jj] = MatC[ii*NSV+jj];
	MatC[ii*NSV+jj] = fabs(MatC[ii*NSV+jj])-1;
      }
    }
  }

  LAPACKE_dgesvd(LAPACK_ROW_MAJOR,'N','N',NSV,NSV,MatC,NSV,SingularValues,U,NSV,V,NSV,superb);


  *CosineError = SingularValues[0];

  free(MatC);
  free(SingularValues);
}

void Read_Bed(std::ifstream &in, double *temp3, struct logistics *logg){
    	int M = logg->M, N = logg->N; double avg, sd; 
	//================================================================
	//  //size of the packed data, in bytes, per SNP
	//================================================================
	uint64_t np = (unsigned long long)ceil((double)M/PACK_DENSITY);
	unsigned char *decbin = new unsigned char[np*PACK_DENSITY];
        unsigned char *readbin = new unsigned char[np];
        double *norm_tmp = new double[M];
	for(unsigned int j = 0; j < N; j++){
                // read raw genotypes
                in.read((char*)readbin, sizeof(char)*np);
                //decode raw genotypes
                decode_plink(decbin,readbin, np);
                standardize(norm_tmp,decbin,M,N, avg, sd);
		for(unsigned int kk = 0; kk < M; kk++){
		   double s = norm_tmp[kk];
		   temp3[j*M + kk] = s;
		}
	}
	delete[] norm_tmp;
	delete[] decbin;
	delete[] readbin;
}
void Read_Bed_Blocks(std::ifstream& in, uint64_t np, 
	uint64_t blk_size, double *temp3, uint64_t startval, 
	struct logistics *logg, unsigned char *decbin, unsigned char *readbin, 
	double *norm_tmp, bool *seen_snp, double *norm_precomp){
	int N = logg->M, M = logg->N; uint64_t idx; 
	double sd=0, avg=0, s;
	for(uint64_t j = 0; j < blk_size; j++){
		// read raw genotypes
		in.read((char*)readbin, sizeof(char)*np);
		idx = startval+j; 
	    if (!seen_snp[idx]){		
			//decode raw genotypes
				decode_plink(decbin,readbin, np);
			//normalize raw genotypes 
			standardize(norm_tmp, decbin, N, M, avg, sd);

			for(unsigned int kk = 0; kk < N; kk++){
				s = norm_tmp[kk];
				temp3[j*N + kk] = s;
			}

 			if(sd > 1e-9){
				norm_precomp[3+(idx*4)] = (0 - avg)/sd;
				norm_precomp[2+(idx*4)] = (1 - avg)/sd;
				norm_precomp[0+(idx*4)] = (2 - avg)/sd;
				norm_precomp[1+(idx*4)] = 0;
			}
			seen_snp[idx] = true; 
		}
		else{
			decode_plink_precomp(decbin, readbin, np);
			for(unsigned int kk = 0; kk < N; kk++){
				int b = (int)decbin[kk];
				s = norm_precomp[b+(idx*4)];
			//	std::cout << s << " ";
	    			temp3[j*N + kk] = s/sqrt(M); 
			}
		}
     }
}

/* decode_plink_precomp is taken from FlashPCA2's data.cpp where they read the Binary PED file 
   and decode it 
  */
  
void decode_plink_precomp(unsigned char * __restrict__ out,
   const unsigned char * __restrict__ in,
   const unsigned int n)
{
   unsigned int i, k;

   for(i = 0 ; i < n ; i++)
   {
      k = PACK_DENSITY * i;

      /* geno is interpreted as a char, however a1 and a2 are bits for allele 1 and
       * allele 2. The final genotype is the sum of the alleles, except for 01
       * which denotes missing.
       */

      out[k] =   (in[i] & MASK0);
      out[k+1] = (in[i] & MASK1) >> 2;
      out[k+2] = (in[i] & MASK2) >> 4;
      out[k+3] = (in[i] & MASK3) >> 6;
   }
}
/* decode_plink is taken from FlashPCA2's data.cpp where they read the Binary PED file 
   and decode it 
  */
void decode_plink(unsigned char * __restrict__ out,
   const unsigned char * __restrict__ in,
   const unsigned int n)
{
   unsigned int i, k;
   unsigned char tmp, geno1, geno2, geno3, geno4;
   unsigned int a1, a2;

   for(i = 0 ; i < n ; i++)
   {
      tmp = in[i];
      k = PACK_DENSITY * i;

      geno1 = (tmp & MASK0);
      if(geno1 == 1)
	 out[k] = 3;
      else
      {
	 a1 = !(geno1 & 1);
	 a2 = !(geno1 >> 1);
	 out[k] = a1 + a2;
      }
      k++;

      geno2 = (tmp & MASK1) >> 2;
      if(geno2 == 1)
	 out[k] = 3;
      else
      {
	 a1 = !(geno2 & 1);
	 a2 = !(geno2 >> 1);
	 out[k] = a1 + a2; 
      }
      k++;

      geno3 = (tmp & MASK2) >> 4;
      if(geno3 == 1)
	 out[k] = 3;
      else
      {
	 a1 = !(geno3 & 1);
	 a2 = !(geno3 >> 1);
	 out[k] = a1 + a2;
      }
      k++;

      geno4 = (tmp & MASK3) >> 6;
      if(geno4 == 1)
	 out[k] = 3;
      else
      {
	 a1 = !(geno4 & 1);
	 a2 = !(geno4 >> 1);
	 out[k] = a1 + a2;
      }
   }
}

/* 
	This is the normalization of each SNP vector across all individuals. 
	It calculates the mean/allele frequency and computes the standard deviation and 
	uses Price et al. (2006) standardization method, just to be consistent with the 
	same algorithm, in order to compare the eigenvalues
*/
void standardize(double *normX, unsigned char *nnX, int N, int M, double &avg, double &sd){
	double sum = 0;
	unsigned int ctr = 0;
	double all_freq, P,stddev; 
	double fact = 2.0;
	for(unsigned int im = 0 ; im < N ; im++)
	{
		double x = (double)nnX[im];
		if(nnX[im] != MISSING)
      		{
      	    		sum += x;
      	    		ctr++;
      	   	}
    	}		
	
    	all_freq = sum / ctr;	
	P = all_freq / 2.0;
	stddev = sqrt(fact * P * (1 - P));
    	avg = all_freq; 
	sd = stddev; 
	double newX;
	for(unsigned int jm = 0; jm < N; jm++)
	{
		double x = (double)nnX[jm];
		
		if (x == MISSING)
			normX[jm] = 0;
		else
			newX = (x - all_freq)/stddev;
			normX[jm] = newX/sqrt(M);
		if (isnan(normX[jm])){
			normX[jm] = 0; 
		}
	}
}

double computeSquare(double x) {return x*x;};

void GetBimInfo( vector< string> lines,struct logistics *logg){
	for(unsigned int i = 0 ; i < lines.size() ; i++)
   {
       stringstream ss(lines[i]);
       string s;
       vector< string> tokens;

      while(ss >> s)
	 tokens.push_back(s);
      logg->snp_ids.push_back(tokens[1]);
   }
}

void GetFamInfo(vector<string> lines,struct logistics *logg){
	for(unsigned int i = 0 ; i < lines.size() ; i++)
   {
       stringstream ss(lines[i]);
       string s;
       vector< string> tokens;
      while(ss >> s)
		tokens.push_back(s);
		logg->fam_ids.push_back(tokens[0]);
		logg->indiv_ids.push_back(tokens[1]);
   }
}
  
string timestamp(struct logistics *logg){
   if(logg->show_timestamp)
   {
      time_t t = time(NULL);
      char *s = asctime(localtime(&t));
      s[strlen(s) - 1] = '\0';
       string str(s);
      str =  string("[") + str +  string("] ");
      return str;
   }
   else
      return  string("");
}
