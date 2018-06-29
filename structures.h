#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <string.h>
#include <sstream>
#include <errno.h>
#include <stdexcept>
#include <iomanip>
#include <inttypes.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <ctype.h>

//==========================================
// Define macros
//==========================================
#define PACK_DENSITY 4
#define MISSING 3
#define MASK0 3	  /* 3 << 2 * 0 */
#define MASK1 12  /* 3 << 2 * 1 */
#define MASK2 48  /* 3 << 2 * 2 */
#define MASK3 192 /* 3 << 2 * 3 */

using namespace std;
using std::vector;

#ifndef LOGISTICS_STRUCT
#define LOGISTICS_STRUCT
struct logistics
{ 
  // Timing variable
  struct tm* timeinfo;

  string filename;
  string prefixname;
  string pure_name;

  // General variables
  int M, N, NSV, NRHS;
  int loops, rows_fetched;
  int PRINT_INFO;
  int power;
  double mem; 
  int    threads, max_threads;
  int    ram_KB;
  double ram_GB;


  // for blockPower
  int    blockPower_total_its;
  int    blockPower_maxiter;
  int    blockPower_conv_crit;
  double blockPower_trace_error; // returned error
  double toll;                   // used in blockPower to stop the iterations
  vector<double> delta_iter;

  // quantities returned
  vector<double> left_sing_vecs;
  vector<double> sing_values;
  vector<double> cos_values;

  // write singular pairs to text files or not
  int filewrite;

  // timers
  double TIME_2_GENERATE_RHS, TIME_2_LOAD_MATRIX, TIME_2_PROJECTED_SVD, TIME_2_TRUE_SVD, TIME_2_MM, TIME_2_GS, TIME_2_OTHER;
  double TIME_2_MM_A_TRANSPOSED, TIME_2_MM_A;


  // to hold the frobenius norm of Uhat^TU-I
  double frob_norm_angle;
  // to hold the cosine error
  double cos_error;
  int trueSVD;

	vector<string> indiv_ids;
	vector<string> fam_ids;
	vector<string> snp_ids;  
	bool show_timestamp ;
		
};
#endif


