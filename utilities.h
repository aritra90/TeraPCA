#include "structures.h"

string toString(int value);

int GetRamInKB(void);

string ExtractFileName(const std::string& fullPath);

string ConstructFilename(struct logistics logg, string fileType);


void print_statistics(struct logistics logg);

void initialize_structure(struct logistics *logg);


void computeCosineError(double *MatA, double *MatB, int M, int NSV, double *CosineValues, double *CosineError);

void Read_Bed(std::ifstream &in, double *temp3, struct logistics *logg);


void Read_Bed_Blocks(std::ifstream& in, uint64_t x, uint64_t y, double *temp3, 
					uint64_t start, struct logistics *logg, unsigned char *dec, 
					unsigned char *read, double *tmp, bool *bb, double *nrm);

void decode_plink(unsigned char * __restrict__ out, const unsigned char * __restrict__ in, const unsigned int n);

void standardize(double* x, unsigned char *nnX, int n, int m, double b, double c);

double computeSquare (double x); 
void GetFamInfo(vector<string> famlines,struct logistics *logg);

void GetBimInfo(vector<string> bimlines,struct logistics *logg);

string timestamp(struct logistics *logg);
