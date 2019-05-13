# TeraPCA
TeraPCA is a multithreaded C++ software suite based on Intel's MKL library (or any other BLAS and/or LAPACK distribution). TeraPCA features no dependencies to external libraries and combines the robustness of subspace iteration with the power of randomization.

To install TeraPCA, one needs a working installation of Intel Math Kernel Library (MKL) (see https://software.intel.com/en-us/mkl), along with Intel’s “icpc“ C++ compiler (this can be downloaded as part of Intel Parallel Studio XE). No other external library/tool is necessary (the BLAS and LAPACK libraries are also used by TeraPCA and an implementation of these can be found in MKL).

Finally, one simply needs to edit the variable “MKL_ROOT“ in the Makefile equal to current MKL  installation directory. After this modification, move to the installation directory of TeraPCA and compile the package by doing: ```make```.

Once compiled run the program as following: ```./TeraPCA.exe -bfile Binary-PED-file ```
TeraPCA can be run with the following parameters.

**Mandatory**

*bfile*: Input Binary PED file whose principal components will be extracted. 

**Optional**  


*nsv*: Number of Principal Components to extract. Default is 10.

*nrhs*: Number of RHS in the sketching matrix. Default is 2*nsv.

*memory*: How much memory to be used by TeraPCA(in GB). Default is 2GB.

  OR
  
*rfetched*: Specify exactly the number of rows you want to extract.   

*power*: For faster computations and avoid fetching the matrix from memory repeatedly. Default is 1.  

*filewrite*: Boolean flag, when set to 1 will write two files with the singular values and the singular vectors, respectively. Default is 0.

*print*: If print is set to 2, print details about convergence and matrix computations. Default is 1.  

*prefix*: Output filename prefix, default is the Binary Ped File name.  

*toll*: Tolerance criteria for   

*blockPower_maxiter*: Maximum iterations to run if convergence criterion is taking longer to achieve. 

*blockPower_conv_crit*: Change the convergence criterion to be used. 

*trueSVD*: Compute true SVD using LAPACK of the matrix A(applies only when the dataset is fully loaded in RAM). 

Note: By default, TeraPCA links to the threaded layer of MKL. To use more than one threads in the MKL-related portions of TeraPCA simply set the environment variable ‘OMP_NUM_THREADS’ equal to the number of threads, i.e., type ‘export OMP_NUM_THREADS=x’ before execution, where ‘x’ denotes the number of threads used.

**Usage**
```
Usage: ./TeraPCA.exe -bfile /path/to/matrix/ [char*] -nsv (default is 10) [int] -nrhs (default 2*nsv) [int] -rfetched [int] or -memory(in GB, default is 2) [int] -power [int] -print [int] -filewrite [int] -toll [double] -blockPower_maxiter [int] -blockPower_conv_crit [int] -prefix [string]

```
An example dataset is given in the ``` example ``` directory with 10 individuals and 50 SNPs randomly chosen from the HapMap dataset. A sample output is also provided for the user to validate results. 

Contributors: Vassilis Kalantzis, Aritra Bose, Eugenia Kontopoulou
Any correspondence/questions about the codes are directed to: ```kalan019@umn.edu``` and/or ```bose6@purdue.edu```
