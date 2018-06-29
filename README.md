# TeraPCA
TeraPCA is a multithreaded C++ software suite based on Intel's MKL library (or any other BLAS and/or LAPACK distribution). TeraPCA features no dependencies to external libraries and combines the robustness of subspace iteration with the power of randomization.

To compile the package, just do: ```make```.

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

*print*: If print is set to 2, print details about convergence and matrix computations. Default is 1.  

*prefix*: Output filename prefix, default is the Binary Ped File name.  

*toll*: Tolerance criteria for   

*blockPower_maxiter*: Maximum iterations to run if convergence criterion is taking longer to achieve. 

*blockPower_conv_crit*: Change the convergence criterion to be used. 

**Usage**
```
Usage: ./TeraPCA.exe -bfile /path/to/matrix/ [char*] -nsv (default is 10) [int] -nrhs (default 2*nsv) [int] -rfetched [int] or -memory(in GB, default is 2) [int] -power [int] -print [int] -nthreads [int] -filewrite [int] -toll [double] -blockPower_maxiter [int] -blockPower_conv_crit [int]
```