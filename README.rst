DhedralPCA
===========

Performs PCA on dihedrals after translating the origin of each to the position with the smallest population. 

Dependecies
-------------

* Uses the Eigen template library for calculation of the eigensystem of the covariance matrix. (should be in deps/eigen)
* Uses ezOptionParser to ... parse options easily. (should be in deps/ezOptionParser)





Changelog
-------------

16-Jun-2013: The calculate dihedrals code works well now, but may not be the best option. Instead, it will accept one (or many) rama files (as generated from g_rama).
