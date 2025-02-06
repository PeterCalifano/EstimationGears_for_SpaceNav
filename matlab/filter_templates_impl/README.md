# How to use these templates and folder?
All the filters in this repository (at least the most updated) uses the same standard functions calls and 
# interfaces. This is done to allow modularity and re-use of the generic filter algorithms.
# The key idea is to provide overriding functionality without having to change the implementation of the 
# filter mechanization. Instead, tailoring of the algorithm is performed through the filter configuration 
# structure, as well as the various input struct (that MUST be specified in a case by case basis). 
# Similarly, the inner implementation of the modules likely requires tailoring (i.e., specifying observation 
# models and dynamical models).
# For instance, to do the latter, you need to tailor the following functions:
# computeDynFcn, computeDynMatrix, computeMeasPred, computeObsMatrix.
# In practice, everything would be simpler using an Object Oriented Programming paradigm. However, the 
# function-based approach was experimented to keep compatibility with MATLAB code generation and SLX.
