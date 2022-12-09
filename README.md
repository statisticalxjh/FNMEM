# FNMEM
Functional nonlinear mixed effects models

 1. If you want to speed up and have compilers such as Rtools or gcc, do not source "Kh.R" and "SCB.R", source "SCB.cpp" instead.
 2. To run as fast as possible, all the variables are declared with fixed sizes. If you want to change the parameters,
    please make corresponding modification in "SCB.cpp". Alternatively you can use the commented code instead, which will declare
    the variables corresponding to the parameters automaticly, but may be a little slower.
 3. If you want to change the nonlinear function, you have to make corresponding modification in "SCB.cpp".
