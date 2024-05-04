Requires CMake version 3.8 or higher, and a C++ compiler supporting C++17 or
higher.

Build (for example) by:
  $ mkdir build
  $ cd build
  $ cmake ..
  $ make -j4

Usage: <executable> <input file> <output file prefix> <ratio[,ratio]*> <threshold>

For example,
  $ ./main ../model/Armadillo.obj ../model/Armadillo_simp 0.5,0.2,0.1 0.1
This would generate ../model/Armadillo_simp_0.5.obj and so on.
