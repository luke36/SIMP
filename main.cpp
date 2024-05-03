#include "mesh.hpp"
#include <fstream>
#include <sstream>

int main(int argc, char *argv[]) {
  if (argc < 4) {
    std::cerr << "Usage: <executable> <input file> <output file> <ratio> <tolerance>" << std::endl;
    exit(1);
  }
  std::ifstream in(argv[1]);
  std::ofstream out(argv[2]);
  std::istringstream ratio_list(argv[3]);
  std::istringstream threshold(argv[4]);

  real ratio;
  ratio_list >> ratio;
  real thres;
  threshold >> thres;
  Mesh m(in);
  m.simplify(ratio, thres);
  m.dump(out);
}
