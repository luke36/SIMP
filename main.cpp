#include "mesh.hpp"
#include <fstream>
#include <sstream>

int main(int argc, char *argv[]) {
  if (argc < 4) {
    std::cerr
      << "Usage: <executable> <input file> <output file prefix> <ratio[,ratio]*> <threshold>"
      << std::endl;
    exit(1);
  }
  std::ifstream in(argv[1]);
  std::istringstream ratio_list(argv[3]);
  std::istringstream threshold(argv[4]);

  std::vector<real> ratios;
  real ratio;
  ratio_list >> ratio;
  ratios.emplace_back(ratio);
  while (ratio_list.good()) {
    char comma;
    ratio_list >> comma >> ratio;
    ratios.emplace_back(ratio);
  }
  real thres;
  threshold >> thres;
  Mesh m(in);
  m.simplify(
             [&argv](Mesh &m, real ratio) {
               std::ostringstream path;
               path << argv[2] << '_' << ratio << ".obj";
               std::ofstream out(path.str());
               m.dump(out);
               out.close();
             },
             ratios, thres);
  in.close();
}
