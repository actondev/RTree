#include <catch2/catch.hpp>
#include <iostream>
#include <vector>

using std::cout;
using std::cerr;
using std::endl;


/**
 * see https://stackoverflow.com/a/49514906/8720686
 * When calling with <3>(2) will give back
  0,0,0,
  0,0,1,
  0,1,0,
  0,1,1,
  1,0,0,
  1,0,1,
  1,1,0,
  1,1,1,
 * 
 */
template<int DIMS>
std::vector<std::array<float, DIMS>> make_grid(int length) {
  std::vector<std::array<float, DIMS>> res;
  int size = pow(length, DIMS);
  res.reserve(size);
  for(int i = 0; i < size; i++) {
    auto& el = res[i];
    for(int dim = 0; dim < DIMS; dim++) {
      int denum = pow(length, dim);
      int x = (int)(i / denum ) % length;
      int mod = x % length;
      el[DIMS - dim - 1] = mod;
    }
    // cerr << "  ";
    // for(int i=0; i< DIMS; i++) {
    //   cerr << el[i] << ",";
    // }
    // cerr << endl;
  }

  return res;
}


TEST_CASE("rtree test init") {
  cerr << "hi there\n";
  auto grid = make_grid<3>(2);
  cerr << "done?\n";
}
