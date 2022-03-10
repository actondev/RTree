#include <array>
#include <catch2/catch.hpp>
#include <iostream>
#include <vector>
#include <RTree.h>
#include <easyprint.hpp>

using std::cerr;
using std::cout;
using std::endl;

template <typename T>
std::string pp(const T& x)
{
  return easyprint::stringify(x);
}

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
template <int DIMS>
std::vector<std::array<double, DIMS>> make_grid_template(int length)
{
  std::vector<std::array<double, DIMS>> res;
  int size = pow(length, DIMS);
  res.reserve(size);
  for (int i = 0; i < size; i++)
  {
    auto &el = res[i];
    for (int dim = 0; dim < DIMS; dim++)
    {
      int denum = pow(length, dim);
      int x = (int)(i / denum) % length;
      int mod = x % length;
      el[DIMS - dim - 1] = mod;
    }
    cerr << "  ";
    for(int i=0; i< DIMS; i++) {
      cerr << el[i] << ",";
    }
    cerr << endl;
  }

  return res;
}

std::vector<std::vector<double>> make_grid(int length, int dims = 2)
{
  std::vector<std::vector<double>> res;
  int size = pow(length, dims);
  res.resize(size);
  for (int i = 0; i < size; i++)
  {
    auto &el = res[i];
    el.resize(dims); // resize or reserve? kind of same?
    for (int dim = 0; dim < dims; dim++)
    {
      int denum = pow(length, dim);
      int x = (int)(i / denum) % length;
      int mod = x % length;
      el[dims - dim - 1] = mod;
    }
    cerr << "  ";
    for(int i=0; i< dims; i++) {
      cerr << el[i] << ",";
    }
    cerr << endl;
  }

  return res;
}

// template<typename T>
// void print_vec(std::ostream &os, const std::vector<T> &vec) {
//   auto size = vec.size();
//   // os << "[" << size << "]";
//   os << "{";
//   for(auto i=0; i < size - 1; i++) {
//     os << vec[i] << ",";
//   }
//   if(size > 0) {
//     os << vec[size-1];
//   }
//   os << "}";
// }

// std::ostream &operator<<(std::ostream &os, const std::vector<int> &vec) {
//   print_vec<int>(os, vec);
//   return os;
// }

// std::ostream &operator<<(std::ostream &os, const std::vector<double> &vec) {
//   print_vec<double>(os, vec);
//   return os;
// }

// std::ostream &operator<<(std::ostream &os, const std::vector<std::vector<double>> &vec) {
//   cerr << "{";
//   for(const auto& x : vec) print_vec<double>(os, x);
//   cerr << "}" << endl;
//   return os;
// }

TEST_CASE("rtree test init")
{

#define DIMS 2
  int length = 4;
  auto grid2 = make_grid(length, DIMS);

  using Type = std::vector<double>;
  RTree<Type, double, DIMS> tree;

  for(auto& el : grid2) {
    double* pos = &el[0];
    cerr << "inserting " << pp(el) << endl;
    tree.Insert(pos, pos, el);
  }
  REQUIRE(tree.Count() == pow(length, DIMS));

  std::vector<Type> found;
  using Callback =  std::function<bool(Type, const double*,const double*)>;
  Callback cb = [&found](Type node, const double* low, const double* high) {
    // cerr << "callback! " << pp(node) << endl;
    found.push_back(node);
    return true;
  };

  double low[DIMS] = {1,1};
  double high[DIMS] = {1,3};
  // we should get 1,1 1,2 1,3
  tree.Search(low, high, cb);
  REQUIRE(found.size() == 3);
  cerr << "found " << pp(found) << endl;
  std::sort(found.begin(), found.end(),
            [](const Type &a, const Type &b) {
              const auto size = std::min(a.size(), b.size());
              for(auto i=0; i < size; i++) {
                if ( a[i] < b[i]) {
                  return true;
                }
              }
              return false;
            });
  cerr << "sorted " << pp(found) << endl;
  std::vector<Type> expected = {{1,1}, {1,2}, {1,3}};
  REQUIRE(found == expected);
  
  cerr << "done?\n";
}
