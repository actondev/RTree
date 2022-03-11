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

using Point = std::vector<double>;
struct Grid {
  int dims;
  std::vector<Point> points;
};

using Callback = std::function<bool(Point, const double *, const double *)>;

/**
 * see https://stackoverflow.com/a/49514906/8720686
 * Given length=3,dims=2 it gives
 * {{0, 0}, {0, 1}, {0, 2},
 *  {1, 0}, {1, 1}, {1, 2},
 *  {2, 0}, {2, 1}, {2, 2}}
 */
Grid make_grid(int length, int dims = 2)
{
  Grid grid;
  grid.dims = 2;
  int size = pow(length, dims);
  grid.points.resize(size);
  for (int i = 0; i < size; i++)
  {
    auto &el = grid.points[i];
    // not optimal, but makes things simple (no template, runtime only)
    el.resize(dims);
    for (int dim = 0; dim < dims; dim++)
    {
      int denum = pow(length, dim);
      int x = (int)(i / denum) % length;
      int mod = x % length;
      el[dims - dim - 1] = mod;
    }
  }

  return grid;
}

RTree<Point> grid_to_rtree(const Grid& grid) {
  RTree<Point, double> tree(grid.dims);
    for (auto &el : grid.points) {
      const double *pos = &el[0];
      tree.Insert(pos, pos, el);
    }
    return tree;
}

TEST_CASE("rtree basic ops")
{
  // this is needed as template param in RTree
  // that's what I want to solve :)
  #define dims 2
  int length = 4;
  auto grid = make_grid(length, dims);
  int size_init = pow(length, dims);

  RTree<Point> tree = grid_to_rtree(grid);

  REQUIRE(tree.Count() == pow(length, dims));

  std::vector<Point> found;
  
  Callback cb = [&found](Point node, const double* low, const double* high) {
    found.push_back(node);
    return true;
  };
  
  double low[dims] = {1,1};
  double high[dims] = {1,3};
  
  tree.Search(low, high, cb);
  REQUIRE(found.size() == 3);
  
  std::vector<Point> expected = {{1,1}, {1,2}, {1,3}};
  REQUIRE_THAT(found, Catch::Matchers::UnorderedEquals(expected));

  int remove_cb_calls = 0;
  Callback remove_cb = [&remove_cb_calls](Point node, const double *low, const double *high) {
    remove_cb_calls++;
    return true;
  };
  { 
    double low[dims] = {1,1};
    double high[dims] = {2,2};
    int remove_count = tree.Remove(low, high, remove_cb);
    REQUIRE(remove_count == 4);
  }
  REQUIRE(remove_cb_calls == 4);

  found.clear();
  tree.Search(low, high, cb);
  REQUIRE(found.size() == 1);

  expected = {{1,3}};
  REQUIRE_THAT(found, Catch::Matchers::UnorderedEquals(expected));
  REQUIRE(tree.Count() == size_init - 4);
  
}


// FIXME
TEST_CASE("bounds after removals & insertions")
{
  int length = 10;
  auto grid = make_grid(length, dims);
  // cout << pp(grid.points) << endl;

  RTree<Point> tree = grid_to_rtree(grid);

  Point expected_min = {0,0};
  Point expected_max = {9,9};

  Point actual_min;
  Point actual_max;
  auto bounds = tree.Bounds();
  actual_min = bounds.m_min;
  actual_max = bounds.m_max;

  REQUIRE(actual_min == expected_min);
  REQUIRE(actual_max == expected_max);

  Callback remove_cb = [](Point node, const double *low, const double *high) {
    return true;
  };
  double low[dims] = {0, 0};
  double high[dims] = {2, 2};

  tree.Remove(low, high, remove_cb);

  bounds = tree.Bounds();
  actual_min = bounds.m_min;
  actual_max = bounds.m_max;

  expected_min = {3,3};
  // expected_max = {9,9}; // same

  // FIXME
  // REQUIRE(actual_min == expected_min);
  // REQUIRE(actual_max == expected_max);
}
