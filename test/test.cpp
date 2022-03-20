#include <drtree/drtree.hpp>
#include <array>
#include <catch2/catch.hpp>
#include <chrono> // for high_resolution_clock
#include <easyprint.hpp>
#include <iostream>
#include <drtree/RTree.h>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;

using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

auto now = high_resolution_clock::now;

#define duration_ms(a)                                                  \
  std::chrono::duration_cast<std::chrono::milliseconds>(a).count()
#define duration_ns(a)                                                  \
  std::chrono::duration_cast<std::chrono::nanoseconds>(a).count()

template <typename T> std::string pp(const T &x) {
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
Grid make_grid(int length, int dims = 2) {
  Grid grid;
  grid.dims = dims;
  int size = pow(length, dims);
  grid.points.resize(size);
  for (int i = 0; i < size; i++) {
    auto &el = grid.points[i];
    // not optimal, but makes things simple (no template, runtime only)
    el.resize(dims);
    for (int dim = 0; dim < dims; dim++) {
      int denum = pow(length, dim);
      int x = (int)(i / denum) % length;
      int mod = x % length;
      el[dims - dim - 1] = mod;
    }
  }

  return grid;
}

drtree<Point> grid_to_rtree(const Grid &grid) {
  drtree<Point, double> tree(grid.dims);
  for (auto &el : grid.points) {
    const double *pos = &el[0];
    tree.Insert(pos, pos, el);
  }
  return tree;
}

template <int DIMS>
RTree<Point, double, DIMS> grid_to_rtree_template(const Grid &grid, int dims = 0) {
  if(dims == 0) dims = grid.dims;
  RTree<Point, double, DIMS> tree(dims);
  for (auto &el : grid.points) {
    const double *pos = &el[0];
    tree.Insert(pos, pos, el);
  }
  return tree;
}

TEST_CASE("basic operations") {
  int dims = 2;
  int length = 4;
  auto grid = make_grid(length, dims);
  int size_init = pow(length, dims);

  drtree<Point> tree = grid_to_rtree(grid);

  REQUIRE(tree.Count() == pow(length, dims));

  std::vector<Point> found;

  Callback cb = [&found](Point node, const double *low, const double *high) {
    found.push_back(node);
    return true;
  };

  double low[2] = {1, 1};
  double high[2] = {1, 3};

  tree.Search(low, high, cb);
  REQUIRE(found.size() == 3);

  std::vector<Point> expected = {{1, 1}, {1, 2}, {1, 3}};
  REQUIRE_THAT(found, Catch::Matchers::UnorderedEquals(expected));

  int remove_cb_calls = 0;
  Callback remove_cb = [&remove_cb_calls](Point node, const double *low,
                                          const double *high) {
    remove_cb_calls++;
    return true;
  };
  {
    double low[2] = {1, 1};
    double high[2] = {2, 2};
    int remove_count = tree.Remove(low, high, remove_cb);
    REQUIRE(remove_count == 4);
  }
  REQUIRE(remove_cb_calls == 4);

  found.clear();
  tree.Search(low, high, cb);
  REQUIRE(found.size() == 1);

  expected = {{1, 3}};
  REQUIRE_THAT(found, Catch::Matchers::UnorderedEquals(expected));
  REQUIRE(tree.Count() == size_init - 4);
}

// FIXME
TEST_CASE("bounds after removals & insertions") {
  int length = 10;
  auto grid = make_grid(length, 2);
  // cout << pp(grid.points) << endl;

  drtree<Point> tree = grid_to_rtree(grid);

  Point expected_min = {0, 0};
  Point expected_max = {9, 9};

  auto bounds = tree.Bounds();

  REQUIRE(bounds.low == expected_min);
  REQUIRE(bounds.high == expected_max);

  Callback remove_cb = [](Point node, const double *low, const double *high) {
    return true;
  };
  double low[2] = {0, 0};
  double high[2] = {2, 2};

  tree.Remove(low, high, remove_cb);

  bounds = tree.Bounds();
  // actual_min.assign(bounds.m_min, bounds.m_min + 2);
  // actual_max.assign(bounds.m_max, bounds.m_max + 2);

  expected_min = {3, 3};
  // expected_max = {9,9}; // same

  // FIXME
  // REQUIRE(actual_min == expected_min);
  // REQUIRE(actual_max == expected_max);
}

// total heap usage: 83,986 allocs, 83,986 frees, 5,511,972 bytes allocated
TEST_CASE("10x4d drtree", "[benchmark][drtree]") {
  auto grid = make_grid(10, 4); // 10k points

  auto t1 = now();
  // BENCHMARK("init") {
  drtree<Point> tree = grid_to_rtree(grid);
  // };
  auto t2 = now();
  // CAPTURE(duration_ms(t2 - t1));
  WARN("init took " << duration_ms(t2 - t1) << "ms");

  double low[4] = {1, 3, 2, 4};
  double high[4] = {2, 3, 4, 4};
  std::vector<Point> found;
  Callback cb = [&](Point node, const double *low, const double *high) {
    found.push_back(node);
    return true;
  };
  t1 = now();
  for (int i = 0; i < 1000; i++) {
    found.clear();
    tree.Search(low, high, cb);
  }
  t2 = now();
  WARN("search x1000 took " << duration_ms(t2 - t1) << " ms ");
  REQUIRE(found.size() == 6);
  std::vector<Point> expected = {{1.0, 3.0, 3.0, 4.0}, {1.0, 3.0, 4.0, 4.0},
                                 {1.0, 3.0, 2.0, 4.0}, {2.0, 3.0, 2.0, 4.0},
                                 {2.0, 3.0, 4.0, 4.0}, {2.0, 3.0, 3.0, 4.0}};
  REQUIRE_THAT(found, Catch::Matchers::UnorderedEquals(expected));
}

// 4d normal: total heap usage: 107,771 allocs, 107,771 frees, 5,694,676 bytes allocated
// 8d MAXDIMS: total heap usage: 107,771 allocs, 107,771 frees, 6,813,908 bytes allocated
TEST_CASE("10x4d rtree", "[benchmark][rtree]") {
#define MAXDIMS 4
  int dims = 4;
  auto grid = make_grid(10, dims);

  auto t1 = now();
  RTree<Point, double, MAXDIMS> tree = grid_to_rtree_template<MAXDIMS>(grid);
  auto t2 = now();
  WARN("init took " << duration_ms(t2 - t1) << "ms");

#undef MAXDIMS
  
  double low[4] = {1, 3, 2, 4};
  double high[4] = {2, 3, 4, 4};
  std::vector<Point> found;
  Callback cb = [&](Point node, const double *low, const double *high) {
    found.push_back(node);
    return true;
  };
  t1 = now();
  for(int i=0; i< 1000; i++) {
    found.clear();
    tree.Search(low, high, cb);
  }
  t2 = now();
  WARN("search x1000 took " << duration_ms(t2 - t1) << " ms ");
  REQUIRE(found.size() == 6);
  std::vector<Point> expected = {{1.0, 3.0, 3.0, 4.0}, {1.0, 3.0, 4.0, 4.0},
                                 {1.0, 3.0, 2.0, 4.0}, {2.0, 3.0, 2.0, 4.0},
                                 {2.0, 3.0, 4.0, 4.0}, {2.0, 3.0, 3.0, 4.0}};
  REQUIRE_THAT(found, Catch::Matchers::UnorderedEquals(expected));
  
}

// total heap usage: 119,490 allocs, 119,490 frees, 10,676,298 bytes allocated
// RECT_ONE_ARRAY: total heap usage: 119,490 allocs, 119,490 frees, 10,073,986 bytes allocated
TEST_CASE("200x2d drtree", "[benchmark][drtree]") {
  auto grid = make_grid(200, 2); // 100x100=>10k points
  // cout << pp(grid.points) << endl;

  auto t1 = now();
  drtree<Point> tree = grid_to_rtree(grid);
  auto t2 = now();
  WARN("init took " << duration_ms(t2 - t1) << "ms");

  // return;
  double low[2] = {5, 2};
  double high[2] = {6, 4};
  std::vector<Point> found;
  Callback cb = [&](Point node, const double *low, const double *high) {
    found.push_back(node);
    return true;
  };
  t1 = now();
  for (int i = 0; i < 1000; i++) {
    found.clear();
    tree.Search(low, high, cb);
  }
  t2 = now();
  WARN("search x1000 took " << duration_ms(t2 - t1) << " ms ");
  REQUIRE(found.size() == 6);
  std::vector<Point> expected = {
    { 5.0, 2.0 },
    { 5.0, 3.0 },
    { 5.0, 4.0 },
    { 6.0, 2.0 },
    { 6.0, 3.0 },
    { 6.0, 4.0 },
  };
  REQUIRE_THAT(found, Catch::Matchers::UnorderedEquals(expected));
}

// total heap usage: 218,693 allocs, 218,693 frees, 9,627,530 bytes allocated
// 4 MAXDIMS: total heap usage: 218,693 allocs, 218,693 frees, 12,036,490 bytes allocated
TEST_CASE("200x2d rtree", "[benchmark][rtree]") {
#define MAXDIMS 2
  int dims = 2;
  
  auto grid = make_grid(200, dims);

  auto t1 = now();
  RTree<Point, double, MAXDIMS> tree = grid_to_rtree_template<MAXDIMS>(grid);
  auto t2 = now();
#undef MAXDIMS
  
  WARN("init took " << duration_ms(t2 - t1) << "ms");
  // return;

  double low[2] = {5, 2};
  double high[2] = {6, 4};
  std::vector<Point> found;
  Callback cb = [&](Point node, const double *low, const double *high) {
    found.push_back(node);
    return true;
  };
  t1 = now();
  for (int i = 0; i < 1000; i++) {
    found.clear();
    tree.Search(low, high, cb);
  }
  t2 = now();
  WARN("search x1000 took " << duration_ms(t2 - t1) << " ms ");
  REQUIRE(found.size() == 6);
  std::vector<Point> expected = {
    { 5.0, 2.0 },
    { 5.0, 3.0 },
    { 5.0, 4.0 },
    { 6.0, 2.0 },
    { 6.0, 3.0 },
    { 6.0, 4.0 },
  };
  REQUIRE_THAT(found, Catch::Matchers::UnorderedEquals(expected));
}

TEST_CASE("allocations", "[.temp]") {
  auto grid = make_grid(2, 2);
  // cout << pp(grid.points) << endl;

  auto t1 = high_resolution_clock::now();
  drtree<Point> tree = grid_to_rtree(grid);
  auto t2 = high_resolution_clock::now();
  cout << "init took " << duration_ms(t2 - t1) << "ms" << endl;
}

TEST_CASE("allocations template", "[.temp]") {
  auto grid = make_grid(2, 2);
  // cout << pp(grid.points) << endl;

  auto t1 = high_resolution_clock::now();
  RTree<Point, double, 2> tree = grid_to_rtree_template<2>(grid);
  auto t2 = high_resolution_clock::now();
  cout << "init templ took " << duration_ms(t2 - t1) << "ms" << endl;
}
