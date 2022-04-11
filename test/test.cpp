#include <drtree/drtree.hpp>
#include <drtree/drtree2.hpp>
#include <drtree/drtree3.hpp>
#include <array>
#include <catch2/catch.hpp>
#include <chrono> // for high_resolution_clock
#include <easyprint.hpp>
#include <iostream>
#include <fstream>
#include <drtree/RTree.h>
#include <aod/rtree.hpp>
#include <vector>
#include <random>


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

void shuffle_deterministic(Grid& grid) {
  std::mt19937 gen(0);
  std::shuffle(std::begin(grid.points), std::end(grid.points), gen);
}

// for debugging, creating the "expected" vectors
void sort_points(std::vector<Point> &points) {
  std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) {
    if( a[0] == b[0] ) {
      return a[1] < b[1];
    } else {
      return a[0] < b[0];
    }
  });
}

/**
 * see https://stackoverflow.com/a/49514906/8720686
 * Given length=3,dims=2 it gives
 * {{0, 0}, {0, 1}, {0, 2},
 *  {1, 0}, {1, 1}, {1, 2},
 *  {2, 0}, {2, 1}, {2, 2}}
 */
Grid make_grid(int length, int dims = 2, bool shuffle = false) {
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
  if(shuffle) {
    shuffle_deterministic(grid);
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

aod::rtree<Point> grid_to_aod_rtree(const Grid &grid) {
  aod::rtree<Point> tree(grid.dims);
  for (auto &el : grid.points) {
    // const double *pos = &el[0];
    tree.insert(el, el, el);
  }
  return tree;
}

drtree2<Point> grid_to_drtree2(const Grid &grid) {
  drtree2<Point, double> tree(grid.dims);
  for (auto &el : grid.points) {
    const double *pos = &el[0];
    tree.Insert(pos, pos, el);
  }
  return tree;
}

drtree3<Point> grid_to_drtree3(const Grid &grid) {
  drtree3<Point> tree(grid.dims);
  for (auto &el : grid.points) {
    const double *pos = &el[0];
    tree.push(pos, pos, el);
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

  // REQUIRE(actual_min == expected_min);
  // REQUIRE(actual_max == expected_max);
}

// total heap usage: 83,986 allocs, 83,986 frees, 5,511,972 bytes allocated
TEST_CASE("10x4d drtree", "[.skip]") {
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
TEST_CASE("10x4d rtree", "[.skip]") {
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

const int benchmark_size = 200;
const int searches = 10000;
const char* search_str = "search x10000 ";
// total heap usage: 119,490 allocs, 119,490 frees, 10,676,298 bytes allocated
// RECT_ONE_ARRAY: total heap usage: 119,490 allocs, 119,490 frees, 10,073,986 bytes allocated
TEST_CASE("200x2d drtree", "[benchmark][drtree]") {
  auto grid = make_grid(benchmark_size, 2); // 100x100=>10k points
  // cout << pp(grid.points) << endl;

  auto t1 = now();
  drtree<Point> tree = grid_to_rtree(grid);
  auto t2 = now();
  WARN("init took " << duration_ms(t2 - t1) << "ms");
  REQUIRE(tree.Count() == benchmark_size*benchmark_size);

  // return;
  double low[2] = {5, 2};
  double high[2] = {6, 4};
  std::vector<Point> found;
  Callback cb = [&](Point node, const double *low, const double *high) {
    found.push_back(node);
    return true;
  };
  t1 = now();
  for (int i = 0; i < searches; i++) {
    found.clear();
    tree.Search(low, high, cb);
  }
  t2 = now();
  WARN(search_str << duration_ms(t2 - t1) << " ms ");
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

  t1 = now();
  drtree<Point> tree2 = tree;
  t2 = now();
  WARN("copy took " << duration_ms(t2 - t1) << " ms ");
  REQUIRE(tree2.Count() == tree.Count());

  found.clear();
  tree2.Search(low, high, cb);
  REQUIRE_THAT(found, Catch::Matchers::UnorderedEquals(expected));
}

// total heap usage: 218,693 allocs, 218,693 frees, 9,627,530 bytes allocated
// 4 MAXDIMS: total heap usage: 218,693 allocs, 218,693 frees, 12,036,490 bytes allocated
TEST_CASE("200x2d rtree", "[benchmark][rtree]") {
#define MAXDIMS 2
  int dims = 2;
  auto grid = make_grid(benchmark_size, dims);

  auto t1 = now();
  RTree<Point, double, MAXDIMS> tree = grid_to_rtree_template<MAXDIMS>(grid);
  auto t2 = now();
  
  WARN("init size " << benchmark_size << " took " << duration_ms(t2 - t1) << "ms");
  // return;

  double low[2] = {5, 2};
  double high[2] = {6, 4};
  std::vector<Point> found;
  Callback cb = [&](Point node, const double *low, const double *high) {
    found.push_back(node);
    return true;
  };
  t1 = now();
  for (int i = 0; i < searches; i++) {
    found.clear();
    tree.Search(low, high, cb);
  }
  t2 = now();
  WARN(search_str << duration_ms(t2 - t1) << " ms ");
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

  t1 = now();
  RTree<Point, double, MAXDIMS> tree2 = tree;
  t2 = now();
  WARN("copy took " << duration_ms(t2 - t1) << " ms ");
  REQUIRE(tree2.Count() == tree.Count());

  #undef MAXDIMS
}

TEST_CASE("200x2d drtree3", "[drtree3][benchmark]") {
  auto grid = make_grid(benchmark_size); // TODO crashes with 3 (works with 2)
  auto t1 = high_resolution_clock::now();
  drtree3<Point> tree = grid_to_drtree3(grid);
  auto t2 = high_resolution_clock::now();
  WARN("init size " << benchmark_size << " took " << duration_ms(t2 - t1) << "ms");
  REQUIRE(tree.size() == benchmark_size*benchmark_size);

  std::vector<Point> expected = {
    { 5.0, 2.0 },
    { 5.0, 3.0 },
    { 5.0, 4.0 },
    { 6.0, 2.0 },
    { 6.0, 3.0 },
    { 6.0, 4.0 },
  };
  std::vector<double> low = {5 , 2};
  std::vector<double> high = {6 , 4};
  auto found = tree.search(low, high);
  REQUIRE_THAT(found, Catch::Matchers::UnorderedEquals(expected));

  t1 = now();
  for (int i = 0; i < searches; i++) {
    tree.search(low, high, found);
  }
  t2 = now();
  WARN(search_str << duration_ms(t2 - t1) << " ms ");

  t1 = now();
  drtree3<Point>tree2 = tree;
  t2 = now();
  WARN("copy took " << duration_ms(t2 - t1) << " ms ");

  REQUIRE(tree2.size() == tree.size());
  auto found2 = tree2.search(low, high);
  REQUIRE_THAT(found2, Catch::Matchers::UnorderedEquals(expected));
}

TEST_CASE("200x2d aod::rtree", "[aod::rtree][benchmark]") {
  int size = benchmark_size;
  auto grid = make_grid(size, 2); // TODO crashes with 3 (works with 2)
  auto t1 = high_resolution_clock::now();
  aod::rtree<Point> tree = grid_to_aod_rtree(grid);
  auto t2 = high_resolution_clock::now();
  WARN("init size " << size << " took " << duration_ms(t2 - t1) << "ms");
  REQUIRE(tree.size() == size*size);

  if (false) {
    std::ofstream ofs("aod-rtree-200x200-init.xml", std::ofstream::out);
    ofs << tree.to_xml();
    ofs.close();
  }

  std::vector<Point> expected = {
    { 5.0, 2.0 },
    { 5.0, 3.0 },
    { 5.0, 4.0 },
    { 6.0, 2.0 },
    { 6.0, 3.0 },
    { 6.0, 4.0 },
  };
  std::vector<double> low = {5 , 2};
  std::vector<double> high = {6 , 4};
  auto found = tree.search(low, high);
  REQUIRE_THAT(found, Catch::Matchers::UnorderedEquals(expected));

  t1 = now();
  for (int i = 0; i < searches; i++) {
    tree.search(low, high, found);
  }
  t2 = now();
  WARN(search_str << duration_ms(t2 - t1) << " ms ");

  t1 = now();
  aod::rtree<Point>tree2 = tree;
  t2 = now();
  WARN("copy took " << duration_ms(t2 - t1) << " ms ");

  REQUIRE(tree2.size() == tree.size());
  auto found2 = tree2.search(low, high);
  REQUIRE_THAT(found2, Catch::Matchers::UnorderedEquals(expected));

  // removing fat cross
  t1 = now();
  // removing a fat cross in the center, leaving corners with side=2
  int retain_side = 2;
  int expected_removed1 = size * (size-2*retain_side);
  int removed1 = tree.remove({0,retain_side},{size-1,size-1-retain_side});
  REQUIRE(removed1 == expected_removed1);
  int expected_removed2 = 2 * (size-2*retain_side)*retain_side;
  int removed2 = tree.remove({retain_side, 0}, {size-1-retain_side, size-1});
  REQUIRE(removed2 == expected_removed2);
  REQUIRE(tree.size() == size*size - removed1 - removed2); // 4 left in each corner (2x2)
  REQUIRE(tree.size() == 16);
  found = tree.search({0, 0}, {size, size});
  REQUIRE(found.size() == 16);
  t2 = now();
  WARN("removing fat cross took " << duration_ms(t2 - t1) << " ms ");
}

TEST_CASE("drtree3 test: 5x5", "[drtree3][basic test]") {
  const int size = 5;
  auto grid = make_grid(size);
  shuffle_deterministic(grid);

  drtree3<Point> tree = grid_to_drtree3(grid);
  REQUIRE(tree.size() == size*size);
  std::vector<Point> expected;
  std::vector<Point> found;
  int removed;

  expected = {
    { 3.0, 2.0 },
    { 3.0, 3.0 },
    { 3.0, 4.0 },
  };

  // leaving only two strips: x=0 and y=4
  removed = tree.remove({1,0}, {4, 3});
  REQUIRE(removed == 16);

  found = tree.search({0,0}, {4, 4});
  REQUIRE(found.size() == 9);
  expected = {
    {0, 0}, {0, 1} , {0, 2}, {0, 3}, {0, 4},
    {1, 4}, {2, 4}, {3, 4}, {4, 4}
  };
  REQUIRE_THAT(found, Catch::Matchers::UnorderedEquals(expected));
}

// getting duplicates with MAXNODES 8
// with MAXNODES 4 it works
TEST_CASE("drtree3 test: 10x10", "[drtree3][basic test][.skip-fix]") {
  const int size = 10;
  auto grid = make_grid(size);
  shuffle_deterministic(grid);
  
  drtree3<Point> tree = grid_to_drtree3(grid);
  REQUIRE(tree.size() == size*size);
  std::vector<Point> expected;
  std::vector<Point> found;
  
  expected = {
    { 5.0, 5.0 }, { 5.0, 6.0 }, { 5.0, 7.0 },
    { 6.0, 5.0 }, { 6.0, 6.0 }, { 6.0, 7.0 }};
  found = tree.search({5, 5}, {6, 7});
  REQUIRE_THAT(found, Catch::Matchers::UnorderedEquals(expected));

  tree.remove({1,1},{8,8});
  found = tree.search({0, 0}, {10, 10});
  sort_points(found);
  expected =  {{0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}, {0, 6}, {0, 7}, {0, 8}, {0, 9},
               {1, 0}, {1, 9},
               {2, 0}, {2, 9},
               {3, 0}, {3, 9},
               {4, 0}, {4, 9},
               {5, 0}, {5, 9},
               {6, 0}, {6, 9},
               {7, 0}, {7, 9},
               {8, 0}, {8, 9},
               {9, 0}, {9, 1}, {9, 2}, {9, 3}, {9, 4}, {9, 5}, {9, 6}, {9, 7}, {9, 8}, {9, 9}};
  REQUIRE_THAT(found, Catch::Matchers::UnorderedEquals(expected));
}

TEST_CASE("drtree3 test: 4x4", "[drtree3][basic test]") {
  const int size = 4;
  auto grid = make_grid(size);
  shuffle_deterministic(grid);
  
  drtree3<Point> tree = grid_to_drtree3(grid);
  REQUIRE(tree.size() == size*size);
  std::vector<Point> expected;
  std::vector<Point> found;
  int removed;

  // removing a fat cross in the center
  // leaving corners with side = 1
  removed = tree.remove({0,1},{size,size-1}); // 12
  REQUIRE(removed == 8);
  removed = tree.remove({1,0},{size-1,size}); // vertical
  // 2 horizontal stripes removed 2 in bottom row & 2 in top row
  REQUIRE(removed == 2 + 2);
  REQUIRE(tree.size() == 1*4); // 1 left in each cornet


  found = tree.search({0, 0}, {10, 10});
  sort_points(found);
  expected =  { { 0.0, 0.0 }, { 0.0, 3.0 }, { 3.0, 0.0 }, { 3.0, 3.0 } };
  REQUIRE_THAT(found, Catch::Matchers::UnorderedEquals(expected));
}

TEST_CASE("aod rtree test: 4x4", "[aod::rtree][basic test]") {
  const int size = 4;
  auto grid = make_grid(size);
  shuffle_deterministic(grid);
  
  aod::rtree<Point> tree = grid_to_aod_rtree(grid);
  REQUIRE(tree.size() == size*size);
  std::vector<Point> expected;
  std::vector<Point> found;
  // int removed;

  // cout << "-- tree " << endl << tree.to_xml() << endl;
  
  std::ofstream ofs("aod-rtree-4x4.xml", std::ofstream::out);
  ofs << tree.to_xml();
  ofs.close();

  found = tree.search({1,1}, {2, 2});
  sort_points(found);
  expected = {{1, 1}, {1, 2}, {2, 1}, {2, 2}};
  REQUIRE_THAT(found, Catch::Matchers::UnorderedEquals(expected));
}

TEST_CASE("drtree3 test: 7x7", "[drtree3][basic test][.skip-fix]") {
  const int size = 7;
  auto grid = make_grid(size);
  shuffle_deterministic(grid);
  
  drtree3<Point> tree = grid_to_drtree3(grid);
  cout << "INIT:" << endl << tree.to_string() << endl;

  REQUIRE(tree.size() == size*size);
  std::vector<Point> expected;
  std::vector<Point> found;
  int removed;

  // removing a fat cross in the center, leaving corners with side=2
  int retain_side = 2;
  // horizontal: 7 * (7-2*2) = 21
  removed = tree.remove({0,retain_side},{size-1,size-1-retain_side});
  REQUIRE(removed == 21);
  REQUIRE(tree.size() == size*size-removed); // 49-21=28
  found = tree.search({0, 0}, {10, 10});
  cout << tree.to_string() << endl;

  sort_points(found);
  cout << "found " << pp(found) << endl;
  REQUIRE(found.size() == tree.size()); // same node appears as child 3 times
  return;

  removed = tree.remove({retain_side, 0}, {size-1-retain_side, size-1}); // vertical
  cout << tree.to_string() << endl;
  // 2 horizontal stripes removed 2 in bottom row (height 2) & 2 in top row
  REQUIRE(removed == 2 * (size-retain_side*2)*retain_side); // 12
  REQUIRE(tree.size() == 16); // 4 left in each corner (2x2)

  found = tree.search({0, 0}, {10, 10});
  cout << tree.to_string() << endl;
  
  sort_points(found);
  expected = {  { 0.0, 0.0 }, { 0.0, 1.0 }, { 1.0, 0.0 }, { 1.0, 1.0 },
                { 0.0, 5.0 }, { 0.0, 6.0 }, { 1.0, 5.0 }, { 1.0, 6.0 },
                { 5.0, 0.0 }, { 6.0, 0.0 }, { 5.0, 1.0 }, { 6.0, 1.0 },
                { 5.0, 5.0}, { 5.0, 6.0 }, { 6.0, 5.0 }, { 6.0, 6.0 } };
  REQUIRE_THAT(found, Catch::Matchers::UnorderedEquals(expected));
}

TEST_CASE("drtree3 test: 8x8", "[drtree3][basic test]") {
  const int size = 8;
  auto grid = make_grid(size);
  shuffle_deterministic(grid);
  
  drtree3<Point> tree = grid_to_drtree3(grid);
  REQUIRE(tree.size() == size*size);
  std::vector<Point> expected;
  std::vector<Point> found;
  int removed;

  // removing a fat cross in the center, leaving corners with side=2
  int retain_side = 2;
  // horizontal: 8 * (8-2*retain_side) = 24
  removed = tree.remove({0,retain_side},{size-1,size-1-retain_side});
  REQUIRE(removed == 24);
  REQUIRE(tree.size() == size*size - 24);

  removed = tree.remove({retain_side, 0}, {size-1-retain_side, size-1}); // vertical
  // 2 horizontal stripes removed 2 rowsof 3 in bottom row & same in top
  REQUIRE(removed == 2 * (size-2*retain_side)*retain_side); // 8
  REQUIRE(tree.size() == 4*4); // 4 left in each corner (2x2)

  found = tree.search({0, 0}, {10, 10});
  sort_points(found);
  expected = {
    { 0.0, 0.0 }, { 0.0, 1.0 }, { 1.0, 0.0 }, { 1.0, 1.0 },
    { 0.0, 4.0 }, { 0.0, 5.0 }, { 1.0, 4.0 }, { 1.0, 5.0 },
    { 4.0, 0.0 }, { 4.0, 1.0 }, { 5.0, 0.0 }, { 5.0, 1.0 },
    { 4.0, 4.0 }, { 4.0, 5.0 }, { 5.0, 4.0 }, { 5.0, 5.0 } };
  REQUIRE_THAT(found, Catch::Matchers::UnorderedEquals(expected));
}

TEST_CASE("aod::rtree 8x8", "[aod::rtree]") {
  const int size = 8;
  auto grid = make_grid(size);
  shuffle_deterministic(grid);
  
  aod::rtree<Point> tree = grid_to_aod_rtree(grid);
  REQUIRE(tree.size() == size*size);

  std::ofstream ofs;
  
  ofs = std::ofstream("aod-rtree-8x8-init.xml", std::ofstream::out);
  ofs << tree.to_xml();
  ofs.close();

  std::vector<Point> expected;
  std::vector<Point> found;
  int removed;

  // removing a fat cross in the center, leaving corners with side=2
  int retain_side = 2;

  found = tree.search({0,retain_side},{size-1,size-1-retain_side});
  REQUIRE(found.size() == 32);
  // horizontal: 8 * (8-2*retain_side) = 32
  // {0, 2}, {7, 5}
  int removed1 = tree.remove({0,retain_side},{size-1,size-1-retain_side});
  REQUIRE(removed1 == 32);
  REQUIRE(tree.size() == size*size - removed1);
  found = tree.search({0, 0}, {10, 10});
  REQUIRE(found.size() == size*size - removed1);

  ofs = std::ofstream("aod-rtree-8x8-removed-1.xml", std::ofstream::out);
  ofs << tree.to_xml();
  ofs.close();

  int removed2 = tree.remove({retain_side, 0}, {size-1-retain_side, size-1}); // vertical
  // 2 horizontal stripes removed 2 rowsof 3 in bottom row & same in top
  // 2 * (size-2*retain_side)*retain_side = 2*4*2 = 16
  REQUIRE(removed2 == 16);
  REQUIRE(tree.size() == size*size - removed1 - removed2); // 4 left in each corner (2x2)
  REQUIRE(tree.size() == 16); // 4 left in each corner (2x2)

  ofs = std::ofstream("aod-rtree-8x8-removed.xml", std::ofstream::out);
  ofs << tree.to_xml();
  ofs.close();

  found = tree.search({0, 0}, {10, 10});
  sort_points(found);
  expected = {
    { 0.0, 0.0 }, { 0.0, 1.0 }, { 1.0, 0.0 }, { 1.0, 1.0 },
    { 0.0, 6.0 }, { 0.0, 7.0 }, { 1.0, 6.0 }, { 1.0, 7.0 },
    { 6.0, 0.0 }, { 6.0, 1.0 }, { 7.0, 0.0 }, { 7.0, 1.0 },
    { 6.0, 6.0 }, { 6.0, 7.0 }, { 7.0, 6.0 }, { 7.0, 7.0 } };
  REQUIRE_THAT(found, Catch::Matchers::UnorderedEquals(expected));
}

TEST_CASE("aod::rtree 10x10", "[aod::rtree]") {
  const int size = 10;
  auto grid = make_grid(size);
  shuffle_deterministic(grid);
  
  aod::rtree<Point> tree = grid_to_aod_rtree(grid);
  REQUIRE(tree.size() == size*size);

  std::ofstream ofs;
  
  ofs = std::ofstream("aod-rtree-10x10-init.xml", std::ofstream::out);
  ofs << tree.to_xml();
  ofs.close();
  
  std::vector<Point> expected;
  std::vector<Point> found;
  int removed;

  // removing a fat cross in the center, leaving corners with side=2
  int retain_side = 2;

  found = tree.search({0,retain_side},{size-1,size-1-retain_side});
  REQUIRE(found.size() == 60);
  // horizontal: 10 * (10-2*retain_side) = 10*(10-2*2) = 60
  int removed1 = tree.remove({0,retain_side},{size-1,size-1-retain_side});
  REQUIRE(removed1 == 60);
  ofs = std::ofstream("aod-rtree-10x10-removed-1.xml", std::ofstream::out);
  ofs << tree.to_xml();
  ofs.close();
  REQUIRE(tree.size() == size*size - removed1);
  found = tree.search({0, 0}, {10, 10});
  REQUIRE(found.size() == size*size - removed1);

  // return;

  int removed2 = tree.remove({retain_side, 0}, {size-1-retain_side, size-1}, true); // vertical
  ofs = std::ofstream("aod-rtree-10x10-removed-2.xml", std::ofstream::out);
  ofs << tree.to_xml();
  ofs.close();
  // 2 horizontal stripes removed 2 rowsof n in bottom row & same in top
  // 2 * (size-2*retain_side)*retain_side = 2*(10-2*2)*2 = 2*6*2
  REQUIRE(removed2 == 24);
  REQUIRE(tree.size() == size*size - removed1 - removed2); // 4 left in each corner (2x2)
  REQUIRE(tree.size() == 16); // 4 left in each corner (2x2)



  found = tree.search({0, 0}, {10, 10});
  sort_points(found);
  expected = {
    { 0.0, 0.0 }, { 0.0, 1.0 }, { 1.0, 0.0 }, { 1.0, 1.0 },
    { 0.0, 8.0 }, { 0.0, 9.0 }, { 1.0, 8.0 }, { 1.0, 9.0 },
    { 8.0, 0.0 }, { 8.0, 1.0 }, { 9.0, 0.0 }, { 9.0, 1.0 },
    { 8.0, 8.0 }, { 8.0, 9.0 }, { 9.0, 8.0 }, { 9.0, 9.0 } };
  REQUIRE_THAT(found, Catch::Matchers::UnorderedEquals(expected));
}

TEST_CASE("aod::rtree temp", "[aod::rtree][.skip]") {
  const int size = 200;
  auto grid = make_grid(size);
  shuffle_deterministic(grid);
  
  aod::rtree<Point> tree = grid_to_aod_rtree(grid);
  REQUIRE(tree.size() == size*size);

  std::ofstream ofs;
  
  ofs = std::ofstream("aod-rtree-temp-init.xml", std::ofstream::out);
  ofs << tree.to_xml();
  ofs.close();
  
  std::vector<Point> expected;
  std::vector<Point> found;
  int removed;

  // removing a fat cross in the center, leaving corners with side=2
  int retain_side = 2;

  // horizontal: size * (size-2*retain_side) = 13*(13-4)
  int expected_removed1 = size * (size-2*retain_side);
  found = tree.search({0,retain_side},{size-1,size-1-retain_side});
  REQUIRE(found.size() == expected_removed1);

  int removed1 = tree.remove({0,retain_side},{size-1,size-1-retain_side});
  ofs = std::ofstream("aod-rtree-temp-removed-1.xml", std::ofstream::out);
  ofs << tree.to_xml();
  ofs.close();
  REQUIRE(removed1 == expected_removed1);
  REQUIRE(tree.size() == size*size - removed1);
  found = tree.search({0, 0}, {size, size});
  REQUIRE(found.size() == size*size - removed1);

  // return;
  //  // vertical
  // 2 horizontal stripes removed 2 rowsof 3 in bottom row & same in top
  // 2 * (size-2*retain_side)*retain_side = 2*(13-2*2)*2 = 36
  int expected_removed2 = 2 * (size-2*retain_side)*retain_side;
  int removed2 = tree.remove({retain_side, 0}, {size-1-retain_side, size-1});
  REQUIRE(removed2 == expected_removed2);
  REQUIRE(tree.size() == size*size - removed1 - removed2); // 4 left in each corner (2x2)
  REQUIRE(tree.size() == 16); // 4 left in each corner (2x2)

  ofs = std::ofstream("aod-rtree-temp-removed-2.xml", std::ofstream::out);
  ofs << tree.to_xml();
  ofs.close();

  found = tree.search({0, 0}, {size, size});
  sort_points(found);
  expected = {
    { 0.0, 0.0 }, { 0.0, 1.0 }, { 1.0, 0.0 }, { 1.0, 1.0 },
    { 0.0, size -2 }, { 0.0, size-1 }, { 1.0, size-2 }, { 1.0, size-1 },
    { size-2, 0.0 }, { size-2, 1.0 }, { size-1, 0.0 }, { size-1, 1.0 },
    { size-2, size-2 }, { size-2, size-1 }, { size-1, size-2 }, { size-1, size-1 } };
  REQUIRE_THAT(found, Catch::Matchers::UnorderedEquals(expected));
}

TEST_CASE("drtree3 test: 100x100", "[drtree3][basic test][.skip-fix]") {
  const int size = 100;
  auto grid = make_grid(size);
  shuffle_deterministic(grid);
  
  drtree3<Point> tree = grid_to_drtree3(grid);
  REQUIRE(tree.size() == size*size);
  std::vector<Point> expected;
  std::vector<Point> found;
  int removed;

  // removing a fat cross in the center:
  // 100x96
  removed = tree.remove({0,2},{99,97}); // horizontal
  REQUIRE(removed == 100*96);
  removed = tree.remove({2,0},{97,99}); // vertical
  // 2 horizontal stripes removed 96 in bottom row & 96 in top row
  REQUIRE(removed == 96 + 96);
  REQUIRE(tree.size() == 4*4);

  // TODO fixme
  found = tree.search({0,0}, {100, 100});
  expected = {
    {0,0}, {0,1}, {1,0}, {1,1},
    {0,98}, {0,99}, {1,98}, {1,99},
    {98,0}, {98,1}, {99,0}, {99,1},
    {98,98}, {98,99}, {99,98}, {99,99},
  };
  cout  << "size " << tree.size() << " found " << pp(found) << endl;
  
  REQUIRE_THAT(found, Catch::Matchers::UnorderedEquals(expected));
}

TEST_CASE("rtree test", "[rtree][basic test]") {
  const int size = 10; // todo with size 9, searching {6,2} to {6,4}, {6,3} not returned
  auto grid = make_grid(size);
  RTree<Point, double, 2> tree = grid_to_rtree_template<2>(grid);
  REQUIRE(tree.Count() == size*size);
  std::vector<Point> expected = {
    { 6.0, 2.0 },
    { 6.0, 3.0 },
    { 6.0, 4.0 },
  };
  // TODO not getting {6,3} (with spherical volume)
  std::vector<Point> results;
  Callback cb = [&](Point node, const double *low, const double *high) {
    results.push_back(node);
    return true;
  };

  double low[] = {6,2};
  double high[] = {6,4};
  tree.Search(low, high, cb);
  REQUIRE_THAT(results, Catch::Matchers::UnorderedEquals(expected));

  int removed = tree.Remove(low, high, cb);

  REQUIRE(removed == 3);
  REQUIRE(tree.Count() == size*size - 3);

  results.clear();
  tree.Search(low, high, cb);
  REQUIRE(results.size() == 0);
  
  double low2[] = {1,1};
  double high2[] = {8,8};
  results.clear();
  return;
  removed = tree.Remove(low2, high2, cb); // crashes, heh
  REQUIRE(removed == 54);
  REQUIRE(tree.Count() == 43);


  double low3[] = {0,0};
  double high3[] = {10,10};
  results.clear();
  tree.Search(low3, high3, cb);
  // // TODO this returns 0
  REQUIRE(results.size() == 6); // [0 0] [0 1] [1 0], [9 9] [9 8] [8 9]
  expected = {};
  REQUIRE_THAT(results, Catch::Matchers::UnorderedEquals(expected));
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

TEST_CASE("drtree2 init", "[drtree2]") {
  auto grid = make_grid(10, 2);

  auto t1 = high_resolution_clock::now();
  drtree2<Point> tree = grid_to_drtree2(grid);
  auto t2 = high_resolution_clock::now();
  cout << "init drtree2 " << duration_ms(t2 - t1) << "ms" << endl;
  // REQUIRE(tree.Count() == 100);

  double low[2] = {5, 2};
  double high[2] = {6, 4};
  std::vector<Point> found;
  Callback cb = [&](Point node, const double *low, const double *high) {
    found.push_back(node);
    return true;
  };

  tree.Search(low, high, cb);
  // REQUIRE(found.size() == 2);
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
  cout << "drtree2 test done " << endl;
    

}
