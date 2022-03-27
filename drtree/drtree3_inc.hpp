#pragma once
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <optional>
#include <vector>
#include <algorithm>
#include <functional>
#include <vector>
#include <assert.h>
#include "./drtree3_template.hpp"

typedef uint32_t id_t;
struct Id {
  // std::numeric_limits<uint32_t>::max()
  // UINT32_MAX
  static const id_t nullid = std::numeric_limits<uint32_t>::max();
  id_t id = nullid; // 
  operator bool() const {
    return id != nullid;
  }
  Id operator+(const Id& other) const {
    return Id{id + other.id};
  };
  Id operator+(const int &offset) const { return Id{id + offset}; };
  Id& operator++() {
    id++;
    return *this;
  }
  Id(): Id(nullid){};
  Id(id_t id): id{id} {}
};

typedef Id Nid;
typedef Id Bid;
typedef Id Rid;

struct Rect {
  // Rid id;
  // Rect() : Rect(Id::nullid){};
  // Rect(Rid id) : id{id} {};
};

struct Node {
  bool is_internal() { return (level > 0); } // Not a leaf, but a internal node
  bool is_leaf() { return (level == 0); }    // A leaf, contains data

  int count = 0; ///< Count
  int level = 0; ///< Leaf is zero, others positive
  Bid branch0;
  Bid get_branch(int idx) {
    return branch0 + idx;
  }
};

struct Branch {
  // Bid id;
  Rid rect_id;
  Nid child;

  // Branch(): Branch(Id::nullid) {};
  // Branch(Bid id): id{id} {};
};
