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
#include <easyprint.hpp>

typedef uint32_t id_t;
struct Id {
  // std::numeric_limits<uint32_t>::max()
  // UINT32_MAX
  static const id_t nullid = std::numeric_limits<uint32_t>::max();
  id_t id = nullid; // 
  operator bool() const {
    return id != nullid;
  }
  // Id operator+(const Id& other) const {
  //   return Id{id + other.id};
  // };
  // Id operator+(const int &offset) const { return Id{id + offset}; };
  // Id& operator++() {
  //   id++;
  //   return *this;
  // }
  Id(): Id(nullid){};
  Id(const Id& other): Id(other.id){};
  Id(id_t id): id{id} {}
};

// typedef Id Nid;
// typedef Id Bid;
// typedef Id Rid;
struct Rid : Id {};

struct Nid : Id {};

struct Bid : Id {
  Bid operator+(const Bid &other) const { return Bid{id + other.id}; };
  Bid operator+(const int &offset) const { return Bid{id + offset}; };
  Bid &operator++() {
    id++;
    return *this;
  }
};

// Data id
struct Did : Id {};

struct Branch {
  Rid rect_id;
  Nid child_id;
  Did data_id;
};


std::ostream &operator<<(std::ostream &os, const Bid &bid) {
  os << "Bid{" << bid.id << "}";
  return os;
}

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
