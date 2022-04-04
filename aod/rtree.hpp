#pragma once

#include <cstdint>
#include <limits>
#include <ostream>
#include <vector>
namespace aod {

using id_t = uint32_t;
struct Id {
  static const id_t nullid = std::numeric_limits<id_t>::max();
  id_t id = nullid; // 
  operator bool() const {
    return id != nullid;
  }
  Id(): Id(nullid){};
  Id(const Id& other): Id(other.id){};
  Id(id_t id): id{id} {}
};

// Rect Id
struct Rid : Id {};

// Node Id
struct Nid : Id {};

struct Bid : Id {};

// Entry Id
struct Eid : Id {};

// Data id
struct Did : Id {};

struct Branch {
  Rid rect_id;
  Nid child_id;
  Did data_id;
};

struct Entry {
  Rid rect_id;
  Nid child_id;
  Did data_id;
};

struct Node {
  bool is_internal() { return (height > 0); } // Not a leaf, but a internal node
  bool is_leaf() { return (height == 0); }    // A leaf, contains data

  int count = 0; ///< Count
  int height = 0; ///< Leaf is zero, others positive
  Rid rect_id;
};

std::ostream &operator<<(std::ostream &os, const Bid &id) {
  os << "Bid{" << id.id << "}";
  return os;
}

std::ostream &operator<<(std::ostream &os, const Nid &id) {
  os << "Nid{" << id.id << "}";
  return os;
}


std::ostream &operator<<(std::ostream &os, const Rid &id) {
  os << "Rid{" << id.id << "}";
  return os;
}

template <class DATATYPE, class ELEMTYPE = double>
class rtree {
  static_assert(std::numeric_limits<ELEMTYPE>::is_iec559,
                "'COORD_TYPE' accepts floating-point types only");

  typedef std::vector<ELEMTYPE> Vec;
  const int MAXNODES = 8;
  const int MINNODES = MAXNODES / 2;

  std::vector<ELEMTYPE> m_rects_low;
  std::vector<ELEMTYPE> m_rects_high;

 private:
  int m_dims;
  size_t m_size;
 public:
  
  rtree() = delete;
  rtree(int dimensions);
  void insert(Vec low, Vec high, const DATATYPE &data);
  size_t size();

 private:
  // Used by insert
  Nid choose_leaf(Rid r);

  /// Ascend from leaf node L to the root, adjusting rectangles and
  /// propagating node splits as necessary
  Nid adjust_tree(Nid n, Nid nn);

  /// Given an R-tree whose root node is T, find the leaf node
  /// containing the index entry E
  Nid find_lead(Nid T, Rid r);

  /// Given a leaf node L, from which an entry has been deleted,
  /// eliminate the node if it has too few entries and relocate its
  /// entries. Propagate node elimination upward as necessary. Adjust
  /// all covering rectagles on the path to the root, making them
  /// smaller if possible
  void condense_tree(Node T);


  // TODO signatures
  
  /// Select two entries to be the first elements of the group
  void pick_seeds();
  /// Select one remaining entry for classification in a group.
  void pick_next();
};


#define PRE template <class DATATYPE, class ELEMTYPE>
#define QUAL rtree<DATATYPE, ELEMTYPE>

PRE QUAL::rtree(int dimensions)
    : m_dims(dimensions)
      {
        m_size = 0;
      }

PRE Nid QUAL::choose_leaf(Rid r) { return Nid{}; }

PRE void QUAL::insert(Vec low, Vec high, const DATATYPE &data) {}

PRE size_t QUAL::size() {
  return m_size;
}


}
