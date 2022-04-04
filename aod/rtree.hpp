#pragma once

#include <cstdint>
#include <limits>
#include <ostream>
#include <vector>
#include <assert.h>
#define ASSERT assert
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

/// Node Id
struct Nid : Id {};

/// Entry Id. Each nodes has `m <= count <= M` children
struct Eid : Id {};

/// Data id
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

  int count = 0; ///< children entries count. Between m and M
  int height = 0; ///< zero for leaf, positive otherwise
  Rid rect_id; // TODO rect (mbr) for Node? not present in RTree code
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

 private:
   /// max entries per node
  const int M = 8; // max entries per node

  /// min entries per node. Could be hardcoded 2 as well, the paper
  /// shows good performance with hardcoded 2 as well
  const int m = M / 2;

  int m_dims;
  
  size_t m_size;
  id_t m_rects_count = 0;
  Nid m_root_id;
  Rid m_temp_rect;

  std::vector<ELEMTYPE> m_rects_low;
  std::vector<ELEMTYPE> m_rects_high;
  std::vector<Node> m_nodes;
  std::vector<Eid> m_node_entries;
  std::vector<Entry> m_entries;

  Rid make_rect_id() {
    Rid res{m_rects_count++};
    m_rects_low.resize(m_rects_count * m_dims, 0);
    m_rects_high.resize(m_rects_count * m_dims, 0);

    return res;
  }

  Eid make_entry_id() { return Eid{}; }

  Eid get_node_entry(Nid n, int idx) { return m_node_entries[n.id * M + idx]; }

  Node& get_node(Nid n) {
    return m_nodes[n.id];
  }

  Entry& get_entry(Eid e) {
    return m_entries[e.id];
  }

  ELEMTYPE rect_volume(Rid);

  inline ELEMTYPE &rect_low_ref(Rid r, int dim) {
    return m_rects_low[r.id * m_dims + dim];
  }
  inline ELEMTYPE &rect_high_ref(Rid r, int dim) {
    return m_rects_high[r.id * m_dims + dim];
  }

 public:
  
  rtree() = delete;
  rtree(int dimensions);
  void insert(Vec low, Vec high, const DATATYPE &data);
  size_t size();

 private:


  // Used by insert
  Nid choose_leaf(Rid r);
  Nid choose_node(Nid n, Rid r, int height);
  Nid choose_subtree(Nid n, Rid r);

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

  void combine_rects(Rid a, Rid b, Rid dst);
};


#define PRE template <class DATATYPE, class ELEMTYPE>
#define QUAL rtree<DATATYPE, ELEMTYPE>

PRE QUAL::rtree(int dimensions)
    : m_dims(dimensions)
      {
        m_size = 0;
      }

PRE Nid QUAL::choose_leaf(Rid r) { return choose_node(m_root_id, r, 0); }

PRE Nid QUAL::choose_node(Nid n, Rid r, int height) {
  Node& node = get_node(n);
  if(node.height == height) {
    return n;
  }
  // CL3. If N is not a leaf, let F be the entry in N whose rectangle
  // FI needs least enlargment to include EI. Resolved ties by
  // choosing the entry with the rectangle of smallest area.
}

PRE Nid QUAL::choose_subtree(Nid n, Rid r) {
  bool firstTime = true;
  ELEMTYPE increase;
  ELEMTYPE bestIncr = (ELEMTYPE)-1;
  ELEMTYPE area;
  ELEMTYPE bestArea;
  Eid best;

  Node& node = get_node(n);
  ASSERT(node.count);
  for (int i = 0; i < node.count; i++) {
    Eid e = get_node_entry(n, i);
    Entry& entry = get_entry(e);
    Rid entry_rect = entry.rect_id;
    area = rect_volume(entry_rect);
    combine_rects(r, entry_rect, m_temp_rect);
    increase = rect_volume(m_temp_rect) - area;
    if ((increase < bestIncr) || firstTime) {
      best = e;
      bestArea = area;
      bestIncr = increase;
      firstTime = false;
    } else if ((increase == bestIncr) && (area < bestArea)) {
      best = e;
      bestArea = area;
      bestIncr = increase;
    }
  }
  Entry& entry = get_entry(best);
  return entry.child_id;
}
PRE void QUAL::combine_rects(Rid a, Rid b, Rid dst) {
  for (int index = 0; index < m_dims; index++) {
    rect_low_ref(dst, index) = Min(rect_low_ref(a, index), rect_low_ref(b, index));
    rect_high_ref(dst, index) = Max(rect_high_ref(a, index), rect_high_ref(b, index));
  }
}

PRE void QUAL::insert(Vec low, Vec high, const DATATYPE &data) {
  Eid e = make_entry_id();
}

PRE ELEMTYPE QUAL::rect_volume(Rid r) {
  // RectVolume (TODO add also RectSphericalVolume ?)
  ELEMTYPE volume = (ELEMTYPE)1;

  for (int i = 0; i < m_dims; ++i) {
    volume *= rect_high_ref(r, i) - rect_low_ref(r, i);
  }

  ASSERT(volume >= (ELEMTYPE)0);

  return volume;
}

PRE size_t QUAL::size() { return m_size; }


}
