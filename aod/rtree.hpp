#pragma once

#include <cstdint>
#include <functional>
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

  using Predicate = std::function<bool(const DATATYPE &)>;
  using SearchCb = std::function<bool(const DATATYPE &)>;

 private:
   /// max entries per node
  const int M = 8; // max entries per node

  /// min entries per node. Could be hardcoded 2 as well, the paper
  /// shows good performance with hardcoded 2 as well
  const int m = M / 2;

  int m_dims;
  
  size_t m_size;
  id_t m_rects_count = 0;
  id_t m_entries_count = 0;
  id_t m_nodes_count = 0;
  id_t m_data_count = 0;
  Nid m_root_id;
  Rid m_temp_rect;

  std::vector<ELEMTYPE> m_rects_low;
  std::vector<ELEMTYPE> m_rects_high;
  std::vector<Node> m_nodes;
  std::vector<Eid> m_node_entries;
  std::vector<Entry> m_entries;
  std::vector<DATATYPE> m_data;

  Rid make_rect_id() {
    Rid res{m_rects_count++};
    m_rects_low.resize(m_rects_count * m_dims, 0);
    m_rects_high.resize(m_rects_count * m_dims, 0);

    return res;
  }

  Nid make_node_id() {
    Nid nid{m_nodes_count++};
    m_nodes.resize(m_nodes_count);
    m_node_entries.resize(m_nodes_count * M);
    for (int i = 0; i < M; i++) {
      // TODO needed?
      set_node_entry(nid, i, make_entry_id());
    }
    return nid;
  }

  Eid make_entry_id() {
    Eid e{m_entries_count++};
    m_entries.resize(m_entries_count);
    Entry &entry = get_entry(e);
    entry.rect_id = make_rect_id();
    return e;
  }

  Did make_data_id() {
    Did d{m_data_count++};
    m_data.resize(m_data_count);
    return d;
  }

  void set_data(Did did, const DATATYPE& data) {
    ASSERT(did);
    m_data[did.id] = data;
  }

  const DATATYPE& get_data(Did did) {
    ASSERT(did);
    return m_data[did.id];
  };

  Eid get_node_entry(Nid n, int idx) { return m_node_entries[n.id * M + idx]; }
  void set_node_entry(Nid n, int idx, Eid e) { m_node_entries[n.id * M + idx] = e; }

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
  void insert(const Vec &low, const Vec &high, const DATATYPE &data);
  size_t size();
  std::vector<DATATYPE> search(const Vec &low, const Vec &high);

 private:
  bool search(Nid, Rid, int &found_count, SearchCb);
  bool rect_contains(Rid bigger, Rid smaller);
  bool rects_overlap(Rid, Rid);

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
        m_root_id = make_node_id();
        m_temp_rect = make_rect_id();
      }

PRE Nid QUAL::choose_leaf(Rid r) { return choose_node(m_root_id, r, 0); }

PRE Nid QUAL::choose_node(Nid n, Rid r, int height) {
  ASSERT(n);
  ASSERT(r);
  Node& node = get_node(n);
  if(node.height == height) {
    return n;
  }
  Nid subtree = choose_subtree(n, r);
  return choose_node(n, r, height);
}

PRE Nid QUAL::choose_subtree(Nid n, Rid r) {
  // CL3. If N is not a leaf, let F be the entry in N whose rectangle
  // FI needs least enlargment to include EI. Resolved ties by
  // choosing the entry with the rectangle of smallest area.

  bool first_run = true;
  ELEMTYPE increase;
  ELEMTYPE best_incr = (ELEMTYPE)-1;
  ELEMTYPE area;
  ELEMTYPE best_area;
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
    if ((increase < best_incr) || first_run) {
      // clearly better (or init)
      best = e;
      best_area = area;
      best_incr = increase;
      first_run = false;
    } else if ((increase == best_incr) && (area < best_area)) {
      // Resolved ties by choosing the entry with the rectangle of
      // smallest area
      best = e;
      best_area = area;
      best_incr = increase;
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

PRE void QUAL::insert(const Vec &low, const Vec &high, const DATATYPE &data) {
  Eid e = make_entry_id(); // also sets the rect
  Entry& entry = get_entry(e);
  entry.data_id = make_data_id();
  set_data(entry.data_id, data);
  Rid r = entry.rect_id;

  for(int i=0; i<m_dims; i++) {
    rect_low_ref(r, i) = low[i];
    rect_high_ref(r, i) = high[i];
  }
  
  Nid n = choose_leaf(entry.rect_id);
  Node& node = get_node(n);
  if(node.count < M) {
    set_node_entry(n, node.count, e);
    ++node.count;
  } else {
    // TODO split
    ASSERT(0);
  }
  ++m_size;
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

PRE std::vector<DATATYPE> QUAL::search(const Vec &low,const Vec &high) {
  ASSERT(low.size() == m_dims);
  ASSERT(high.size() == m_dims);
  std::vector<DATATYPE> res;
  for(int i=0; i<m_dims; ++i) {
    rect_low_ref(m_temp_rect, i) = low[i];
    rect_high_ref(m_temp_rect, i) = high[i];
  }
  SearchCb cb = [&res](const DATATYPE &data) {
    res.push_back(data);
    return true;
  };

  int found_count = 0;
  search(m_root_id, m_temp_rect, found_count, cb);

  return res;
}

PRE bool QUAL::search(Nid n, Rid r, int &found_count, SearchCb cb) {
  Node& node = get_node(n);
  if(node.is_internal()) {
    for(int i=0; i<node.count; ++i) {
      Eid e = get_node_entry(n, i);
      Entry& entry = get_entry(e);
      if(rect_contains(r, entry.rect_id)) {
        if(!search(entry.child_id, r, found_count, cb)) {
          // stop searching
          return false;
        }
      }
    }
  } else {
    // leaf
    for (int i = 0; i < node.count; ++i) {
      Eid e = get_node_entry(n, i);
      Entry &entry = get_entry(e);
      // TODO define some algorithm for search: overlap vs contain, etc..
      if (rects_overlap(r, entry.rect_id)) {
        if (!cb(get_data(entry.data_id))) {
          // stop searching
          return false;
        } else {
          ++found_count;
        }
      }
    }
  }
  return true; // continue searching
}

PRE bool QUAL::rect_contains(Rid bigger, Rid smaller) {
  for (unsigned int index = 0; index < m_dims; ++index) {
    if (rect_low_ref(bigger, index) > rect_low_ref(smaller, index) ||
        rect_high_ref(bigger, index) < rect_high_ref(smaller, index)) {
      return false;
    }
  }
  return true;
}

PRE bool QUAL::rects_overlap(Rid a, Rid b) {
  for (unsigned int index = 0; index < m_dims; ++index) {
    if (rect_low_ref(a, index) > rect_high_ref(b, index) ||
        rect_low_ref(b, index) > rect_high_ref(a, index)) {
      return false;
    }
  }
  return true;
}

} // aod namespace
