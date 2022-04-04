#pragma once

#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <ostream>
#include <vector>
#include <assert.h>
#define ASSERT assert
#define pp(x) easyprint::stringify(x)
namespace aod {

using std::cout;
using std::endl;

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

  std::vector<Nid> m_traversal;

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
  /// Splits node, places the new M+1 entries (old & new one "e"), & returns the new node id
  Nid split_and_insert(Nid n, Eid e);
  void adjust_tree(Nid n, Nid nn, Eid e, const std::vector<Nid>& traversal);
  std::array<Eid, 2> pick_seeds(Nid n, Rid r);
  bool search(Nid, Rid, int &found_count, SearchCb);
  bool rect_contains(Rid bigger, Rid smaller);
  bool rects_overlap(Rid, Rid);

  /// insert entry into leaf node: if a split occured, returns a valid new node
  Nid insert(Nid, Eid, const std::vector<Nid> &traversal);
  void plain_insert(Nid, Eid);
  Nid choose_leaf(Rid r, std::vector<Nid> &traversal);
  Nid choose_node(Nid n, Rid r, int height, std::vector<Nid> &traversal);
  Nid choose_subtree(Nid n, Rid r);

  /// Ascend from leaf node L to the root, adjusting rectangles and
  /// propagating node splits as necessary
  void adjust_tree(Nid n, Nid nn);

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

PRE Nid QUAL::choose_leaf(Rid r, std::vector<Nid> &traversal) { return choose_node(m_root_id, r, 0, traversal); }

PRE Nid QUAL::choose_node(Nid n, Rid r, int height, std::vector<Nid> &traversal) {
  ASSERT(n);
  ASSERT(r);
  Node& node = get_node(n);
  if(node.height == height) {
    return n;
  }
  traversal.push_back(n);
  Nid subtree = choose_subtree(n, r);
  return choose_node(n, r, height, traversal);
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

  m_traversal.clear();
  Nid n = choose_leaf(entry.rect_id, m_traversal);
  insert(n, e, m_traversal);
  ++m_size;
}

PRE Nid QUAL::insert(Nid n, Eid e, const std::vector<Nid>& traversal) {
  Nid nn; // new node (falsy)
  Node &node = get_node(n);
  if (node.count < M) {
    plain_insert(n, e);
  } else {
    nn = split_and_insert(n, e);
    ASSERT(nn);
    cout << "split : " << nn << "traversal " << pp(traversal) << endl;
    // TODO
    //
    // We need to divided the M+1 entries between two nodes.  We
    // need to partition the set of M+1 entries into two groups, one
    // for each node.
    // ASSERT(0);
  }
  adjust_tree(n, nn, e, traversal);
}

PRE void QUAL::plain_insert(Nid n, Eid e) {
  Node &node = get_node(n);
  ASSERT(node.count < M);
  set_node_entry(n, node.count, e);
  ++node.count;
}

PRE Nid QUAL::split_and_insert(Nid n, Eid e) {
  // std::vector<Eid> entries;
  Nid nn = make_node_id();
  // seends are the 2 entries that will be placed first into the 2 new nodes

  // QS1: Apply algorithm PickSeeds to choose two entries to be the first elemetns fo the groups
  // std::array<Eid, 2> seeds = pick_seeds(n, r);

  // QS1: Sssign each to a group
  // set_node_entry(groups[0], 0, seeds[0]);
  // set_node_entry(groups[1], 0, seeds[1]);

  // TODO picknext
  // naive atm
  // Node& node = get_node(n);
  // for(int i=0; i< half; ++i) {
  //   set_node_entry(groups[0], i+1, get_node_entry(n, i));
  // }

  // for (int i = half; i < node.count; ++i) {
  //   set_node_entry(groups[1], (i-half)+1, get_node_entry(n, i));
  // }

  // naively move half the entries into the new branch

  Node& node = get_node(n);
  int half = node.count / 2;
  for(int i=half; i< node.count; ++i) {
    plain_insert(nn, get_node_entry(n, i));
  }
  plain_insert(nn, e); // adding also the new entry
  node.count = half;

  return nn;
}

PRE void QUAL::adjust_tree(Nid n, Nid nn, Eid e, const std::vector<Nid>& traversal) {
  for (auto it = traversal.rbegin(); it != traversal.rend(); ++it) {
    Nid parent = *it;
    cout << "adjusting parent " << parent << endl;
    // TODO adjusting MBR
  }
  if(n == m_root_id && nn) {
    // static std::vector<Nid> no_traversal;
    // Nid new_root = insert(m_root, Eid e)
    // root split
    // 1 create entry with rect covering nn (new node)
    Eid new_entry = make_entry_id(); // TODO make_entry_from_node
    // 2
    Nid new_root = split_and_insert(n, new_entry);
    // TODO adjust MBR
    
  }
}

/// Linear Pick Seeds. Pick first entry for each group. We get 2
/// groups after splittinga node.
PRE std::array<Eid, 2> QUAL::pick_seeds(Nid n, Rid r) {
  std::array<Eid, 2> res;


  // TODO
  res[0] = get_node_entry(n, 0);
  res[1] = get_node_entry(n, 1);

  // [Find extreme rectangles along all dimensions] Along each dimension, find the entry whose rectangle has the highest low side, and the one with the lowest high side. Record the separation.

  // [Adjust for shape of the rectangle cluster] Normalize the
  // separations by dividing by the width of the entire set along the
  // corresponding dimension

  /**
     From RTree

  int seed0 = 0, seed1 = 0;
  ELEMTYPEREAL worst, waste;
  ELEMTYPEREAL area[MAXNODES + 1];

  for (int index = 0; index < a_parVars->m_total; ++index) {
    area[index] = CalcRectVolume(&a_parVars->m_branchBuf[index].m_rect);
  }

  worst = -a_parVars->m_coverSplitArea - 1;
  for (int indexA = 0; indexA < a_parVars->m_total - 1; ++indexA) {
    for (int indexB = indexA + 1; indexB < a_parVars->m_total; ++indexB) {
      Rect oneRect = CombineRect(&a_parVars->m_branchBuf[indexA].m_rect,
                                 &a_parVars->m_branchBuf[indexB].m_rect);
      waste = CalcRectVolume(&oneRect) - area[indexA] - area[indexB];
      if (waste > worst) {
        worst = waste;
        seed0 = indexA;
        seed1 = indexB;
      }
    }
  }
   */

  return res;
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
#undef pp
