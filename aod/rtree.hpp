#pragma once

#include <assert.h>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <ostream>
#include <vector>
#define ASSERT assert
#define pp(x) easyprint::stringify(x)
namespace aod {

using std::cout;
using std::endl;

using id_t = uint32_t;
struct Id {
  static const id_t nullid = std::numeric_limits<id_t>::max();
  id_t id = nullid; //
  operator bool() const { return id != nullid; }
  Id() : Id(nullid){};
  Id(const Id &other) : Id(other.id){};
  Id(id_t id) : id{id} {}
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

  int count = 0;  ///< children entries count. Between m and M
  int height = 0; ///< zero for leaf, positive otherwise
  Rid rect_id;    // TODO rect (mbr) for Node? not present in RTree code
};

std::ostream &operator<<(std::ostream &os, const Eid &id) {
  os << "Eid{" << id.id << "}";
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

template <class DATATYPE, class ELEMTYPE = double> class rtree {
  static_assert(std::numeric_limits<ELEMTYPE>::is_iec559,
                "'COORD_TYPE' accepts floating-point types only");

  typedef std::vector<ELEMTYPE> Vec;

  using Predicate = std::function<bool(const DATATYPE &)>;
  using SearchCb = std::function<bool(const DATATYPE &)>;

private:
  /// Helper struct used in partitioning M+1 entries into two groups
  struct Partition {
    std::vector<Eid> entries; // up to M+1
    Rid groups_rects[2];             // each group's MBR
    ELEMTYPE groups_areas[2];        // area covered by each group's MBR
    Rid temp_rects[2]; // temporary rect for each group, comparing growth
    std::vector<ELEMTYPE> entries_areas;
  };

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
  Partition m_partition;

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

  void set_data(Did did, const DATATYPE &data) {
    ASSERT(did);
    m_data[did.id] = data;
  }

  const DATATYPE &get_data(Did did) {
    ASSERT(did);
    return m_data[did.id];
  };

  Eid get_node_entry(Nid n, int idx) { return m_node_entries[n.id * M + idx]; }
  void set_node_entry(Nid n, int idx, Eid e) {
    m_node_entries[n.id * M + idx] = e;
  }

  Node &get_node(Nid n) { return m_nodes[n.id]; }

  Entry &get_entry(Eid e) { return m_entries[e.id]; }

  ELEMTYPE rect_volume(Rid);

  inline ELEMTYPE &rect_low_ref(const Rid &r, int dim) {
    return m_rects_low[r.id * m_dims + dim];
  }
  inline ELEMTYPE &rect_high_ref(const Rid &r, int dim) {
    return m_rects_high[r.id * m_dims + dim];
  }
  std::string rect_to_string(Rid);
  std::string node_to_string(Nid nid, int level = 0);

  void copy_rect(Rid src, Rid dst) {
    for (auto i = 0; i < m_dims; ++i) {
      rect_low_ref(dst, i) = rect_low_ref(src, i);
      rect_high_ref(dst, i) = rect_high_ref(src, i);
    }
  }

public:
  rtree() = delete;
  rtree(int dimensions);
  void insert(const Vec &low, const Vec &high, const DATATYPE &data);
  size_t size();
  std::vector<DATATYPE> search(const Vec &low, const Vec &high);
  std::string to_string();

private:
  /// Splits node, places the new M+1 entries (old & new one "e"), & returns the
  /// new node id
  Nid split_and_insert(Nid n, Eid e);
  void adjust_tree(Nid n, Nid nn, Eid e, const std::vector<Nid> &traversal);
  std::array<Eid, 2> pick_seeds(const std::vector<Eid> &entries);
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

  void combine_rects(Rid a, Rid b, Rid dst);
  void update_entry_rect(Eid e);
};

#define PRE template <class DATATYPE, class ELEMTYPE>
#define QUAL rtree<DATATYPE, ELEMTYPE>

PRE QUAL::rtree(int dimensions) : m_dims(dimensions) {
  m_size = 0;
  m_root_id = make_node_id();
  m_temp_rect = make_rect_id();
  m_partition.entries.resize(M + 1);

  m_partition.groups_rects[0] = make_rect_id();
  m_partition.groups_rects[1] = make_rect_id();

  m_partition.temp_rects[0] = make_rect_id();
  m_partition.temp_rects[1] = make_rect_id();

  m_partition.entries_areas.resize(M+1);
}

PRE Nid QUAL::choose_leaf(Rid r, std::vector<Nid> &traversal) {
  return choose_node(m_root_id, r, 0, traversal);
}

PRE Nid QUAL::choose_node(Nid n, Rid r, int height,
                          std::vector<Nid> &traversal) {
  ASSERT(n);
  ASSERT(r);
  Node &node = get_node(n);
  if (node.height == height) {
    return n;
  }
  traversal.push_back(n);
  Nid subtree = choose_subtree(n, r);
  ASSERT(subtree);
  return choose_node(subtree, r, height, traversal);
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

  Node &node = get_node(n);
  ASSERT(node.count);

  for (int i = 0; i < node.count; i++) {
    Eid e = get_node_entry(n, i);
    Entry &entry = get_entry(e);
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
  Entry &entry = get_entry(best);
  return entry.child_id;
}

PRE void QUAL::combine_rects(Rid a, Rid b, Rid dst) {
  for (int i = 0; i < m_dims; i++) {
    rect_low_ref(dst, i) = Min(rect_low_ref(a, i), rect_low_ref(b, i));
    rect_high_ref(dst, i) = Max(rect_high_ref(a, i), rect_high_ref(b, i));
  }
}

PRE void QUAL::update_entry_rect(Eid e) {
  const Entry &the_entry = get_entry(e);
  Rid r = the_entry.rect_id;
  Nid n = the_entry.child_id;
  if (!n)
    return;
  Node &node = get_node(n);
  if (!node.count)
    return;
  // going through the child node entries
  Eid child_e = get_node_entry(n, 0);
  const Entry &child_entry = get_entry(child_e);
  for (int i = 0; i < m_dims; i++) {
    rect_low_ref(r, i) = rect_low_ref(child_entry.rect_id, i);
    rect_high_ref(r, i) = rect_high_ref(child_entry.rect_id, i);
  }
  for (int i = 1; i < node.count; ++i) {
    Eid e = get_node_entry(n, i);
    combine_rects(get_entry(e).rect_id, r, r);
  }
}

PRE void QUAL::insert(const Vec &low, const Vec &high, const DATATYPE &data) {
  // cout << "insert " << pp(data) << endl;
  Eid e = make_entry_id(); // also sets the rect
  Entry &entry = get_entry(e);
  entry.data_id = make_data_id();
  set_data(entry.data_id, data);
  Rid r = entry.rect_id;

  for (int i = 0; i < m_dims; ++i) {
    rect_low_ref(r, i) = low[i];
    rect_high_ref(r, i) = high[i];
  }

  m_traversal.clear();
  Nid n = choose_leaf(entry.rect_id, m_traversal);
  insert(n, e, m_traversal);
  ++m_size;
}

PRE Nid QUAL::insert(Nid n, Eid e, const std::vector<Nid> &traversal) {
  Nid nn; // new node (falsy)
  Node &node = get_node(n);
  const DATATYPE& data = get_data(get_entry(e).data_id);
  if (node.count < M) {
    cout << "insert plain : " << pp(data) << " " << n << endl;
    plain_insert(n, e);
  } else {
    nn = split_and_insert(n, e);
    ASSERT(nn);
    cout << "insert split : " << pp(data) << " " << n << ", " << nn << " traversal " << pp(traversal) << endl;
    cout << node_to_string(n) << endl << node_to_string(nn) << endl;
    // exit(0);
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
  Node &node = get_node(n);
  Node &new_node = get_node(nn);
  ASSERT(node.count == M); // if node is not full, makes no sense being here

  std::vector<Eid> &split_entries = m_partition.entries;
  split_entries.clear();
  split_entries.resize(M + 1);

  for (int i = 0; i < node.count; ++i) {
    split_entries[i] = get_node_entry(n, i);
  }
  split_entries[M] = e; // the new entry
  const Entry &entry = get_entry(e);
  std::array<Eid, 2> seeds = pick_seeds(split_entries);
  node.count = 0; // TODO free up entries
  new_node.count = 0;

  plain_insert(n, seeds[0]);
  plain_insert(nn, seeds[1]);

  // initializing helper rects
  copy_rect(get_entry(seeds[0]).rect_id, m_partition.groups_rects[0]);
  copy_rect(get_entry(seeds[1]).rect_id, m_partition.groups_rects[1]);
  m_partition.groups_areas[0] = rect_volume(m_partition.groups_rects[0]);
  m_partition.groups_areas[1] = rect_volume(m_partition.groups_rects[1]);

  // TODO implement Eid == operator
  split_entries.erase(
      std::remove_if(split_entries.begin(), split_entries.end(),
                     [&seeds](const Eid &o) { return o.id == seeds[0].id; }));
  split_entries.erase(
      std::remove_if(split_entries.begin(), split_entries.end(),
                     [&seeds](const Eid &o) { return o.id == seeds[1].id; }));

  cout << "--- seeds " << seeds[0] << ", " << seeds[1] << endl << node_to_string(n) << endl << node_to_string(nn) << endl << " split entries " << pp(split_entries) << endl;;

  ASSERT(split_entries.size() == M-1); // M+1 -2

  ELEMTYPE biggestDiff;
  int group, entry_index = 0, betterGroup = 0;
  // existing node: group[0], new node: group[1]
  Nid groups[2] = {n, nn};
  int* groups_count[2] = {
    &node.count,
    &new_node.count
  };
  int max_fill = M+1-m;
  while(!split_entries.empty() &&
        node.count < max_fill &&
        new_node.count < max_fill
        ) {
    for(int i=0; i<split_entries.size(); ++i) {
      const Entry& entry = get_entry(split_entries[i]);
      Rid r = entry.rect_id;
      combine_rects(r, m_partition.groups_rects[0], m_partition.temp_rects[0]);
      combine_rects(r, m_partition.groups_rects[1], m_partition.temp_rects[1]);
      ELEMTYPE growth0 = rect_volume(m_partition.temp_rects[0]) - m_partition.groups_areas[0];
      ELEMTYPE growth1 = rect_volume(m_partition.temp_rects[1]) - m_partition.groups_areas[1];
      ELEMTYPE diff = growth1 - growth0;
      if(diff > 0 ) {
        group = 0;
      } else {
        group = 1;
        diff = -diff;
      }
      if (diff > biggestDiff) {
        biggestDiff = diff;
        entry_index = i;
        betterGroup = group;
      } else if ((diff == biggestDiff) && (*groups_count[group] <
                                           *groups_count[betterGroup])) {
        entry_index = i;
        betterGroup = group;
      }
    }
    // put chose into the betterGroup
    combine_rects(get_entry(split_entries[entry_index]).rect_id, m_partition.groups_rects[betterGroup], m_partition.groups_rects[betterGroup]);
    m_partition.groups_areas[betterGroup] = rect_volume(m_partition.groups_rects[betterGroup]);
    plain_insert(groups[betterGroup], split_entries[entry_index]);
    split_entries.erase(split_entries.begin() + entry_index);
  }

  // exited early cause one node filled too much
  if(node.count + new_node.count < M+1 ) {
    Nid to_fill;
    if(node.count >= max_fill) {
      to_fill = nn;
    } else {
      to_fill = n;
    }
    for(const Eid& e : split_entries) {
      plain_insert(to_fill, e);
    }
  }

  ASSERT(node.count + new_node.count == M+1);
  ASSERT(node.count >=m);
  ASSERT(new_node.count >=m);

  return nn;
}

PRE void QUAL::adjust_tree(Nid n, Nid nn, Eid e,
                           const std::vector<Nid> &traversal) {
  for (auto it = traversal.rbegin(); it != traversal.rend(); ++it) {
    Nid parent = *it;
    cout << "adjusting parent " << parent << endl;
    // TODO adjusting MBR
    // TODO if nn, try to insert to parent
  }
  if (n == m_root_id && nn) {
    cout << "root split" << endl;
    Node &root = get_node(n);
    Node &new_node = get_node(nn);
    ASSERT(root.height == new_node.height);

    Eid old_root_e = make_entry_id();
    Entry &old_root_entry = get_entry(old_root_e);
    old_root_entry.child_id = m_root_id;
    update_entry_rect(old_root_e);

    Eid new_node_e = make_entry_id();
    Entry &new_node_entry = get_entry(new_node_e);
    new_node_entry.child_id = nn;
    update_entry_rect(new_node_e);

    Nid new_root_id = make_node_id();
    Node &new_root = get_node(new_root_id);
    new_root.height = root.height + 1;
    plain_insert(new_root_id, old_root_e);
    plain_insert(new_root_id, new_node_e);

    m_root_id = new_root_id;
    cout << "root split done, new root " << m_root_id << ", children " << n << ", " << nn
         << endl;

    // static std::vector<Nid> no_traversal;
    // Nid new_root = insert(m_root, Eid e)
    // root split
    // 1 create entry with rect covering nn (new node)
    // Eid new_entry = make_entry_id(); // TODO make_entry_from_node
    // 2
    // Nid new_root = split_and_insert(n, new_entry);
    // TODO adjust MBR
  }
}

/// Linear Pick Seeds. Pick first entry for each group. We get 2
/// groups after splittinga node.
PRE std::array<Eid, 2> QUAL::pick_seeds(const std::vector<Eid> &entries) {
  std::array<Eid, 2> res;
  ASSERT(entries.size() == M+1);

  int seed0 = 0, seed1 = 0;
  ELEMTYPE worst, waste;
  worst = std::numeric_limits<ELEMTYPE>::min();
  for(int i=0; i< entries.size(); ++i) {
    m_partition.entries_areas[i] = rect_volume(get_entry(entries[i]).rect_id);
  }

  for (int i = 0; i < entries.size() - 1; ++i) {
    for (int j = i + 1; j < entries.size(); ++j) {
      combine_rects(get_entry(entries[i]).rect_id, get_entry(entries[j]).rect_id, m_temp_rect);
      waste = rect_volume(m_temp_rect) - m_partition.entries_areas[i] - m_partition.entries_areas[j];
      if (waste > worst) {
        worst = waste;
        seed0 = i;
        seed1 = j;
      }
    }
  }
  ASSERT(seed0 != seed1);
  res[0] = entries[seed0];
  res[1] = entries[seed1];

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

PRE std::vector<DATATYPE> QUAL::search(const Vec &low, const Vec &high) {
  ASSERT(low.size() == m_dims);
  ASSERT(high.size() == m_dims);
  std::vector<DATATYPE> res;
  for (int i = 0; i < m_dims; ++i) {
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
  Node &node = get_node(n);
  if (node.is_internal()) {
    for (int i = 0; i < node.count; ++i) {
      Eid e = get_node_entry(n, i);
      Entry &entry = get_entry(e);
      // hm, though I could use rect_contains, but I may be asking for a big rect (search), not a specific rect I have inserted
      // thus, the entry rect (in this case) can not contain the query rect
      // TODO different modes for search? When asking for only one rect (or eg querying a point rect), then I should use rect_contains
      if (rects_overlap(entry.rect_id ,r)) {
        if (!search(entry.child_id, r, found_count, cb)) {
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
PRE std::string QUAL::node_to_string(Nid nid, int level) {
  std::ostringstream os;
  Node node = get_node(nid);
  std::string indent;
  for (int i = 0; i < level; ++i) {
    indent += "  ";
  }
  os << indent << level << ":" << nid << endl;
  for (int i = 0; i < node.count; i++) {
    Eid e = get_node_entry(nid, i);
    Entry entry = get_entry(e);
    Rid r = entry.rect_id;
    os << indent << e << " " << rect_to_string(r);
    if (node.height == 0) {
      os << " data: " << pp(get_data(entry.data_id)) << endl;
    } else {
      os << " child: " << entry.child_id << endl;
      os << node_to_string(entry.child_id, level + 1);
    }
  }
  return os.str();
};

PRE std::string QUAL::to_string() {
  return node_to_string(m_root_id, 0);
}

PRE std::string QUAL::rect_to_string(Rid rid) {
  std::ostringstream os;

  ASSERT(rid);
  for(int i=0; i<m_dims; i++) {
    if(i==0) {
      os << "{";
    }
    os << rect_low_ref(rid, i);
    if(i != m_dims-1) {
      os << ", ";
    } else {
      os << "}";
    }
  }

  os << "...";
  
  for (int i = 0; i < m_dims; i++) {
    if (i == 0) {
      os << "{";
    }
    os << rect_high_ref(rid, i);
    if (i != m_dims - 1) {
      os << ", ";
    } else {
      os << "}";
    }
  }

  return os.str();
}

} // namespace aod
#undef pp
