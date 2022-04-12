#pragma once

#include <assert.h>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <ostream>
#include <set>
#include <vector>
#include <easyprint.hpp>
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
  friend bool operator==(const Id &lhs, const Id &rhs) {
    return lhs.id == rhs.id;
  }
  // for ordered set
  friend bool operator<(const Id &lhs, const Id &rhs) {
    return lhs.id < rhs.id;
  }
};

// Rect Id
struct Rid : Id {};
/// Node Id
struct Nid : Id {};
/// Entry Id. Each node has `m <= count <= M` children
struct Eid : Id {};
/// Data id
struct Did : Id {};

struct Entry {
  Rid rect_id;
  Nid child_id;
  Did data_id;
};

struct Parent {
  Eid entry;
  Nid node;
};

struct Node {
  bool is_internal() const { return (height > 0); }
  bool is_leaf() const { return (height == 0); } // A leaf contains data

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

std::ostream &operator<<(std::ostream &os, const Did &id) {
  os << "Did{" << id.id << "}";
  return os;
}

std::ostream &operator<<(std::ostream &os, const Parent &p) {
  os << "Parent{" << p.entry << "->" << p.node << "}";
  return os;
}

template <class DATATYPE, class ELEMTYPE = double> class rtree {
  static_assert(std::numeric_limits<ELEMTYPE>::is_iec559,
                "'ELEMTYPE' accepts floating-point types only");
  using Vec = std::vector<ELEMTYPE>;

  using Traversal = std::vector<Parent>;
  using Traversals = std::vector<Traversal>;

  using Predicate = std::function<bool(const DATATYPE &)>;
  using SearchCb = std::function<bool(const DATATYPE &)>;

public:
  struct Xml {
    rtree *tree;
    int spaces = 4;
    Xml() = delete;
    Xml(rtree *tree) : tree{tree} {}
    friend std::ostream &operator<<(std::ostream &os, const Xml &o) {
      os << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\" ?>"
         << endl;
      o.tree->to_string(o.spaces, os);
      return os;
    }
  };

private:
  /// Helper struct used in partitioning M+1 entries into two groups
  struct Partition {
    std::vector<Eid> entries; // up to M+1
    Rid groups_rects[2];      // each group's MBR
    ELEMTYPE groups_areas[2]; // area covered by each group's MBR
    Rid temp_rects[2];        // temporary rect for each group, comparing growth
    std::vector<ELEMTYPE> entries_areas;
    Rid cover_rect;
    ELEMTYPE cover_area = 0;
  };

  /// max entries per node (best 8?)
  const int M = 8;

  /// min entries per node. Could be hardcoded 2 as well, the paper
  /// shows good performance with hardcoded 2 as well
  const int m = M / 2;

  int m_dims;

  size_t m_size = 0;
  id_t m_rects_count = 0;
  id_t m_entries_count = 0;
  id_t m_nodes_count = 0;
  id_t m_data_count = 0;
  Nid m_root_id;
  Rid m_temp_rect;
  Partition m_partition;

  Traversal m_traversal;

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

  inline ELEMTYPE rect_volume(Rid);

  inline ELEMTYPE &rect_low_ref(const Rid &r, int dim) {
    return m_rects_low[r.id * m_dims + dim];
  }
  inline ELEMTYPE &rect_high_ref(const Rid &r, int dim) {
    return m_rects_high[r.id * m_dims + dim];
  }
  std::string rect_to_string(Rid);
  void rect_to_string(Rid, std::ostream &os);
  std::string traversal_to_string(const Traversal &traversal);
  std::string node_to_string(Nid nid, int level = 0, int spaces = 4);
  void node_to_string(Nid nid, int level, int spaces, std::ostream &);
  std::string entry_to_string(Eid, int level = 0, int spaces = 4);
  void entry_to_string(Eid, int level, int spaces, std::ostream &);
  bool has_duplicate_nodes();

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
  void search(const Vec &low, const Vec &high, std::vector<DATATYPE> &results);
  int remove(const Vec &low, const Vec &high);
  std::string to_string();
  void to_string(int spaces, std::ostream &);
  Xml to_xml() { return Xml(this); }

private:
  void reinsert_entry(Eid e);
  /// Splits node, places the new M+1 entries (old & new one "e"), & returns the
  /// new node id
  Nid split_and_insert(Nid n, Eid e);

  /// Ascend from leaf node L to the root, adjusting rectangles and
  /// propagating node splits as necessary
  void adjust_tree(const Traversal &traversal, Eid e, Nid nn);
  void adjust_rects(const Traversal &traversal);
  std::array<int, 2> pick_seeds(const std::vector<Eid> &entries);
  void distribute_entries(Nid n, Nid nn, std::vector<Eid> entries,
                          const std::array<int, 2> &seeds);
  void distribute_entries_naive(Nid n, Nid nn, std::vector<Eid> entries);

  bool search(Nid, Rid, int &found_count, const SearchCb&);
  bool rect_contains(Rid bigger, Rid smaller);
  inline bool rects_overlap(Rid, Rid);

  /// insert entry into leaf node: if a split occured, returns a valid new node
  Nid insert(Nid, Eid);
  void plain_insert(Nid, Eid);
  Nid choose_leaf(Rid r, Traversal &traversal);
  Nid choose_node(Nid n, Rid r, int height, Traversal &traversal);
  Eid choose_subtree(Nid n, Rid r);

  /// Given an R-tree whose root node is T, find the leaf node
  /// containing the index entry E
  Nid find_lead(Nid T, Rid r);

  /// Given a leaf node L, from which an entry has been deleted,
  /// eliminate the node if it has too few entries and relocate its
  /// entries. Propagate node elimination upward as necessary. Adjust
  /// all covering rectagles on the path to the root, making them
  /// smaller if possible

  void remove(Nid n, Rid r, int &removed, Traversals &,
              const Traversal &cur_traversal, Predicate cb);

  bool remove_node_entry(Nid n, Eid e);
  void remove_node_entry(Nid n, int idx);

  /// Called after removal of Node(s), by passing the traversal(s) of
  /// the tree that took place for the removed entries
  void condense_tree(const Traversals &);

  int count(Nid n);
  void count(Nid n, int &); // recursive

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
  m_partition.cover_rect = make_rect_id();

  m_partition.entries_areas.resize(M + 1);
}

PRE Nid QUAL::choose_leaf(Rid r, Traversal &traversal) {
  // Eid subtree = choose_subtree(n, r);
  // Nid child = get_entry(subtree).child_id;
  Parent p;
  p.node = m_root_id;
  traversal.push_back(p);
  return choose_node(m_root_id, r, 0, traversal);
}

PRE Nid QUAL::choose_node(Nid n, Rid r, int height, Traversal &traversal) {
  ASSERT(n);
  ASSERT(r);
  Node &node = get_node(n);
  if (node.height == height) {
    return n;
  }
  Eid subtree = choose_subtree(n, r);
  Nid child = get_entry(subtree).child_id;
  Parent p;
  p.entry = subtree;
  p.node = child;
  traversal.push_back(p);

  ASSERT(subtree);
  return choose_node(child, r, height, traversal);
}

PRE Eid QUAL::choose_subtree(Nid n, Rid r) {
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
  return best;
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
  Eid child0 = get_node_entry(n, 0);
  const Entry &child_entry = get_entry(child0);
  for (int i = 0; i < m_dims; i++) {
    rect_low_ref(r, i) = rect_low_ref(child_entry.rect_id, i);
    rect_high_ref(r, i) = rect_high_ref(child_entry.rect_id, i);
  }
  for (int i = 1; i < node.count; ++i) {
    Eid child = get_node_entry(n, i);
    combine_rects(get_entry(child).rect_id, r, r);
  }
}

PRE void QUAL::insert(const Vec &low, const Vec &high, const DATATYPE &data) {
  ++m_size;
  const Eid e = make_entry_id(); // also sets the rect
  Entry &entry = get_entry(e);
  entry.data_id = make_data_id();
  set_data(entry.data_id, data);

  const Rid r = entry.rect_id;

  for (int i = 0; i < m_dims; ++i) {
    rect_low_ref(r, i) = low[i];
    rect_high_ref(r, i) = high[i];
  }

  m_traversal.clear();
  const Nid n = choose_leaf(entry.rect_id, m_traversal);
  const Nid nn = insert(n, e);

  adjust_tree(m_traversal, e, nn);
}

PRE void QUAL::reinsert_entry(Eid e) {
  const Entry &entry = get_entry(e);
  if(entry.child_id) {
    const Node& entry_node = get_node(entry.child_id);
    const Node &root = get_node(m_root_id);
    const bool should_split_reinsert = root.height <= entry_node.height;

    if (should_split_reinsert) {
      const int children_count = get_node(entry.child_id).count;
      for (int i = 0; i < children_count; ++i) {
        reinsert_entry(get_node_entry(entry.child_id, i));
      }
      return;
    }
  }

  m_traversal.clear();
  // TODO traversal needs also root, I'm manually doing it atm
  Parent p_root;
  p_root.node = m_root_id;
  m_traversal.push_back(p_root);
  const int height = entry.child_id ? get_node(entry.child_id).height + 1 : 0;
  const Nid n = choose_node(m_root_id, entry.rect_id, height, m_traversal);
  const Nid nn = insert(n, e);

  adjust_tree(m_traversal, e, nn);
}

PRE Nid QUAL::insert(Nid n, Eid e) {
  Nid nn; // new node (falsy)
  Node &node = get_node(n);
  if (node.count < M) {
    plain_insert(n, e);
  } else {
    nn = split_and_insert(n, e);
    ASSERT(nn);
  }
  // NB: nn (new node - after split_and_insert) isn't yet a part of
  // the tree!
  return nn;
}

PRE void QUAL::plain_insert(Nid n, Eid e) {
  Node &node = get_node(n);
  ASSERT(node.count < M);
  set_node_entry(n, node.count, e);
  ++node.count;
}

PRE Nid QUAL::split_and_insert(Nid n, Eid e) {
  Nid nn = make_node_id();
  Node &node = get_node(n);
  Node &new_node = get_node(nn);
  new_node.height = node.height; // important!
  ASSERT(node.count == M); // if node is not full, makes no sense being here

  std::vector<Eid> &entries = m_partition.entries;
  entries.clear();
  entries.resize(M + 1);

  for (int i = 0; i < node.count; ++i) {
    entries[i] = get_node_entry(n, i);
  }
  entries[M] = e; // the new entry

  std::array<int, 2> seeds = pick_seeds(entries);

  node.count = 0; // TODO free up entries (entry ids) to be reused
  new_node.count = 0;

  distribute_entries(n, nn, entries, seeds);
  // debugging: just placing first half entries to n, rest to nn
  // distribute_entries_naive(n, nn, entries);

  // cout << " split and insert " << n << nn << endl;
  // cout << ":: distributed entries" << endl;
  // cout << node_to_string(n) << endl << node_to_string(nn) << endl;
  ASSERT(node.count + new_node.count == M + 1);
  ASSERT(node.count >= m);
  ASSERT(new_node.count >= m);

  return nn;
}

PRE void QUAL::distribute_entries(Nid n, Nid nn, std::vector<Eid> entries,
                                  const std::array<int, 2> &seeds) {
  const Node &node = get_node(n);
  const Node &new_node = get_node(nn);

  plain_insert(n, entries[seeds[0]]);
  plain_insert(nn, entries[seeds[1]]);

  // initializing helper rects
  copy_rect(get_entry(entries[seeds[0]]).rect_id, m_partition.groups_rects[0]);
  copy_rect(get_entry(entries[seeds[1]]).rect_id, m_partition.groups_rects[1]);
  m_partition.groups_areas[0] = rect_volume(m_partition.groups_rects[0]);
  m_partition.groups_areas[1] = rect_volume(m_partition.groups_rects[1]);

  // quick remove seeds from the entries. First removing the 2nd seed, as it
  // will be bigger, cause it might be the last one
  entries[seeds[1]] = entries.back();
  entries.pop_back();
  entries[seeds[0]] = entries.back();
  entries.pop_back();

  ASSERT(entries.size() == M - 1); // (M+1) -2

  ELEMTYPE biggestDiff;
  int group, entry_index = 0, betterGroup = 0;
  // existing node: group[0], new node: group[1]
  const Nid groups[2] = {n, nn};
  const Node *nodes[2] = {&node, &new_node};
  const int max_fill = M + 1 - m;
  while (!entries.empty() && node.count < max_fill &&
         new_node.count < max_fill) {
    biggestDiff = (ELEMTYPE)-1;
    for (int i = 0; i < entries.size(); ++i) {
      const Entry &entry = get_entry(entries[i]);
      const Rid r = entry.rect_id;
      combine_rects(r, m_partition.groups_rects[0], m_partition.temp_rects[0]);
      combine_rects(r, m_partition.groups_rects[1], m_partition.temp_rects[1]);
      const ELEMTYPE growth0 =
          rect_volume(m_partition.temp_rects[0]) - m_partition.groups_areas[0];
      const ELEMTYPE growth1 =
          rect_volume(m_partition.temp_rects[1]) - m_partition.groups_areas[1];
      ELEMTYPE diff = growth1 - growth0;
      if (diff > 0) {
        group = 0;
      } else {
        group = 1;
        diff = -diff;
      }
      if (diff > biggestDiff) {
        biggestDiff = diff;
        entry_index = i;
        betterGroup = group;
      } else if ((diff == biggestDiff) &&
                 (nodes[group]->count < nodes[betterGroup]->count)) {
        entry_index = i;
        betterGroup = group;
      }
    }
    combine_rects(get_entry(entries[entry_index]).rect_id,
                  m_partition.groups_rects[betterGroup],
                  m_partition.groups_rects[betterGroup]);
    m_partition.groups_areas[betterGroup] =
        rect_volume(m_partition.groups_rects[betterGroup]);

    plain_insert(groups[betterGroup], entries[entry_index]);

    // quick remove:
    entries[entry_index] = entries.back();
    entries.pop_back();
  }
  // exited early cause one node filled too much
  if (node.count + new_node.count < M + 1) {
    Nid to_fill;
    if (node.count >= max_fill) {
      to_fill = nn;
    } else {
      to_fill = n;
    }
    for (const Eid &e : entries) {
      plain_insert(to_fill, e);
    }
  }
}

PRE void QUAL::distribute_entries_naive(Nid n, Nid nn,
                                        std::vector<Eid> entries) {
  int half = entries.size() / 2 + 1;
  for (int i = 0; i < half; ++i) {
    plain_insert(n, entries[i]);
  }
  for (int i = half; i < entries.size(); ++i) {
    plain_insert(nn, entries[i]);
  }
}

PRE void QUAL::adjust_rects(const Traversal &traversal) {
  for (int i = traversal.size() - 1; i >= 0; --i) {
    const Parent &parent = traversal[i];
    if (parent.entry) {
      update_entry_rect(parent.entry);
    }
  }
}

PRE void QUAL::adjust_tree(const Traversal &traversal, Eid e, Nid nn) {
  const bool dbg = false; // e.id == 159;
  if (dbg) {
    // cout << ">>> adjust tree " << pp(traversal) << " " << e << " nn " << nn
         // << endl;
  }
  Nid n;
  for (int level = traversal.size() - 1; level >= 0; --level) {
    const Parent &parent = traversal[level];
    n = parent.node;
    const Node &parent_node = get_node(parent.node);
    // a parent is: entry -> node
    if (parent.entry) {
      // aka not root node. root node doesn't belong to an entry
      if (nn) {
        // if split (new node), need to readjust the parent entry
        // rect. Some leaf nodes now probably not part of the traversal!
        // cout << "    updating parent rect: " << parent << endl;
        update_entry_rect(parent.entry);
      } else {
        // making room for newly inserted entry's rect
        // cout << "    absorbing newly inserted entry into parent" << parent <<
        // e << endl;
        const Entry &parent_entry = get_entry(parent.entry);
        combine_rects(parent_entry.rect_id, get_entry(e).rect_id,
                      parent_entry.rect_id);
      }
    }

    // propagating new node (split) upwards.
    if (nn && level < traversal.size() - 1) {
      const Node &new_node = get_node(nn);
      ASSERT(parent_node.height == new_node.height + 1);

      const Eid ne = make_entry_id();
      Entry &new_entry = get_entry(ne);
      // cout << ".. inserting nn, made entry " << ne << " child nn " << nn <<
      // endl;
      new_entry.child_id = nn;
      update_entry_rect(ne);

      nn = insert(parent.node, ne);
    }
  }
  ASSERT(n == m_root_id);
  if (nn) {
    // Root split!

    // Important: make_node_id before getting nodes! Otherwise,
    // references might be invalid.
    const Nid new_root_id = make_node_id();

    Node &new_root = get_node(new_root_id); // non const - setting height
    const Node &existing_root = get_node(m_root_id);
    const Node &new_node = get_node(nn);

    new_root.height = existing_root.height + 1;
    ASSERT(existing_root.height == new_node.height);

    const Eid old_root_e = make_entry_id();
    const Eid new_node_e = make_entry_id();

    Entry &old_root_entry = get_entry(old_root_e);
    old_root_entry.child_id = m_root_id;
    update_entry_rect(old_root_e);

    Entry &new_node_entry = get_entry(new_node_e);
    new_node_entry.child_id = nn;
    update_entry_rect(new_node_e);

    plain_insert(new_root_id, old_root_e);
    plain_insert(new_root_id, new_node_e);

    m_root_id = new_root_id;
  }
  // cout << "<<< adjust tree" << endl;
}

/// Pick first entry for each group. We get 2 groups after splitting a
/// node. Returns the entries indexes of those seed entries.  Code
/// from Superliminal rtree
PRE std::array<int, 2> QUAL::pick_seeds(const std::vector<Eid> &entries) {
  std::array<int, 2> res;
  ASSERT(entries.size() == M + 1);

  copy_rect(get_entry(entries[0]).rect_id, m_partition.cover_rect);
  for (int i = 0; i < entries.size(); ++i) {
    combine_rects(m_partition.cover_rect, get_entry(entries[i]).rect_id,
                  m_partition.cover_rect);
  }
  m_partition.cover_area = rect_volume(m_partition.cover_rect);

  int seed0 = 0, seed1 = 0;
  ELEMTYPE worst, waste;
  worst = -m_partition.cover_area - 1;
  for (int i = 0; i < entries.size(); ++i) {
    m_partition.entries_areas[i] = rect_volume(get_entry(entries[i]).rect_id);
  }

  for (int i = 0; i < entries.size() - 1; ++i) {
    for (int j = i + 1; j < entries.size(); ++j) {
      combine_rects(get_entry(entries[i]).rect_id,
                    get_entry(entries[j]).rect_id, m_temp_rect);
      waste = rect_volume(m_temp_rect) - m_partition.entries_areas[i] -
              m_partition.entries_areas[j];
      if (waste > worst) {
        worst = waste;
        seed0 = i;
        seed1 = j;
      }
    }
  }
  ASSERT(seed0 != seed1);
  res[0] = seed0;
  res[1] = seed1;

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
  std::vector<DATATYPE> res;

  search(low, high, res);
  return res;
}

PRE void QUAL::search(const Vec &low, const Vec &high,
                      std::vector<DATATYPE> &results) {
  ASSERT(low.size() == m_dims);
  ASSERT(high.size() == m_dims);
  for (int i = 0; i < m_dims; ++i) {
    rect_low_ref(m_temp_rect, i) = low[i];
    rect_high_ref(m_temp_rect, i) = high[i];
  }
  SearchCb cb = [&results](const DATATYPE &data) {
    results.push_back(data);
    return true;
  };

  int found_count = 0;
  search(m_root_id, m_temp_rect, found_count, cb);
}

PRE bool QUAL::search(Nid n, Rid r, int &found_count, const SearchCb &cb) {
  const Node &node = get_node(n);
  if (node.is_internal()) {
    for (int i = 0; i < node.count; ++i) {
      const Eid e = get_node_entry(n, i);
      const Entry &entry = get_entry(e);
      // hm, though I could use rect_contains, but I may be asking for a big
      // rect (search), not a specific rect I have inserted thus, the entry rect
      // (in this case) can not contain the query rect
      // TODO different modes for search? When asking for only one rect (or eg
      // querying a point rect), then I should use rect_contains
      if (rects_overlap(r, entry.rect_id)) {
        if (!search(entry.child_id, r, found_count, cb)) {
          // stop searching
          return false;
        }
      }
    }
  } else {
    // leaf
    for (int i = 0; i < node.count; ++i) {
      const Eid e = get_node_entry(n, i);
      const Entry &entry = get_entry(e);
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

PRE int QUAL::remove(const Vec &low, const Vec &high) {
  ASSERT(low.size() == m_dims);
  ASSERT(high.size() == m_dims);
  for (int i = 0; i < m_dims; ++i) {
    rect_low_ref(m_temp_rect, i) = low[i];
    rect_high_ref(m_temp_rect, i) = high[i];
  }
  int removed = 0;
  Traversals remove_traverals;
  Predicate cb = [](const DATATYPE &data) { return true; };
  Traversal cur_traversal;
  Parent p;
  p.node = m_root_id;
  cur_traversal.push_back(p);
  remove(m_root_id, m_temp_rect, removed, remove_traverals, cur_traversal, cb);

  // sorting remove_traverals: longer traverals first, thus leaf node entries (data) removed first
  std::sort(remove_traverals.begin(), remove_traverals.end(),
            [](const Traversal &a, const Traversal &b) {
              return a.size() > b.size();
            });
  condense_tree(remove_traverals);
  m_size -= removed;
  return removed;
}

PRE void QUAL::remove(Nid n, Rid r, int &counter, Traversals &traversals,
                      const Traversal &cur_traversal, Predicate cb) {
  // cout << "cur traversal " << pp(cur_traversal) << endl;
  Node &node = get_node(n);
  if (node.is_internal()) {
    for (int i = 0; i < node.count; ++i) {
      Eid e = get_node_entry(n, i);
      const Entry &entry = get_entry(e);
      if (rects_overlap(r, entry.rect_id)) {
        Traversal sub_traversal = cur_traversal;
        Parent p;
        p.entry = e;
        p.node = entry.child_id;
        sub_traversal.push_back(p);
        // cout << "here , cur_traversal " << pp(cur_traversal) << " sub " <<
        // pp(sub_traversal) << endl;
        remove(entry.child_id, r, counter, traversals, sub_traversal, cb);
      }
    }
  } else {
    // leaf
    bool removed = false;
    for (int i = 0; i < node.count; ++i) {
      Eid e = get_node_entry(n, i);
      const Entry &entry = get_entry(e);
      if (rects_overlap(r, entry.rect_id)) {
        if (cb(get_data(entry.data_id))) {
          removed = true;
          ++counter;
          remove_node_entry(n, i);
          --i; // need to revisit i again (remove_node_entry shuffles entries)
        }
      }
    }
    if (removed) {
      traversals.push_back(cur_traversal);
    }
  }
}

PRE bool QUAL::remove_node_entry(Nid n, Eid e) {
  Node &node = get_node(n);
  int idx = -1;
  for (int i = 0; i < node.count; ++i) {
    if (get_node_entry(n, i) == e) {
      idx = i;
      break;
    }
  }
  if (idx == -1)
    return false;

  set_node_entry(n, idx, get_node_entry(n, node.count - 1));
  --node.count;
  return true;
}

PRE void QUAL::remove_node_entry(Nid n, int idx) {
  Node &node = get_node(n);
  ASSERT(idx < node.count);
  set_node_entry(n, idx, get_node_entry(n, node.count - 1));
  --node.count;
}

PRE void QUAL::condense_tree(const Traversals &traversals) {
  // cout << ">>> condense tree" << endl;
  std::set<Eid> entries_to_reinsert;
  for (const Traversal &traversal : traversals) {
    // cout << "    " << pp(traversal) << endl;
    for (int i = traversal.size() - 1; i >= 1; --i) {
      const Parent &x = traversal[i];
      Node &node = get_node(x.node);
      if (node.count < m) {
        const Parent &p = traversal[i - 1];
        if (node.count == 0) {
          // should I just go on? seems like this was already processed
        }
        // Remove X from its parent p
        for (int j = 0; j < node.count; ++j) {
          // Add all children of X to S
          entries_to_reinsert.insert(get_node_entry(x.node, j));
        }
        node.count = 0;
        // this might return false, but it's fine
        // cout << " removing node entry " << p.node << x.entry << endl;
        remove_node_entry(p.node, x.entry);
      }
    }
    adjust_rects(traversal);
  }

  // if root has no children, make it a leaf node?
  // what if a node with height 1 was removed & then will be reinserted?
  Node &root = get_node(m_root_id);
  if (root.count == 0) {
    // cout << "root empty!, setting height to 0 " << endl;
    root.height = 0;
  }
  // skipping reinserting entries that point to empty node
  std::set<Eid>::iterator it = entries_to_reinsert.begin();
  while (it != entries_to_reinsert.end()) {
    auto current = it++;
    const Entry &entry = get_entry(*current);
    const Nid n = entry.child_id;
    if (n && get_node(n).count == 0) {
      // hmm not sure about this
      // cout << "remove from reinsert list?: " << *current << endl;
      entries_to_reinsert.erase(current);
    }
  }

  // cout << "    reinsert entries " << pp(entries_to_reinsert) << endl;
  for (Eid e : entries_to_reinsert) {
    // cout << "    reinserting " << e << endl; // << entry_to_string(e) << endl;
    reinsert_entry(e);
  }
  // cout << "<<< condense tree" << endl;
}

PRE int QUAL::count(Nid n) {
  int counter = 0;
  count(n, counter);
  return counter;
}

PRE void QUAL::count(Nid n, int &counter) {
  const Node &node = get_node(n);
  if (node.is_internal()) {
    for (int i = 0; i < node.count; ++i) {
      const Entry &entry = get_entry(get_node_entry(n, i));
      ASSERT(entry.child_id);
      count(entry.child_id, counter);
    }
  } else {
    // leaf
    counter += node.count;
  }
}

PRE bool QUAL::has_duplicate_nodes() {
  std::set<id_t> nodes;

  std::function<bool(Nid)> traverse = [&](Nid n) {
    Node &node = get_node(n);
    if (nodes.count(n.id) > 0) {
      cout << "Node " << n << " already in the container" << endl;
      return true;
    }
    nodes.insert(n.id);
    if (!node.is_internal()) {
      // leaf
      return false;
    }
    for (int i = 0; i < node.count; ++i) {
      Eid e = get_node_entry(n, i);
      Entry &entry = get_entry(e);
      if (traverse(entry.child_id)) {
        return true;
      }
    }
    return false;
  };

  return traverse(m_root_id);
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

PRE void QUAL::entry_to_string(Eid e, int level, int spaces, std::ostream &os) {
  if (!e)
    return;
  auto indent = [&]() {
    for (int i = 0; i < level * spaces; ++i) {
      os << " ";
    }
  };
  Entry entry = get_entry(e);
  Rid r = entry.rect_id;
  indent();
  os << " <Entry id=\"" << e.id << "\"";
  // os << " rect_id=\"" << r.id << "\"";
  os << " mbr=\"";
  rect_to_string(r, os);
  os << "\"";
  if (entry.data_id) {
    // self closing entry, adding data-id
    // os << " data-id=\"" << entry.data_id.id << "\"";
    os << " data=\"" << pp(get_data(entry.data_id)) << "\"";
    os << " />" << endl;
  } else if (entry.child_id) {
    // child node
    os << " >" << endl;
    node_to_string(entry.child_id, level + 1, spaces, os);
    indent();
    os << " </Entry>" << endl;
  } else {
    // should not reach this
    os << " />" << endl;
  }
}

PRE std::string QUAL::entry_to_string(Eid e, int level, int spaces) {
  std::ostringstream os;
  entry_to_string(e, level, spaces, os);
  return os.str();
};

PRE void QUAL::node_to_string(Nid n, int level, int spaces, std::ostream &os) {
  if (!n)
    return;
  Node node = get_node(n);
  auto indent = [&]() {
    for (int i = 0; i < level * spaces; ++i) {
      os << " ";
    }
  };
  indent();
  os << "<Node id=\"" << n.id << "\" height=\"" << node.height
     << "\" children=\"" << node.count << "\" >" << endl;
  for (int i = 0; i < node.count; i++) {
    Eid e = get_node_entry(n, i);
    entry_to_string(e, level, spaces, os);
  }
  indent();
  os << "</Node>" << endl;
};

PRE std::string QUAL::node_to_string(Nid nid, int level, int spaces) {
  std::ostringstream os;
  node_to_string(nid, level, spaces, os);
  return os.str();
};

PRE std::string QUAL::to_string() {
  std::ostringstream os;
  to_string(os);
  return os.str();
}

PRE void QUAL::to_string(int spaces, std::ostream &os) {
  node_to_string(m_root_id, 0, spaces, os);
}

PRE void QUAL::rect_to_string(Rid rid, std::ostream &os) {
  ASSERT(rid);
  for (int i = 0; i < m_dims; i++) {
    if (i == 0) {
      os << "{";
    }
    os << rect_low_ref(rid, i);
    if (i != m_dims - 1) {
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

}

PRE std::string QUAL::rect_to_string(Rid rid) {
  std::ostringstream os;
  rect_to_string(rid, os);
  return os.str();
}

PRE std::string QUAL::traversal_to_string(const Traversal &traversal) {
  std::ostringstream os;
  os << "Traversal{" << pp(traversal) << "}";

  return os.str();
}

} // namespace aod
#undef pp
#undef PRE
#undef QUAL
