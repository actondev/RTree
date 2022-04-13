#include "./rtree_base.hpp"
// #include <easyprint.hpp>
// #define pp(x) easyprint::stringify(x)

#include <iostream>
#include <sstream>
#include <set>

#include <assert.h>
#define ASSERT assert

#ifndef Min
#define Min std::min
#endif // Min
#ifndef Max
#define Max std::max
#endif // Max


namespace aod {

using std::cout;
using std::endl;

std::ostream &operator<<(std::ostream &os, const rtree_base::Eid &id) {
  os << "Eid{" << id.id << "}";
  return os;
}

std::ostream &operator<<(std::ostream &os, const rtree_base::Nid &id) {
  os << "Nid{" << id.id << "}";
  return os;
}

std::ostream &operator<<(std::ostream &os, const rtree_base::Rid &id) {
  os << "Rid{" << id.id << "}";
  return os;
}

std::ostream &operator<<(std::ostream &os, const rtree_base::Did &id) {
  os << "Did{" << id.id << "}";
  return os;
}

std::ostream &operator<<(std::ostream &os, const rtree_base::Parent &p) {
  os << "Parent{" << p.entry << "->" << p.node << "}";
  return os;
}

std::ostream &operator<<(std::ostream &os, const rtree_base::Traversal& tr) {
  os << "Traversal{ " << endl;
  for(auto const &parent : tr) {
    os << parent << ", ";
  }
  os << " }";
  return os;
}

rtree_base::Rid rtree_base::make_rect_id() {
  Rid res{m_rects_count++};
  m_rects_low.resize(m_rects_count * m_dims, 0);
  m_rects_high.resize(m_rects_count * m_dims, 0);

  return res;
}
rtree_base::Nid rtree_base::make_node_id() {
  Nid nid{m_nodes_count++};
  m_nodes.resize(m_nodes_count);
  m_node_entries.resize(m_nodes_count * M);
  return nid;
}
rtree_base::Eid rtree_base::make_entry_id() {
  Eid e{m_entries_count++};
  m_entries.resize(m_entries_count);
  Entry &entry = get_entry(e);
  entry.rect_id = make_rect_id();
  return e;
}
rtree_base::Did rtree_base::make_data_id() {
  Did d{m_data_count++};
  // m_data.resize(m_data_count);
  return d;
}

inline rtree_base::Node &rtree_base::get_node(Nid n) { return m_nodes[n.id]; }
inline rtree_base::Entry &rtree_base::get_entry(Eid e) {
  return m_entries[e.id];
}
inline rtree_base::Eid rtree_base::get_node_entry(Nid n, int idx) const {
  return m_node_entries[n.id * M + idx];
}
inline void rtree_base::set_node_entry(Nid n, int idx, Eid e) {
  m_node_entries[n.id * M + idx] = e;
}

// TODO noexcept? inline has a great impact in performance!
// Especially when used as a (statically) linked lib vs just
// compiling together the sources.
inline rtree_base::ELEMTYPE& rtree_base::rect_low_rw(const Rid &r, int dim) {
  return m_rects_low[r.id * m_dims + dim];
}
inline rtree_base::ELEMTYPE& rtree_base::rect_high_rw(const Rid &r, int dim) {
  return m_rects_high[r.id * m_dims + dim];
}
// const ref vs value? makes no difference I guess
inline const rtree_base::ELEMTYPE rtree_base::rect_low_ro(const Rid &r, int dim) const {
  return m_rects_low[r.id * m_dims + dim];
}
inline const rtree_base::ELEMTYPE rtree_base::rect_high_ro(const Rid &r, int dim) const {
  return m_rects_high[r.id * m_dims + dim];
}

inline rtree_base::ELEMTYPE rtree_base::rect_volume(Rid r) const {
  // RectVolume (TODO add also RectSphericalVolume ?)
  ELEMTYPE volume = (ELEMTYPE)1;

  for (int i = 0; i < m_dims; ++i) {
    volume *= rect_high_ro(r, i) - rect_low_ro(r, i);
  }

  ASSERT(volume >= (ELEMTYPE)0);

  return volume;
}
inline bool rtree_base::rect_contains(Rid bigger, Rid smaller) const {
  for (int index = 0; index < m_dims; ++index) {
    if (rect_low_ro(bigger, index) > rect_low_ro(smaller, index) ||
        rect_high_ro(bigger, index) < rect_high_ro(smaller, index)) {
      return false;
    }
  }
  return true;
}
inline bool rtree_base::rects_overlap(Rid a, Rid b) const {
  for (int index = 0; index < m_dims; ++index) {
    if (rect_low_ro(a, index) > rect_high_ro(b, index) ||
        rect_low_ro(b, index) > rect_high_ro(a, index)) {
      return false;
    }
  }
  return true;
}

inline void rtree_base::copy_rect(Rid src, Rid dst) {
  for (int i = 0; i < m_dims; ++i) {
    rect_low_rw(dst, i) = rect_low_ro(src, i);
    rect_high_rw(dst, i) = rect_high_ro(src, i);
  }
}
inline void rtree_base::combine_rects(Rid a, Rid b, Rid dst) {
  for (int i = 0; i < m_dims; i++) {
    rect_low_rw(dst, i) = Min(rect_low_ro(a, i), rect_low_ro(b, i));
    rect_high_rw(dst, i) = Max(rect_high_ro(a, i), rect_high_ro(b, i));
  }
}

rtree_base::rtree_base(int dimensions) : m_dims(dimensions) {
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

rtree_base::Nid rtree_base::choose_leaf(Rid r, Traversal &traversal) {
  Parent p;
  p.node = m_root_id;
  traversal.push_back(p);
  return choose_node(m_root_id, r, 0, traversal);
}

rtree_base::Nid rtree_base::choose_node(Nid n, Rid r, int height, Traversal &traversal) {
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

rtree_base::Eid rtree_base::choose_subtree(Nid n, Rid r) {
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

void rtree_base::update_entry_rect(Eid e) {
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
    rect_low_rw(r, i) = rect_low_ro(child_entry.rect_id, i);
    rect_high_rw(r, i) = rect_high_ro(child_entry.rect_id, i);
  }
  for (int i = 1; i < node.count; ++i) {
    Eid child = get_node_entry(n, i);
    combine_rects(get_entry(child).rect_id, r, r);
  }
}

void rtree_base::insert(const Vec &low, const Vec &high, Did did) {
  ++m_size;
  const Eid e = make_entry_id(); // also sets the rect
  Entry &entry = get_entry(e);
  entry.data_id = did;
  // set_data(entry.data_id, data);

  const Rid r = entry.rect_id;

  for (int i = 0; i < m_dims; ++i) {
    rect_low_rw(r, i) = low[i];
    rect_high_rw(r, i) = high[i];
  }

  m_traversal.clear();
  const Nid n = choose_leaf(entry.rect_id, m_traversal);
  const Nid nn = insert(n, e);

  adjust_tree(m_traversal, e, nn);
}

void rtree_base::reinsert_entry(Eid e) {
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

rtree_base::Nid rtree_base::insert(Nid n, Eid e) {
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

void rtree_base::plain_insert(Nid n, Eid e) {
  Node &node = get_node(n);
  ASSERT(node.count < M);
  set_node_entry(n, node.count, e);
  ++node.count;
}

rtree_base::Nid rtree_base::split_and_insert(Nid n, Eid e) {
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

  Seeds seeds = pick_seeds(entries);

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

void rtree_base::distribute_entries(Nid n, Nid nn, std::vector<Eid> entries,
                                    const Seeds &seeds) {
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

  ASSERT(entries.size() == static_cast<uint>(M - 1)); // (M+1) -2

  ELEMTYPE biggestDiff;
  uint group, entry_index = 0, betterGroup = 0;
  // existing node: group[0], new node: group[1]
  const Nid groups[2] = {n, nn};
  const Node *nodes[2] = {&node, &new_node};
  const int max_fill = M + 1 - m;
  while (!entries.empty() && node.count < max_fill &&
         new_node.count < max_fill) {
    biggestDiff = (ELEMTYPE)-1;
    for (uint i = 0; i < entries.size(); ++i) {
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

void rtree_base::distribute_entries_naive(Nid n, Nid nn,
                                          std::vector<Eid> entries) {
  uint half = entries.size() / 2 + 1;
  for (uint i = 0; i < half; ++i) {
    plain_insert(n, entries[i]);
  }
  for (uint i = half; i < entries.size(); ++i) {
    plain_insert(nn, entries[i]);
  }
}

void rtree_base::adjust_rects(const Traversal &traversal) {
  // important! no uint here (decreasing, want to reach 0)
  for (int i = static_cast<int>(traversal.size()) - 1; i >= 0; --i) {
    const Parent &parent = traversal[i];
    if (parent.entry) {
      update_entry_rect(parent.entry);
    }
  }
}

void rtree_base::adjust_tree(const Traversal &traversal, Eid e, Nid nn) {
  const bool dbg = false; // e.id == 159;
  if (dbg) {
    // cout << ">>> adjust tree " << traversal << " " << e << " nn " << nn
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

    bool is_above_new_node = level < (int)traversal.size() - 1;
    // propagating new node (split) upwards.
    if (nn && is_above_new_node) {
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
rtree_base::Seeds rtree_base::pick_seeds(const std::vector<Eid> &entries) {
  Seeds res;
  ASSERT(entries.size() == static_cast<uint>(M + 1));

  copy_rect(get_entry(entries[0]).rect_id, m_partition.cover_rect);
  for (size_t i = 0; i < entries.size(); ++i) {
    combine_rects(m_partition.cover_rect, get_entry(entries[i]).rect_id,
                  m_partition.cover_rect);
  }
  m_partition.cover_area = rect_volume(m_partition.cover_rect);

  int seed0 = 0, seed1 = 0;
  ELEMTYPE worst, waste;
  worst = -m_partition.cover_area - 1;
  for (uint i = 0; i < entries.size(); ++i) {
    m_partition.entries_areas[i] = rect_volume(get_entry(entries[i]).rect_id);
  }

  for (uint i = 0; i < entries.size() - 1; ++i) {
    for (uint j = i + 1; j < entries.size(); ++j) {
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

size_t rtree_base::size() { return m_size; }

std::vector<rtree_base::Did> rtree_base::search(const Vec &low, const Vec &high) {
  std::vector<Did> res;

  search(low, high, res);
  return res;
}

void rtree_base::search(const Vec &low, const Vec &high,
                        std::vector<Did> &results) {
  ASSERT(low.size() == static_cast<uint>(m_dims));
  ASSERT(high.size() == static_cast<uint>(m_dims));
  for (int i = 0; i < m_dims; ++i) {
    rect_low_rw(m_temp_rect, i) = low[i];
    rect_high_rw(m_temp_rect, i) = high[i];
  }
  SearchCb cb = [&results](const Did &data) {
    results.push_back(data);
    return true;
  };

  int found_count = 0;
  search(m_root_id, m_temp_rect, found_count, cb);
}

bool rtree_base::search(Nid n, Rid r, int &found_count, const SearchCb &cb) {
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
        if (!cb(entry.data_id)) {
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

int rtree_base::remove(const Vec &low, const Vec &high) {
  ASSERT(low.size() == static_cast<uint>(m_dims));
  ASSERT(high.size() == static_cast<uint>(m_dims));
  for (int i = 0; i < m_dims; ++i) {
    rect_low_rw(m_temp_rect, i) = low[i];
    rect_high_rw(m_temp_rect, i) = high[i];
  }
  int removed = 0;
  Traversals remove_traverals;
  Predicate cb = [](Did ) { return true; };
  Traversal cur_traversal;
  Parent p;
  p.node = m_root_id;
  cur_traversal.push_back(p);
  remove(m_root_id, m_temp_rect, removed, remove_traverals, cur_traversal, cb);

  // sorting remove_traverals: longer traverals first, thus leaf node
  // entries (data) removed first
  std::sort(remove_traverals.begin(), remove_traverals.end(),
            [](const Traversal &a, const Traversal &b) {
              return a.size() > b.size();
            });
  condense_tree(remove_traverals);
  m_size -= removed;
  return removed;
}

void rtree_base::remove(Nid n, Rid r, int &counter, Traversals &traversals,
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
        if (cb(entry.data_id)) {
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

bool rtree_base::remove_node_entry(Nid n, Eid e) {
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

void rtree_base::remove_node_entry(Nid n, int idx) {
  Node &node = get_node(n);
  ASSERT(idx < node.count);
  set_node_entry(n, idx, get_node_entry(n, node.count - 1));
  --node.count;
}

void rtree_base::condense_tree(const Traversals &traversals) {
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

int rtree_base::count(Nid n) {
  int counter = 0;
  count(n, counter);
  return counter;
}

void rtree_base::count(Nid n, int &counter) {
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

bool rtree_base::has_duplicate_nodes() {
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

void rtree_base::entry_to_string(Eid e, int level, int spaces, std::ostream &os) {
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
    os << " data-id=\"" << entry.data_id.id << "\"";
    // os << " data=\"" << pp(get_data(entry.data_id)) << "\"";
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

std::string rtree_base::entry_to_string(Eid e, int level, int spaces) {
  std::ostringstream os;
  entry_to_string(e, level, spaces, os);
  return os.str();
};

void rtree_base::node_to_string(Nid n, int level, int spaces, std::ostream &os) {
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

std::string rtree_base::node_to_string(Nid nid, int level, int spaces) {
  std::ostringstream os;
  node_to_string(nid, level, spaces, os);
  return os.str();
};

std::string rtree_base::to_string() {
  std::ostringstream os;
  to_string(4, os);
  return os.str();
}

void rtree_base::to_string(int spaces, std::ostream &os) {
  node_to_string(m_root_id, 0, spaces, os);
}

void rtree_base::rect_to_string(Rid rid, std::ostream &os) {
  ASSERT(rid);
  for (int i = 0; i < m_dims; i++) {
    if (i == 0) {
      os << "{";
    }
    os << rect_low_ro(rid, i);
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
    os << rect_high_rw(rid, i);
    if (i != m_dims - 1) {
      os << ", ";
    } else {
      os << "}";
    }
  }

}

std::string rtree_base::rect_to_string(Rid rid) {
  std::ostringstream os;
  rect_to_string(rid, os);
  return os.str();
}

std::ostream &operator<<(std::ostream &os, const rtree_base::Xml& xml) {
  os << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\" ?>"
     << std::endl;
  xml.tree->to_string(xml.spaces, os);
  return os;
}

rtree_base::Xml rtree_base::to_xml() {
  return Xml(this);
}

// 
} // aod
