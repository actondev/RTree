#pragma once
#include "./drtree3_inc.hpp"

#define ASSERT assert

#define pp(x) easyprint::stringify(x)

template <class ELEMTYPE> struct TPartitionVars {
  enum { NOT_TAKEN = -1 }; // indicates that position
  int m_total;
  int m_minFill;
  int m_branchCount;

  ELEMTYPE m_coverSplitArea;
  Rid m_coverSplit;
  ELEMTYPE m_area[2];
  Rid m_cover[2];
  int m_count[2];

  std::vector<int> m_partition;
  std::vector<Bid> m_branchBuf;

  TPartitionVars() = delete;
  TPartitionVars(size_t size) // size: MAXNODES + 1
      : m_partition(size),
        m_branchBuf(size)
  { }
};

template <class DATATYPE, class ELEMTYPE = double> class drtree3 {
  static_assert(std::numeric_limits<ELEMTYPE>::is_iec559,
                "'COORD_TYPE' accepts floating-point types only");

  typedef std::function<bool(const DATATYPE &, const ELEMTYPE *,
                             const ELEMTYPE *)>
      Callback;
  typedef std::function<bool(const DATATYPE &)>
  RemovePredicate;
  typedef TPartitionVars<ELEMTYPE> PartitionVars;
  typedef std::vector<ELEMTYPE> VEC;

private:
  const int MAXNODES = 8;
  const int MINNODES = MAXNODES / 2;
  unsigned int m_dims; // set by the constructor
  std::vector<ELEMTYPE> m_rects_min;
  std::vector<ELEMTYPE> m_rects_max;

  std::vector<Branch> m_branches;
  std::vector<DATATYPE> m_data_store;

  std::vector<Node> m_nodes;

  ELEMTYPE m_unitSphereVolume;

  PartitionVars m_partition_vars;
  id_t m_rects_count = 0;
  id_t m_branches_count = 0;
  id_t m_nodes_count = 0;

  Bid m_insert_branch;
  Bid m_insert_rect_rec_branch;
  Bid m_insert_rect_branch;
  Rid m_remove_rect;
  Rid m_search_rect;
  //
  Rid m_pick_branch_rect;
  Rid m_pick_seeds_rect;
  Rid m_choose_partition_rect0;
  Rid m_choose_partition_rect1;

  size_t m_size = 0;

  drtree3() = delete;

  Rid make_rect_id() {
    Rid res{m_rects_count++};
    m_rects_min.resize(m_rects_count * m_dims, 0);
    m_rects_max.resize(m_rects_count * m_dims, 0);

    return res;
  }

  Bid make_branch_id() {
    Bid bid{m_branches_count++};
    m_branches.resize(m_branches_count);
    Branch &branch = m_branches[bid.id];
    branch.rect_id = make_rect_id();

    return bid;
  }

  Nid m_root_id;

  Branch &get_branch(Bid id) { return m_branches[id.id]; }

  // each node has MAXNODES children

  Nid make_node_id(bool attach_branches = true) {
    Nid id{m_nodes_count++};
    m_nodes.resize(m_nodes_count);

    if (attach_branches) {
      node_make_branches(id);
    }
    return id;
  }
  Node &get_node(Nid id) { return m_nodes[id.id]; }
  void node_make_branches(Nid id) {
    Node &node = get_node(id);
    node.branch0 = make_branch_id();
    for (int i = 0; i < MAXNODES - 1; i++) {
      make_branch_id();
    }
  }

  Branch &node_get_branch(Node &node, int idx) {
    ASSERT(idx < MAXNODES);
    Bid bid = node.branch0 + idx;
    return get_branch(bid);
  }

  Node &make_node(bool attach_branches = true) {
    Nid nid = make_node_id(attach_branches);

    return m_nodes[nid.id];
  }

  Did store_data(const DATATYPE& data) {
    m_data_store.push_back(data);
    Did did{m_data_store.size() - 1};
    return did;
  }

  const DATATYPE &branch_data(Bid bid) {
    Branch& br = get_branch(bid);
    Did did = br.data_id;
    ASSERT(did);
    return m_data_store[did.id];
  }

  ELEMTYPE &rect_min_ref(Rid rid, int dim) {
    return m_rects_min[rid.id * m_dims + dim];
  }
  ELEMTYPE &rect_max_ref(Rid rid, int dim) {
    return m_rects_max[rid.id * m_dims + dim];
  }
  ELEMTYPE *rect_min(Rid id) { return m_rects_min.data() + id.id * m_dims; }
  ELEMTYPE *rect_max(Rid id) { return m_rects_max.data() + id.id * m_dims; }

  // logic

  bool InsertRect(Bid, Nid&, int level);
  bool InsertRectRec(Bid branch, Nid node, Nid &new_node, int level);
  bool AddBranch(Bid, Nid node, Nid& new_node);
  void SplitNode(Nid, Bid, Nid& new_nid);
  void ChoosePartition(PartitionVars&, int min_fill);
  void InitParVars(PartitionVars&, int a_maxRects, int a_minFill);
  void PickSeeds(PartitionVars& a_parVars);
  void Classify(int a_index, int a_group, PartitionVars&);
  void GetBranches(Nid, Bid, PartitionVars&);
  void LoadNodes(Nid a, Nid b, PartitionVars&);

  int PickBranch(Rid rid, Nid nid);
  void node_cover(Rid dst, Nid nid);
  void combine_rects(Rid dst, Rid a, Rid b);
  void copy_rect(Rid dst, Rid src);
  void copy_branch(Bid dst, Bid src);
  bool Overlap(Rid, Rid);
  ELEMTYPE RectSphericalVolume(Rid);
  ELEMTYPE RectVolume(Rid);
  ELEMTYPE CalcRectVolume(Rid);
  bool Search(Nid, Rid, int &found_count, Callback);

  bool RemoveRect(Nid &a_root, Rid &a_rect, int &a_foundCount);
  bool RemoveRectRec(Nid nid, Rid a_rect, int &a_removedCount, std::vector<Nid> &a_listNode);

public:
  drtree3(unsigned int dims)
      : m_dims{dims},
        m_partition_vars(MAXNODES + 1),    //
        m_insert_branch{make_branch_id()}, //
        m_insert_rect_rec_branch{make_branch_id()},
        m_insert_rect_branch{make_branch_id()}, //
        m_remove_rect{make_rect_id()},          //
        m_search_rect{make_rect_id()},           //
        m_pick_branch_rect{make_rect_id()},     //
        m_pick_seeds_rect{make_rect_id()},
        m_choose_partition_rect0{make_rect_id()},
        m_choose_partition_rect1{make_rect_id()}
  {
    m_root_id = make_node_id();
    for (int i = 0; i < MAXNODES + 1; i++) {
      m_partition_vars.m_branchBuf[i] = make_branch_id();
    }
    m_partition_vars.m_cover[0] = make_rect_id();
    m_partition_vars.m_cover[1] = make_rect_id();
    m_partition_vars.m_coverSplit = make_rect_id();

    // Precomputed volumes of the unit spheres for the first few dimensions
    const float UNIT_SPHERE_VOLUMES[] = {
        0.000000f, 2.000000f, 3.141593f, // Dimension  0,1,2
        4.188790f, 4.934802f, 5.263789f, // Dimension  3,4,5
        5.167713f, 4.724766f, 4.058712f, // Dimension  6,7,8
        3.298509f, 2.550164f, 1.884104f, // Dimension  9,10,11
        1.335263f, 0.910629f, 0.599265f, // Dimension  12,13,14
        0.381443f, 0.235331f, 0.140981f, // Dimension  15,16,17
        0.082146f, 0.046622f, 0.025807f, // Dimension  18,19,20
    };

    m_unitSphereVolume = (ELEMTYPE)UNIT_SPHERE_VOLUMES[m_dims];
  }

  void push(const ELEMTYPE *a_min, const ELEMTYPE *a_max,
            const DATATYPE &a_dataId);
  size_t size() { return m_size; }

  int search(VEC low, VEC high, Callback cb);
  std::vector<DATATYPE> search(VEC low, VEC high);
  void search(VEC low, VEC high, std::vector<DATATYPE>& found);
  int remove(VEC low, VEC high);
};

#include "./drtree3_impl.hpp"

#undef pp
