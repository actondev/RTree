#include "./drtree3_inc.hpp"
#include "drtree/drtree3_template.hpp"

DRTREE_TEMPLATE
std::vector<std::reference_wrapper<const DATATYPE>> QUAL::search(VEC low, VEC high) {
  std::vector<std::reference_wrapper<const DATATYPE>> found;
    Callback cb = [&](const DATATYPE& data, const double *low, const double *high) {
      found.push_back(data);
    return true;
  };
  search(low, high, cb);
  return found;
}

DRTREE_TEMPLATE
int QUAL::search(VEC low, VEC high, Callback callback) {
  ASSERT(low.size() == m_dims);
  ASSERT(high.size() == m_dims);
  for (unsigned int axis = 0; axis < m_dims; ++axis) {
    rect_min_ref(m_search_rect, axis) = low.data()[axis];
    rect_max_ref(m_search_rect, axis) = high.data()[axis];
  }

  int found_count = 0;
  Search(m_root_id, m_search_rect, found_count, callback);

  return found_count;
}

DRTREE_TEMPLATE
bool QUAL::Search(Nid nid, Rid rid, int &found_count,
                         Callback callback) {
  // ASSERT(a_node);
  // ASSERT(a_node->m_level >= 0);
  // ASSERT(a_rect);
  Node& node = get_node(nid);
  Bid nbid = node.get_branch(0);
  if (node.is_internal()) {
    // This is an internal node in the tree
    for (int index = 0; index < node.count; ++index, ++nbid) {
      Branch& branch = get_branch(nbid);
      if (Overlap(rid, branch.rect_id)) {
        if (!Search(branch.child, rid, found_count,
                    callback)) {
          // The callback indicated to stop searching
          return false;
        }
      }
    }
  } else {
    // This is a leaf node
    for (int index = 0; index < node.count; ++index, ++nbid) {
      Branch& branch = get_branch(nbid);
      if (Overlap(rid, branch.rect_id)) {
        const DATATYPE &data = branch_data(nbid);
        ++found_count;

        if (!callback(data, rect_min(branch.rect_id), rect_max(branch.rect_id))) {
          return false; // Don't continue searching
        }
      }
    }
  }

  return true; // Continue searching
}


DRTREE_TEMPLATE
void QUAL::push(const ELEMTYPE *low, const ELEMTYPE *high,
                         const DATATYPE &data) {
#ifdef _DEBUG
  for (int index = 0; index < dims; ++index) {
    ASSERT(a_min[index] <= a_max[index]);
  }
#endif //_DEBUG

  // this can be reused (m_insert_branch)
  Bid bid = make_branch_id();
  set_branch_data(bid, data);
  Branch& branch = get_branch(bid);
  
  for (unsigned int axis = 0; axis < m_dims; ++axis) {
    rect_min_ref(branch.rect_id, axis) = low[axis];
    rect_max_ref(branch.rect_id, axis) = high[axis];
  }

  InsertRect(bid, m_root_id, 0);
  m_size++;
}


DRTREE_TEMPLATE
int QUAL::PickBranch(Rid rid, Nid nid) {
  // ASSERT(a_rect && a_node);

  bool firstTime = true;
  ELEMTYPE increase;
  ELEMTYPE bestIncr = (ELEMTYPE)-1;
  ELEMTYPE area;
  ELEMTYPE bestArea;
  int best = 0;

  Node& node = get_node(nid);
  Bid bid = node.branch0;
  for (int index = 0; index < node.count; ++index, ++bid) {
    Branch& branch = get_branch(bid);
    Rid cur_rect = branch.rect_id;
    area = CalcRectVolume(cur_rect);
    combine_rects(m_pick_branch_rect, rid, cur_rect);
    increase = CalcRectVolume(m_pick_branch_rect) - area;
    if ((increase < bestIncr) || firstTime) {
      best = index;
      bestArea = area;
      bestIncr = increase;
      firstTime = false;
    } else if ((increase == bestIncr) && (area < bestArea)) {
      best = index;
      bestArea = area;
      bestIncr = increase;
    }
  }
  return best;
}

DRTREE_TEMPLATE
bool QUAL::InsertRect(Bid bid, Nid& root_id, int level) {
  ASSERT(root_id);
  Node& root = get_node(root_id);
  ASSERT(level >= 0 && level <= root.level);
#ifdef _DEBUG
  for (int index = 0; index < dims; ++index) {
    ASSERT(rect_min_ref(a_branch.m_rect, index) <=
           a_branch.m_rect.m_max[index]);
  }
#endif //_DEBUG

  Nid new_node;
  Nid new_root_id; // falsy

  if (InsertRectRec(bid, root_id, new_node, level)) // Root split
  {
    // Grow tree taller and new root
    // Node *newRoot = AllocNode();
    new_root_id = make_node_id();
    Node &new_root = get_node(new_root_id);
    new_root.level = root.level + 1;

    Nid falsy_nid; // TODO passing NULL to addbranch in original

    Branch &insert_rect_branch = get_branch(m_insert_rect_branch);
    // add old root node as a child of the new root
    node_cover(insert_rect_branch.rect_id, root_id);
    insert_rect_branch.child = root_id;
    AddBranch(m_insert_rect_branch, new_root_id, falsy_nid);

    // add the split node as a child of the new root
    node_cover(insert_rect_branch.rect_id, new_node);
    insert_rect_branch.child = new_node;
    AddBranch(m_insert_rect_branch, new_root_id, falsy_nid);

    // set the new root as the root node
    // *a_root = newRoot;
    root_id = new_root_id;
  }
  // truthy if new_root_id = make_node_id();
  return new_root_id;
}

DRTREE_TEMPLATE
bool QUAL::InsertRectRec(Bid bid, Nid nid,
                                Nid& new_nid, int level) {
  // ASSERT(a_node && a_newNode);
  // ASSERT(a_level >= 0);
  // ASSERT(a_level <= a_node->m_level);

  // recurse until we reach the correct level for the new record. data records
  // will always be called with a_level == 0 (leaf)
  Node& node = get_node(nid);
  Branch& branch = get_branch(bid);
  if (node.level > level) {
    // Still above level for insertion, go down tree recursively
    // Node *otherNode;
    Nid other_nid;

    // find the optimal branch for this record
    int index = PickBranch(branch.rect_id, nid);
    Bid node_bid = node.get_branch(index);
    Branch& node_branch = get_branch(node_bid);

    // recursively insert this record into the picked branch
    bool childWasSplit = InsertRectRec(bid, node_branch.child, other_nid, level);

    if (!childWasSplit) {
      // Child was not split. Merge the bounding box of the new record with the
      // existing bounding box
      combine_rects(node_branch.rect_id, branch.rect_id,
                    node_branch.rect_id);
      return false;
    } else {
      // Child was split. The old branches are now re-partitioned to two nodes
      // so we have to re-calculate the bounding boxes of each node
      node_cover(node_branch.rect_id,
                node_branch.child);

      Branch& insert_rect_rec_branch = get_branch(m_insert_rect_rec_branch);
      insert_rect_rec_branch.child = other_nid;
      node_cover(insert_rect_rec_branch.rect_id, other_nid);

      // The old node is already a child of a_node. Now add the newly-created
      // node to a_node as well. a_node might be split because of that.
      return AddBranch(m_insert_rect_rec_branch, nid, new_nid);
    }
  } else if (node.level == level) {
    // We have reached level for insertion. Add rect, split if necessary
    return AddBranch(bid, nid, new_nid);
  } else {
    // Should never occur
    ASSERT(0);
    return false;
  }
}

DRTREE_TEMPLATE
bool QUAL::AddBranch(Bid bid, Nid nid, Nid& new_nid) {
  Node &node = get_node(nid);
  if (node.count < MAXNODES) // Split won't be necessary
  {
    Branch& node_branch = get_branch(node.get_branch(node.count));
    Branch& branch = get_branch(bid);
    node_branch = branch; // copy assignment (rect & child?)
    // a_node->m_branch[a_node->m_count] = *a_branch;
    ++node.count;

    return false;
  } else {
    // ASSERT(a_newNode);

    // TODO
    SplitNode(nid, bid, new_nid);
    return true;
  }
}

DRTREE_TEMPLATE
void QUAL::SplitNode(Nid nid, Bid bid, Nid &new_nid) {
  ASSERT(nid);
  ASSERT(bid);

  // Load all the branches into a buffer, initialize old node
  // GetBranches(a_node, a_branch, &m_split_parition_vars);

  // Find partition
  ChoosePartition(m_partition_vars, MINNODES);

  // Create a new node to hold (about) half of the branches
  new_nid = make_node_id();
  Node& new_node = get_node(new_nid);
  Node& node = get_node(nid);
  new_node.level = node.level;

  // Put branches from buffer into 2 nodes according to the chosen partition
  node.count = 0;
  // LoadNodes(a_node, *a_newNode, &m_split_parition_vars);

  ASSERT((node.count + new_node.count) ==
         m_partition_vars.m_total);
}

DRTREE_TEMPLATE
void QUAL::ChoosePartition(PartitionVars& a_parVars, int a_minFill) {
  ELEMTYPE biggestDiff;
  int group, chosen = 0, betterGroup = 0;

  InitParVars(a_parVars, a_parVars.m_branchCount, a_minFill);
  PickSeeds(a_parVars);

  while (
      ((a_parVars.m_count[0] + a_parVars.m_count[1]) < a_parVars.m_total) &&
      (a_parVars.m_count[0] < (a_parVars.m_total - a_parVars.m_minFill)) &&
      (a_parVars.m_count[1] < (a_parVars.m_total - a_parVars.m_minFill))) {
    biggestDiff = (ELEMTYPE)-1;
    Bid bid = a_parVars.m_branchBuf[0];
    for (int index = 0; index < a_parVars.m_total; ++index, ++bid) {
      Branch& branch = get_branch(bid);
      if (PartitionVars::NOT_TAKEN == a_parVars.m_partition[index]) {
        // Rect *curRect = &a_parVars->m_branchBuf[index].m_rect;
        Rid cur_rect = branch.rect_id;
        combine_rects(m_choose_partition_rect0, cur_rect,
                      a_parVars.m_cover[0]);
        combine_rects(m_choose_partition_rect1, cur_rect,
                      a_parVars.m_cover[1]);

        ELEMTYPE growth0 =
            CalcRectVolume(m_choose_partition_rect0) - a_parVars.m_area[0];
        ELEMTYPE growth1 =
            CalcRectVolume(m_choose_partition_rect1) - a_parVars.m_area[1];
        ELEMTYPE diff = growth1 - growth0;
        if (diff >= 0) {
          group = 0;
        } else {
          group = 1;
          diff = -diff;
        }

        if (diff > biggestDiff) {
          biggestDiff = diff;
          chosen = index;
          betterGroup = group;
        } else if ((diff == biggestDiff) && (a_parVars.m_count[group] <
                                             a_parVars.m_count[betterGroup])) {
          chosen = index;
          betterGroup = group;
        }
      }
    }
    Classify(chosen, betterGroup, a_parVars);
  }

  // If one group too full, put remaining rects in the other
  if ((a_parVars.m_count[0] + a_parVars.m_count[1]) < a_parVars.m_total) {
    if (a_parVars.m_count[0] >= a_parVars.m_total - a_parVars.m_minFill) {
      group = 1;
    } else {
      group = 0;
    }
    for (int index = 0; index < a_parVars.m_total; ++index) {
      if (PartitionVars::NOT_TAKEN == a_parVars.m_partition[index]) {
        Classify(index, group, a_parVars);
      }
    }
  }

  ASSERT((a_parVars.m_count[0] + a_parVars.m_count[1]) == a_parVars.m_total);
  ASSERT((a_parVars.m_count[0] >= a_parVars.m_minFill) &&
         (a_parVars.m_count[1] >= a_parVars.m_minFill));
}

DRTREE_TEMPLATE
void QUAL::InitParVars(PartitionVars& a_parVars, int a_maxRects,
                              int a_minFill) {
  a_parVars.m_count[0] = a_parVars.m_count[1] = 0;
  a_parVars.m_area[0] = a_parVars.m_area[1] = (ELEMTYPE)0;
  a_parVars.m_total = a_maxRects;
  a_parVars.m_minFill = a_minFill;
  for (int index = 0; index < a_maxRects; ++index) {
    a_parVars.m_partition[index] = PartitionVars::NOT_TAKEN;
  }
}

DRTREE_TEMPLATE
void QUAL::node_cover(Rid dst, Nid nid) {
  ASSERT(nid);
  // *dst = a_node->m_branch[0].m_rect;
  Node& node = get_node(nid);
  Bid bid0 = node.get_branch(0);
  Branch& branch0 = get_branch(bid0);
  copy_rect(branch0.rect_id, dst);
  Bid bid = bid0;
  for (int index = 1; index < node.count; ++index, ++bid) {
    // we just increase breanch id
    Branch& branch = get_branch(bid);
    combine_rects(dst, dst, branch.rect_id);
  }
}

DRTREE_TEMPLATE
void QUAL::copy_rect(Rid src, Rid dst) {
  for (auto i = 0; i < m_dims; i++) {
    rect_min_ref(dst, i) = rect_min_ref(src, i);
    rect_max_ref(dst, i) = rect_max_ref(src, i);
  }
}

DRTREE_TEMPLATE
void QUAL::combine_rects(Rid dst, Rid a, Rid b) {
  ASSERT(dst && a && b);

  for (unsigned int index = 0; index < m_dims; index++) {
    rect_min_ref(dst, index) = Min(rect_min_ref(a, index), rect_min_ref(b, index));
    rect_max_ref(dst, index) = Max(rect_max_ref(a, index), rect_max_ref(b, index));
  }
}


// Calculate the n-dimensional volume of a rectangle
DRTREE_TEMPLATE
ELEMTYPE QUAL::RectVolume(Rid rid) {
  return 0;
}

// The exact volume of the bounding sphere for the given Rect
DRTREE_TEMPLATE
ELEMTYPE QUAL::RectSphericalVolume(Rid rid) {
  ELEMTYPE sumOfSquares = (ELEMTYPE)0;
  ELEMTYPE radius;
  for (unsigned int index = 0; index < m_dims; ++index) {
    const ELEMTYPE halfExtent =
        (rect_max_ref(rid, index) - rect_min_ref(rid, index)) *
        (ELEMTYPE)0.5;

    sumOfSquares += halfExtent * halfExtent;
  }

  radius = (ELEMTYPE)sqrt(sumOfSquares);

  // Pow maybe slow, so test for common dims like 2,3 and just use x*x, x*x*x.
  switch(m_dims) {
  case 2:
    return (radius * radius * m_unitSphereVolume);
  case 3:
    return (radius * radius * radius * m_unitSphereVolume);
    default:
      return (ELEMTYPE)(pow(radius, m_dims) * m_unitSphereVolume);
  }
}

// Use one of the methods to calculate retangle volume
DRTREE_TEMPLATE
ELEMTYPE QUAL::CalcRectVolume(Rid rid) {
#ifdef RTREE_USE_SPHERICAL_VOLUME
  return RectSphericalVolume(rid); // Slower but helps certain merge cases
#else                                 // RTREE_USE_SPHERICAL_VOLUME
  return RectVolume(rid); // Faster but can cause poor merges
#endif                                // RTREE_USE_SPHERICAL_VOLUME
}

// Decide whether two rectangles overlap.
DRTREE_TEMPLATE
bool QUAL::Overlap(Rid a, Rid b) {
  ASSERT(a && b);

  for (unsigned int index = 0; index < m_dims; ++index) {
    if (rect_min_ref(a, index) > rect_max_ref(b, index) ||
        rect_min_ref(b, index) > rect_max_ref(a, index)) {
      return false;
    }
  }
  return true;
}
