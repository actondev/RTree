#include "./drtree3_inc.hpp"
#include "drtree/drtree3_template.hpp"

DRTREE_TEMPLATE
std::vector<DATATYPE> QUAL::search(VEC low, VEC high) {
  std::vector<DATATYPE> found;
    Callback cb = [&](const DATATYPE& data, const double *low, const double *high) {
      found.push_back(data);
    return true;
  };
  search(low, high, cb);
  return found;
}

DRTREE_TEMPLATE
void QUAL::search(VEC low, VEC high, std::vector<DATATYPE>& found) {
  found.clear();
  Callback cb = [&](const DATATYPE& data, const double *low, const double *high) {
    found.push_back(data);
    return true;
  };
  search(low, high, cb);
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
  if (node.is_internal()) {
    // This is an internal node in the tree
    for (int index = 0; index < node.count; ++index) {
      Bid nbid = get_node_bid(nid, index);
      Branch& branch = get_branch(nbid);
      if (Overlap(rid, branch.rect_id)) {
        if (!Search(branch.child_id, rid, found_count,
                    callback)) {
          // The callback indicated to stop searching
          return false;
        }
      }
    }
  } else {
    // This is a leaf node
    for (int index = 0; index < node.count; ++index) {
      Bid nbid = get_node_bid(nid, index);
      Branch& branch = get_branch(nbid);
      if (Overlap(rid, branch.rect_id)) {
        const DATATYPE &data = branch_data(nbid);
        // cout << "search res: bid " << nbid << " data " << pp(data) << endl;
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
int QUAL::remove(VEC low, VEC high) {
  ASSERT(low.size() == m_dims);
  ASSERT(high.size() == m_dims);
  for (unsigned int axis = 0; axis < m_dims; ++axis) {
    rect_min_ref(m_remove_rect, axis) = low.data()[axis];
    rect_max_ref(m_remove_rect, axis) = high.data()[axis];
  }

  int removed_count = 0;

  RemoveRect(m_root_id, m_remove_rect, removed_count);
  m_size -= removed_count;
  return removed_count;
}

DRTREE_TEMPLATE
bool QUAL::RemoveRect(Nid nid, Rid rid, int &found_count) {
  ASSERT(nid && rid);

  std::vector<Nid> reinsert_list;

  if (!RemoveRectRec(nid, rid, found_count, reinsert_list)) {
    return false; // nothing removed
  }

  // removed
  while (!reinsert_list.empty()) {
    Nid temp_nid = reinsert_list.back();
    // cout << "reinsert list " << temp_nid << endl;
    reinsert_list.pop_back();
    Node& temp_node = get_node(temp_nid);
    for (int index = 0; index < temp_node.count; ++index) {
      Bid node_bid = get_node_bid(temp_nid, index);
      InsertRect(node_bid, nid, temp_node.level);
    }
    // TODO free node
    // FreeNode(remLNode->m_node);
    // FreeListNode(remLNode);
  }

  Node& root = get_node(nid);
  if (root.count == 1 && root.is_internal()) {
    // TODO handle this, free node
    // Node *tempNode = (*a_root)->m_branch[0].m_child;
    
    // ASSERT(tempNode);
    // FreeNode(*a_root);
    // *a_root = tempNode;
  }
  return true;
}

DRTREE_TEMPLATE
bool QUAL::RemoveRectRec(Nid nid, Rid rid, int &removed_count, std::vector<Nid> &reinsert) {
  ASSERT(nid && rid);
  Node& node = get_node(nid);
  ASSERT(node.level >= 0);
  if (node.is_internal()) {
    // not a leaf node
    bool removed = false;
    for (int index = 0; index < node.count; ++index) {
      Bid node_bid = get_node_bid(nid, index);
      Branch& branch = get_branch(node_bid);
      Rid branch_rid = branch.rect_id;
      Nid branch_child_id = branch.child_id;
      if (Overlap(rid, branch_rid)) {
        if (RemoveRectRec(branch.child_id, rid,
                          removed_count, reinsert)) {
          Node& branch_child = get_node(branch_child_id);
          if (branch_child.count >= MINNODES) {
            // child removed, just resize parent rect
            node_cover(branch_rid, branch_child_id);
            // node_cover(&a_node->m_branch[index].m_rect,
            //            a_node->m_branch[index].m_child);
          } else {
            // child removed, not enough entries in node, eliminate node
            // TODO
            // ReInsert(a_node->m_branch[index].m_child, reinsert);
            reinsert.push_back(branch_child_id);
            DisconnectBranch(nid, index);
            // NB: Before remove refactor this was returning
            index--; // have to revisit same index, as now it's swapped with the
                     // last item
          }
          removed = true;
        }
      }
    }
    return removed;
  } else {
    // A leaf node
    bool removed = false;
        
    for (int index = 0; index < node.count; ++index) {
      Bid node_bid = get_node_bid(nid, index);
      Branch &branch = get_branch(node_bid);
      if (Overlap(rid, branch.rect_id)) {
        // if (predicate(branch_data(branch))) {
          removed = true;
          // TODO
          DisconnectBranch(nid, index);
          // NB: Before remove refactor this was returning
          index--; // have to revisit same index, as now it's swapped with the
                   // last item
          removed_count++;
        // }
      }
    }
    return removed;
  }
}

DRTREE_TEMPLATE
void QUAL::DisconnectBranch(Nid nid, int index) {
  ASSERT(nid);
  Node& node = get_node(nid);
  ASSERT((index >= 0) && (index < MAXNODES));
  ASSERT(node.count > 0);

  // Remove element by swapping with the last element to prevent gaps in array
  set_node_bid(nid, index, get_node_bid(nid, node.count-1));

  --node.count;
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
  Bid bid = m_insert_branch;
  Branch& branch = get_branch(bid);
  Did did = store_data(data);
  branch.data_id = did;
  // cout << "Insert: bid id " << bid << " has data " << pp(data) << endl;
  
  for (unsigned int axis = 0; axis < m_dims; ++axis) {
    rect_min_ref(branch.rect_id, axis) = low[axis];
    rect_max_ref(branch.rect_id, axis) = high[axis];
  }

  InsertRect(bid, m_root_id, 0);
  m_size++;
}


DRTREE_TEMPLATE
int QUAL::PickBranch(Rid rid, Nid nid) {
  ASSERT(rid && nid);

  bool firstTime = true;
  ELEMTYPE increase;
  ELEMTYPE bestIncr = (ELEMTYPE)-1;
  ELEMTYPE area;
  ELEMTYPE bestArea;
  int best = 0;

  Node& node = get_node(nid);
  for (int index = 0; index < node.count; ++index) {
    Bid bid = get_node_bid(nid, index);
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
    insert_rect_branch.child_id = root_id;
    AddBranch(m_insert_rect_branch, new_root_id, falsy_nid);

    // add the split node as a child of the new root
    node_cover(insert_rect_branch.rect_id, new_node);
    insert_rect_branch.child_id = new_node;
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
  ASSERT(nid);
  // ASSERT(new_nid);
  ASSERT(level >= 0);

  // recurse until we reach the correct level for the new record. data records
  // will always be called with a_level == 0 (leaf)
  Node& node = get_node(nid);
  ASSERT(level <= node.level);
  Branch& branch = get_branch(bid);
  Rid rid = branch.rect_id;
  if (node.level > level) {
    // Still above level for insertion, go down tree recursively
    // Node *otherNode;
    Nid other_nid;

    // find the optimal branch for this record
    int index = PickBranch(rid, nid);
    Bid node_bid = get_node_bid(nid, index);
    Branch& node_branch = get_branch(node_bid);
    Nid nb_child = node_branch.child_id;
    Rid nb_rect = node_branch.rect_id;
    // recursively insert this record into the picked branch
    bool childWasSplit = InsertRectRec(bid, nb_child, other_nid, level);

    if (!childWasSplit) {
      // Child was not split. Merge the bounding box of the new record with the
      // existing bounding box
      combine_rects(nb_rect, rid, nb_rect);
      return false;
    } else {
      // Child was split. The old branches are now re-partitioned to two nodes
      // so we have to re-calculate the bounding boxes of each node

      node_cover(nb_rect, nb_child);

      Branch& insert_rect_rec_branch = get_branch(m_insert_rect_rec_branch);
      insert_rect_rec_branch.child_id = other_nid;
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
    // TODO should addbranch perhaps just return a bid? (which branch was chosen)
    Bid node_insert_bid = get_node_bid(nid, node.count);
    // cout << "add branch copying: " << bid << " => node" << nid << "branch#" << node.count << " bid " << node_insert_bid << endl;
    copy_branch(bid, node_insert_bid);
    ++node.count;
    return false;
  } else {
    SplitNode(nid, bid, new_nid);
    return true;
  }
}

DRTREE_TEMPLATE
void QUAL::SplitNode(Nid nid, Bid bid, Nid &new_nid) {
  ASSERT(nid);
  ASSERT(bid);
  // if(!new_nid) {
  //   new_nid = make_node_id();
  // }

  // Load all the branches into a buffer, initialize old node
  GetBranches(nid, bid, m_partition_vars);

  // Find partition
  ChoosePartition(m_partition_vars, MINNODES);

  // Create a new node to hold (about) half of the branches
  new_nid = make_node_id();
  Node& new_node = get_node(new_nid);
  Node& node = get_node(nid);
  new_node.level = node.level;

  // Put branches from buffer into 2 nodes according to the chosen partition
  node.count = 0;
  LoadNodes(nid, new_nid, m_partition_vars);

  ASSERT((node.count + new_node.count) ==
         m_partition_vars.m_total);
}

DRTREE_TEMPLATE
void QUAL::LoadNodes(Nid a, Nid b,
                            PartitionVars &a_parVars) {
  ASSERT(a);
  ASSERT(b);

  for (int index = 0; index < a_parVars.m_total; ++index) {
    ASSERT(a_parVars.m_partition[index] == 0 ||
           a_parVars.m_partition[index] == 1);

    int targetNodeIndex = a_parVars.m_partition[index];
    Nid targetNodes[2] = {a, b};

    // It is assured that AddBranch here will not cause a node split.
    Nid nullnid;
    bool nodeWasSplit = AddBranch(a_parVars.m_branchBuf[index],
                                  targetNodes[targetNodeIndex], nullnid);
    ASSERT(!nodeWasSplit);
  }
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
void QUAL::PickSeeds(PartitionVars& a_parVars) {
  int seed0 = 0, seed1 = 0;
  ELEMTYPE worst, waste;
  ELEMTYPE area[MAXNODES + 1];

  for (int index = 0; index < a_parVars.m_total; ++index) {
    Branch& branch = get_branch(a_parVars.m_branchBuf[index]);
    area[index] = CalcRectVolume(branch.rect_id);
  }

  worst = -a_parVars.m_coverSplitArea - 1;

  for (int indexA = 0; indexA < a_parVars.m_total - 1; ++indexA) {
    for (int indexB = indexA + 1; indexB < a_parVars.m_total; ++indexB) {
      Branch& branch_a = get_branch(a_parVars.m_branchBuf[indexA]);
      Branch& branch_b = get_branch(a_parVars.m_branchBuf[indexB]);
      combine_rects(m_pick_seeds_rect, branch_a.rect_id,
                    branch_b.rect_id);
      waste = CalcRectVolume(m_pick_seeds_rect) - area[indexA] - area[indexB];
      if (waste > worst) {
        worst = waste;
        seed0 = indexA;
        seed1 = indexB;
      }
    }
  }

  Classify(seed0, 0, a_parVars);
  Classify(seed1, 1, a_parVars);
}

DRTREE_TEMPLATE
void QUAL::Classify(int a_index, int a_group, PartitionVars& a_parVars) {
  ASSERT(PartitionVars::NOT_TAKEN == a_parVars.m_partition[a_index]);

  a_parVars.m_partition[a_index] = a_group;
  Branch& branch = get_branch(a_parVars.m_branchBuf[a_index]);
  // Calculate combined rect
  if (a_parVars.m_count[a_group] == 0) {
    // a_parVars.m_cover[a_group] = a_parVars.m_branchBuf[a_index].m_rect;
    copy_rect(branch.rect_id,
              a_parVars.m_cover[a_group]);
  } else {
    combine_rects(a_parVars.m_cover[a_group], a_parVars.m_cover[a_group],
                  branch.rect_id);
  }

  // Calculate volume of combined rect
  a_parVars.m_area[a_group] = CalcRectVolume(a_parVars.m_cover[a_group]);

  ++a_parVars.m_count[a_group];
}

DRTREE_TEMPLATE
void QUAL::GetBranches(Nid nid, Bid bid, PartitionVars &a_parVars) {
  ASSERT(nid);
  ASSERT(bid);

  Node& node = get_node(nid);
  ASSERT(node.count == MAXNODES);

  // Load the branch buffer
  // Bid nbid = node.get_branch(0);
  for (int index = 0; index < MAXNODES; ++index) {
    // Branch& node_branch = get_branch(nbid);
    // a_parVars.m_branchBuf[index] = a_node.m_branch[index];
    copy_branch(get_node_bid(nid, index), a_parVars.m_branchBuf[index]);
  }
  // a_parVars.m_branchBuf[MAXNODES] = *a_branch;
  copy_branch(bid, a_parVars.m_branchBuf[MAXNODES]);
  a_parVars.m_branchCount = MAXNODES + 1;

  // Calculate rect containing all in the set
  // Note: without implementing Rect copy assignment,
  // when modifying m_coverSplit, it would also modify m_branchBuf[0].m_rect
  // a_parVars.m_coverSplit = a_parVars.m_branchBuf[0].m_rect;
  copy_rect(get_branch(a_parVars.m_branchBuf[0]).rect_id, a_parVars.m_coverSplit);
  for (int index = 1; index < MAXNODES + 1; ++index) {
    combine_rects(a_parVars.m_coverSplit, a_parVars.m_coverSplit,
                  get_branch(a_parVars.m_branchBuf[index]).rect_id);
  }
  a_parVars.m_coverSplitArea = CalcRectVolume(a_parVars.m_coverSplit);
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
  Node& node = get_node(nid);
  Bid bid0 = get_node_bid(nid, 0);
  Branch& branch0 = get_branch(bid0);
  copy_rect(branch0.rect_id, dst);
  for (int index = 1; index < node.count; ++index) {
    Bid bid = get_node_bid(nid, index);
    // we just increase breanch id
    Branch& branch = get_branch(bid);
    combine_rects(dst, dst, branch.rect_id);
  }
}

DRTREE_TEMPLATE
void QUAL::copy_rect(Rid src, Rid dst) {
  // callgrind says std::copy is slower?
  // const ELEMTYPE* min_src = rect_min(src);
  // const ELEMTYPE* max_src = rect_max(src);
  // std::copy(min_src, min_src+m_dims, rect_min(dst));
  // std::copy(max_src, max_src+m_dims, rect_max(dst));
  for (auto i = 0; i < m_dims; ++i) {
    rect_min_ref(dst, i) = rect_min_ref(src, i);
    rect_max_ref(dst, i) = rect_max_ref(src, i);
  }
}

DRTREE_TEMPLATE
void QUAL::copy_branch(Bid src, Bid dst) {
  ASSERT(src);
  ASSERT(dst);
  Branch& src_branch = get_branch(src);
  Branch& dst_branch = get_branch(dst);
  dst_branch.child_id = src_branch.child_id;
  dst_branch.data_id = src_branch.data_id;
  copy_rect(src_branch.rect_id, dst_branch.rect_id);
  // cout << "copy branch: " << src << "->" << dst << ": " << pp(branch_data(src)) << endl;
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
  ELEMTYPE volume = (ELEMTYPE)1;

  for (int index = 0; index < m_dims; ++index) {
    volume *= rect_max_ref(rid, index) - rect_min_ref(rid, index);
  }

  ASSERT(volume >= (ELEMTYPE)0);

  return volume;
}

// The exact volume of the bounding sphere for the given Rect
DRTREE_TEMPLATE
ELEMTYPE QUAL::RectSphericalVolume(Rid rid) {
  ELEMTYPE sumOfSquares = (ELEMTYPE)0;
  ELEMTYPE radius;
  for (unsigned int index = 0; index < m_dims; ++index) {
    const ELEMTYPE halfExtent =
        ((ELEMTYPE)rect_max_ref(rid, index) - (ELEMTYPE)rect_min_ref(rid, index)) *
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
  return RectSphericalVolume(rid);
  // return RectVolume(rid);
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
