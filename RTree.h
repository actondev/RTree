#include <vector>
#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#ifndef RTREE_H
#define RTREE_H

// NOTE This file compiles under MSVC 6 SP5 and MSVC .Net 2003 it may not work
// on other compilers without modification.

// NOTE These next few lines may be win32 specific, you may need to modify them
// to compile on other platform
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <functional>
#include <vector>

#define ASSERT assert // RTree uses ASSERT( condition )
#ifndef Min
#define Min std::min
#endif // Min
#ifndef Max
#define Max std::max
#endif // Max

//
// RTree.h
//

// datatype is what is stored in the tree (for example a gui element)
// elemtype is double or float
#define RTREE_TEMPLATE                                                         \
  template <class DATATYPE, class ELEMTYPE, class ELEMTYPEREAL,   \
            int TMAXNODES, int TMINNODES>
#define RTREE_QUAL                                                             \
  RTree<DATATYPE, ELEMTYPE, ELEMTYPEREAL, TMAXNODES, TMINNODES>

#define RTREE_USE_SPHERICAL_VOLUME // Better split classification, may be slower
                                   // on some systems

// Fwd decl
class RTFileStream; // File I/O helper class, look below for implementation and
                    // notes.

/// \class RTree
/// Implementation of RTree, a multidimensional bounding rectangle tree.
/// Example usage: For a 3-dimensional tree use RTree<Object*, float, 3> myTree;
///
/// This modified, templated C++ version by Greg Douglas at Auran
/// (http://www.auran.com)
///
/// DATATYPE Referenced data, should be int, void*, obj* etc. no larger than
/// sizeof<void*> and simple type ELEMTYPE Type of element such as int or float
/// NUMDIMS Number of dimensions such as 2 or 3
/// ELEMTYPEREAL Type of element that allows fractional and large values such as
/// float or double, for use in volume calcs
///
/// NOTES: Inserting and removing data requires the knowledge of its constant
/// Minimal Bounding Rectangle.
///        This version uses new/delete for nodes, I recommend using a fixed
///        size allocator for efficiency. Instead of using a callback function
///        for returned results, I recommend and efficient pre-sized, grow-only
///        memory array similar to MFC CArray or STL Vector for returning search
///        query result.
///
template <class DATATYPE, class ELEMTYPE = double,
          class ELEMTYPEREAL = ELEMTYPE, int TMAXNODES = 8,
          int TMINNODES = TMAXNODES / 2>
class RTree {
  static_assert(std::numeric_limits<ELEMTYPEREAL>::is_iec559,
                "'ELEMTYPEREAL' accepts floating-point types only");

  typedef std::function<bool(const DATATYPE &, const ELEMTYPE*, const ELEMTYPE*)> Callback;

protected:
  uint32_t dims; // set by the constructor
  struct Node;   // Fwd decl.  Used by other internal structs and iterator
  uint32_t count = 0;
public:
  // These constant must be declared after Branch and before Node struct
  // Stuck up here for MSVC 6 compiler.  NSVC .NET 2003 is much happier.
  enum {
    MAXNODES = TMAXNODES, ///< Max elements in node
    MINNODES = TMINNODES, ///< Min elements in node
  };

public:
  RTree() = delete;
  RTree(int dims);
  RTree(const RTree &other);
  virtual ~RTree();

  /// Insert entry
  /// \param a_min Min of bounding rect
  /// \param a_max Max of bounding rect
  /// \param a_dataId Positive Id of data.  Maybe zero, but negative numbers not
  /// allowed.
  void Insert(const ELEMTYPE *a_min, const ELEMTYPE *a_max,
              const DATATYPE &a_dataId);

  /// Remove entry \param a_min Min of bounding rect \param a_max Max
  /// of bounding rect \param predicate passing the stored type:
  /// should return true if node should be removed \return Returns the
  /// number of entries removed
  /// TODO this does not update the bounding rects (will later insertions be fine?)
  int Remove(const ELEMTYPE *a_min, const ELEMTYPE *a_max, Callback predicate);

  /// Find all within search rectangle
  /// \param a_min Min of search bounding rect
  /// \param a_max Max of search bounding rect
  /// \param a_searchResult Search result array. Caller should set grow size.
  /// Function will reset, not append to array. \param callback Callback
  /// function to return result. Callback should return 'true' to continue
  /// searching \return Returns the number of entries found
  int Search(const ELEMTYPE *a_min, const ELEMTYPE *a_max,
             Callback callback) const;

  /// Remove all entries from tree
  void RemoveAll();

  /// Move all nodes by an offset
  void MoveAll(const ELEMTYPE *a_offset);

  /// Count the data elements in this container.  This is slow as no internal
  /// counter is maintained.
  int Count();

  /// Load tree contents from file
  bool Load(const char *a_fileName);
  /// Load tree contents from stream
  bool Load(RTFileStream &a_stream);

  /// Save tree contents to file
  bool Save(const char *a_fileName);
  /// Save tree contents to stream
  bool Save(RTFileStream &a_stream);

protected:
  /// Minimal bounding rectangle (n-dimensional)
  struct Rect {
    std::vector<ELEMTYPE> m_min; ///< Min dimensions of bounding box
    std::vector<ELEMTYPE> m_max; ///< Max dimensions of bounding box
    Rect() = delete;
    Rect(const RTREE_QUAL* tree) {
      // cout << "rect ctor " <<tree->Dimensions() << endl;
      m_min.resize(tree->Dimensions());
      m_max.resize(tree->Dimensions());
    }
    Rect(int dims) {
      // cout << "rect dims ctor " << dims << endl;
      m_min.resize(dims);
      m_max.resize(dims);
    }
    // Rect(Rect&& other) = default;
  };

  /// May be data or may be another subtree
  /// The parents level determines this.
  /// If the parents level is 0, then this is data
  struct Branch {
    Rect m_rect;     ///< Bounds
    Node *m_child;   ///< Child node
    DATATYPE m_data; ///< Data Id
    Branch() = delete;
    Branch(const RTREE_QUAL* tree)
        : m_rect(tree){

    };
    Branch(int dims)
        : m_rect(dims){

          };
    Branch(const Branch& other)
        : m_rect(other.m_rect.m_min.size()){
      // cout << "Branch copy ctor" << endl;
    };

    // TODO learn a thing or two about move semantics !!
    // Branch(Branch&& other)
    //     : m_rect(other.m_rect.size())
    //     // m_rect(std::move(other.m_rect))
    // {
    //   cout << "Branch move ctor" << endl;
    // };

  };

  /// Node for each branch level
  struct Node {
    bool IsInternalNode() {
      return (m_level > 0);
    }                                        // Not a leaf, but a internal node
    bool IsLeaf() { return (m_level == 0); } // A leaf, contains data

    int m_count;               ///< Count
    int m_level;               ///< Leaf is zero, others positive
    std::vector<Branch> m_branch;
    Node() = delete;
    Node(const RTREE_QUAL* tree) {
      // cout << "Node ctor" <<endl;
      m_branch.resize(MAXNODES, Branch(tree));
    };
  };

  /// A link list of nodes for reinsertion after a delete operation
  struct ListNode {
    ListNode *m_next; ///< Next in list
    Node *m_node;     ///< Node
  };

  /// Variables for finding a split partition
  struct PartitionVars {
    enum { NOT_TAKEN = -1 }; // indicates that position

    int m_partition[MAXNODES + 1];
    int m_total;
    int m_minFill;
    int m_count[2];
    Rect m_cover[2];
    ELEMTYPEREAL m_area[2];
    std::vector<Branch> m_branchBuf;
    int m_branchCount;
    Rect m_coverSplit;
    ELEMTYPEREAL m_coverSplitArea;
    PartitionVars() = delete;
    PartitionVars(const RTREE_QUAL* tree)
        : m_cover{{tree}, {tree}},
          m_coverSplit{tree}
    {
      m_branchBuf.resize(MAXNODES +1, Branch{tree});
    };
  };

  Node *AllocNode();
  void FreeNode(Node *a_node);
  void InitNode(Node *a_node);
  void InitRect(Rect *a_rect);
  bool InsertRectRec(const Branch &a_branch, Node *a_node, Node **a_newNode,
                     int a_level);
  bool InsertRect(const Branch &a_branch, Node **a_root, int a_level);
  Rect NodeCover(Node *a_node);
  bool AddBranch(const Branch *a_branch, Node *a_node, Node **a_newNode);
  void DisconnectBranch(Node *a_node, int a_index);
  int PickBranch(const Rect *a_rect, Node *a_node);
  Rect CombineRect(const Rect *a_rectA, const Rect *a_rectB);
  void SplitNode(Node *a_node, const Branch *a_branch, Node **a_newNode);
  ELEMTYPEREAL RectSphericalVolume(Rect *a_rect);
  ELEMTYPEREAL RectVolume(Rect *a_rect);
  ELEMTYPEREAL CalcRectVolume(Rect *a_rect);
  void GetBranches(Node *a_node, const Branch *a_branch,
                   PartitionVars *a_parVars);
  void ChoosePartition(PartitionVars *a_parVars, int a_minFill);
  void LoadNodes(Node *a_nodeA, Node *a_nodeB, PartitionVars *a_parVars);
  void InitParVars(PartitionVars *a_parVars, int a_maxRects, int a_minFill);
  void PickSeeds(PartitionVars *a_parVars);
  void Classify(int a_index, int a_group, PartitionVars *a_parVars);
  bool RemoveRect(Node **a_root, Rect *a_rect, int &a_removedCount, Callback predicate);
  bool RemoveRectRec(Node *a_node, Rect *a_rect, int &a_removedCount, Callback predicate,
                     ListNode **a_listNode);
  ListNode *AllocListNode();
  void FreeListNode(ListNode *a_listNode);
  bool Overlap(Rect *a_rectA, Rect *a_rectB) const;
  void ReInsert(Node *a_node, ListNode **a_listNode);
  bool Search(Node *a_node, Rect *a_rect, int &a_foundCount,
              Callback callback) const;
  void RemoveAllRec(Node *a_node);
  void Reset();
  void MoveChildren(Node* a_node, const ELEMTYPE *a_offset);

  void CopyRec(Node *current, Node *other);

  Node *m_root;                    ///< Root of tree
  ELEMTYPEREAL m_unitSphereVolume; ///< Unit sphere constant for required number
                                   ///< of dimensions

public:
  // return all the AABBs that form the RTree
  std::vector<Rect> ListTree() const;

  Rect Bounds() const;
  int Dimensions() const;
};


RTREE_TEMPLATE
RTREE_QUAL::RTree(int dims) {
  ASSERT(MAXNODES > MINNODES);
  ASSERT(MINNODES > 0);
  this->dims = dims;

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

  m_root = AllocNode();
  m_root->m_level = 0;
  m_unitSphereVolume = (ELEMTYPEREAL)UNIT_SPHERE_VOLUMES[dims];
}

// RTREE_TEMPLATE
// RTREE_QUAL::RTree(const RTree &other) : RTree() {
//   CopyRec(m_root, other.m_root);
// }

RTREE_TEMPLATE
RTREE_QUAL::~RTree() {
  Reset(); // Free, or reset node memory
}

RTREE_TEMPLATE
void RTREE_QUAL::Insert(const ELEMTYPE *a_min, const ELEMTYPE *a_max,
                        const DATATYPE &a_dataId) {
#ifdef _DEBUG
  for (int index = 0; index < dims; ++index) {
    ASSERT(a_min[index] <= a_max[index]);
  }
#endif //_DEBUG

  Branch branch(this);
  branch.m_data = a_dataId;
  branch.m_child = NULL;

  for (unsigned int axis = 0; axis < dims; ++axis) {
    branch.m_rect.m_min[axis] = a_min[axis];
    branch.m_rect.m_max[axis] = a_max[axis];
  }

  InsertRect(branch, &m_root, 0);
  count++;
}

RTREE_TEMPLATE
int RTREE_QUAL::Remove(const ELEMTYPE *a_min, const ELEMTYPE *a_max,
                        Callback predicate) {
#ifdef _DEBUG
  for (int index = 0; index < dims; ++index) {
    ASSERT(a_min[index] <= a_max[index]);
  }
#endif //_DEBUG

  Rect rect(this);

  for (unsigned int axis = 0; axis < dims; ++axis) {
    rect.m_min[axis] = a_min[axis];
    rect.m_max[axis] = a_max[axis];
  }

  int removedCount = 0;
  RemoveRect(&m_root, &rect, removedCount, predicate);
  count -= removedCount;
  return removedCount;
}

RTREE_TEMPLATE
int RTREE_QUAL::Search(const ELEMTYPE *a_min, const ELEMTYPE *a_max,
                       Callback callback) const {
#ifdef _DEBUG
  for (int index = 0; index < dims; ++index) {
    ASSERT(a_min[index] <= a_max[index]);
  }
#endif //_DEBUG

  Rect rect(this);

  for (unsigned int axis = 0; axis < dims; ++axis) {
    rect.m_min[axis] = a_min[axis];
    rect.m_max[axis] = a_max[axis];
  }

  // NOTE: May want to return search result another way, perhaps returning the
  // number of found elements here.

  int foundCount = 0;
  Search(m_root, &rect, foundCount, callback);

  return foundCount;
}

RTREE_TEMPLATE
int RTREE_QUAL::Count() {
  return count;
}


RTREE_TEMPLATE
void RTREE_QUAL::MoveChildren(Node* a_node, const ELEMTYPE *a_offset) {
  ASSERT(a_node);
  ASSERT(a_node->m_level >= 0);
  ASSERT(a_offset);

  for (int index = 0; index < a_node->m_count; ++index) {
    Branch &branch = a_node->m_branch[index];
    Rect &rect = branch.m_rect;
    for (unsigned int i = 0; i < dims; i++) {
      rect.m_min[i] += a_offset[i];
      rect.m_max[i] += a_offset[i];
    }
    if (a_node->IsInternalNode()) {
      MoveChildren(branch.m_child, a_offset);
    }
  }
}

RTREE_TEMPLATE
void RTREE_QUAL::RemoveAll() {
  // Delete all existing nodes
  Reset();
  count = 0;

  m_root = AllocNode();
  m_root->m_level = 0;
}

RTREE_TEMPLATE
void RTREE_QUAL::MoveAll(const ELEMTYPE *a_offset) {
  MoveChildren(m_root, a_offset);
}

RTREE_TEMPLATE
void RTREE_QUAL::Reset() {
  RemoveAllRec(m_root);
}

RTREE_TEMPLATE
void RTREE_QUAL::RemoveAllRec(Node *a_node) {
  ASSERT(a_node);
  ASSERT(a_node->m_level >= 0);

  if (a_node->IsInternalNode()) // This is an internal node in the tree
  {
    for (int index = 0; index < a_node->m_count; ++index) {
      RemoveAllRec(a_node->m_branch[index].m_child);
    }
  }
  FreeNode(a_node);
}

RTREE_TEMPLATE
typename RTREE_QUAL::Node *RTREE_QUAL::AllocNode() {
  Node *newNode;
  newNode = new Node(this);
  InitNode(newNode);
  return newNode;
}

RTREE_TEMPLATE
void RTREE_QUAL::FreeNode(Node *a_node) {
  ASSERT(a_node);

  delete a_node;
}

// Allocate space for a node in the list used in DeletRect to
// store Nodes that are too empty.
RTREE_TEMPLATE
typename RTREE_QUAL::ListNode *RTREE_QUAL::AllocListNode() {
  return new ListNode;
}

RTREE_TEMPLATE
void RTREE_QUAL::FreeListNode(ListNode *a_listNode) {
  delete a_listNode;
}

RTREE_TEMPLATE
void RTREE_QUAL::InitNode(Node *a_node) {
  a_node->m_count = 0;
  a_node->m_level = -1;
}

RTREE_TEMPLATE
void RTREE_QUAL::InitRect(Rect *a_rect) {
  for (int index = 0; index < dims; ++index) {
    a_rect->m_min[index] = (ELEMTYPE)0;
    a_rect->m_max[index] = (ELEMTYPE)0;
  }
}

// Inserts a new data rectangle into the index structure.
// Recursively descends tree, propagates splits back up.
// Returns 0 if node was not split.  Old node updated.
// If node was split, returns 1 and sets the pointer pointed to by
// new_node to point to the new node.  Old node updated to become one of two.
// The level argument specifies the number of steps up from the leaf
// level to insert; e.g. a data rectangle goes in at level = 0.
RTREE_TEMPLATE
bool RTREE_QUAL::InsertRectRec(const Branch &a_branch, Node *a_node,
                               Node **a_newNode, int a_level) {
  ASSERT(a_node && a_newNode);
  ASSERT(a_level >= 0 && a_level <= a_node->m_level);

  // recurse until we reach the correct level for the new record. data records
  // will always be called with a_level == 0 (leaf)
  if (a_node->m_level > a_level) {
    // Still above level for insertion, go down tree recursively
    Node *otherNode;

    // find the optimal branch for this record
    int index = PickBranch(&a_branch.m_rect, a_node);

    // recursively insert this record into the picked branch
    bool childWasSplit = InsertRectRec(
        a_branch, a_node->m_branch[index].m_child, &otherNode, a_level);

    if (!childWasSplit) {
      // Child was not split. Merge the bounding box of the new record with the
      // existing bounding box
      a_node->m_branch[index].m_rect =
          CombineRect(&a_branch.m_rect, &(a_node->m_branch[index].m_rect));
      return false;
    } else {
      // Child was split. The old branches are now re-partitioned to two nodes
      // so we have to re-calculate the bounding boxes of each node
      a_node->m_branch[index].m_rect =
          NodeCover(a_node->m_branch[index].m_child);
      Branch branch(this);
      branch.m_child = otherNode;
      branch.m_rect = NodeCover(otherNode);

      // The old node is already a child of a_node. Now add the newly-created
      // node to a_node as well. a_node might be split because of that.
      return AddBranch(&branch, a_node, a_newNode);
    }
  } else if (a_node->m_level == a_level) {
    // We have reached level for insertion. Add rect, split if necessary
    return AddBranch(&a_branch, a_node, a_newNode);
  } else {
    // Should never occur
    ASSERT(0);
    return false;
  }
}

// Insert a data rectangle into an index structure.
// InsertRect provides for splitting the root;
// returns 1 if root was split, 0 if it was not.
// The level argument specifies the number of steps up from the leaf
// level to insert; e.g. a data rectangle goes in at level = 0.
// InsertRect2 does the recursion.
//
RTREE_TEMPLATE
bool RTREE_QUAL::InsertRect(const Branch &a_branch, Node **a_root,
                            int a_level) {
  ASSERT(a_root);
  ASSERT(a_level >= 0 && a_level <= (*a_root)->m_level);
#ifdef _DEBUG
  for (int index = 0; index < dims; ++index) {
    ASSERT(a_branch.m_rect.m_min[index] <= a_branch.m_rect.m_max[index]);
  }
#endif //_DEBUG

  Node *newNode;

  if (InsertRectRec(a_branch, *a_root, &newNode, a_level)) // Root split
  {
    // Grow tree taller and new root
    Node *newRoot = AllocNode();
    newRoot->m_level = (*a_root)->m_level + 1;

    Branch branch(this);

    // add old root node as a child of the new root
    branch.m_rect = NodeCover(*a_root);
    branch.m_child = *a_root;
    AddBranch(&branch, newRoot, NULL);

    // add the split node as a child of the new root
    branch.m_rect = NodeCover(newNode);
    branch.m_child = newNode;
    AddBranch(&branch, newRoot, NULL);

    // set the new root as the root node
    *a_root = newRoot;

    return true;
  }

  return false;
}

// Find the smallest rectangle that includes all rectangles in branches of a
// node.
RTREE_TEMPLATE
typename RTREE_QUAL::Rect RTREE_QUAL::NodeCover(Node *a_node) {
  ASSERT(a_node);

  Rect rect = a_node->m_branch[0].m_rect;
  for (int index = 1; index < a_node->m_count; ++index) {
    rect = CombineRect(&rect, &(a_node->m_branch[index].m_rect));
  }

  return rect;
}

// Add a branch to a node.  Split the node if necessary.
// Returns 0 if node not split.  Old node updated.
// Returns 1 if node split, sets *new_node to address of new node.
// Old node updated, becomes one of two.
RTREE_TEMPLATE
bool RTREE_QUAL::AddBranch(const Branch *a_branch, Node *a_node,
                           Node **a_newNode) {
  ASSERT(a_branch);
  ASSERT(a_node);

  if (a_node->m_count < MAXNODES) // Split won't be necessary
  {
    a_node->m_branch[a_node->m_count] = *a_branch;
    ++a_node->m_count;

    return false;
  } else {
    ASSERT(a_newNode);

    SplitNode(a_node, a_branch, a_newNode);
    return true;
  }
}

// Disconnect a dependent node.
// Caller must return (or stop using iteration index) after this as count has
// changed
RTREE_TEMPLATE
void RTREE_QUAL::DisconnectBranch(Node *a_node, int a_index) {
  ASSERT(a_node && (a_index >= 0) && (a_index < MAXNODES));
  ASSERT(a_node->m_count > 0);

  // Remove element by swapping with the last element to prevent gaps in array
  a_node->m_branch[a_index] = a_node->m_branch[a_node->m_count - 1];

  --a_node->m_count;
}

// Pick a branch.  Pick the one that will need the smallest increase
// in area to accomodate the new rectangle.  This will result in the
// least total area for the covering rectangles in the current node.
// In case of a tie, pick the one which was smaller before, to get
// the best resolution when searching.
RTREE_TEMPLATE
int RTREE_QUAL::PickBranch(const Rect *a_rect, Node *a_node) {
  ASSERT(a_rect && a_node);

  bool firstTime = true;
  ELEMTYPEREAL increase;
  ELEMTYPEREAL bestIncr = (ELEMTYPEREAL)-1;
  ELEMTYPEREAL area;
  ELEMTYPEREAL bestArea;
  int best = 0;
  Rect tempRect(this);

  for (int index = 0; index < a_node->m_count; ++index) {
    Rect *curRect = &a_node->m_branch[index].m_rect;
    area = CalcRectVolume(curRect);
    tempRect = CombineRect(a_rect, curRect);
    increase = CalcRectVolume(&tempRect) - area;
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

// Combine two rectangles into larger one containing both
RTREE_TEMPLATE
typename RTREE_QUAL::Rect RTREE_QUAL::CombineRect(const Rect *a_rectA,
                                                  const Rect *a_rectB) {
  ASSERT(a_rectA && a_rectB);

  Rect newRect(this);

  for (unsigned int index = 0; index < dims; ++index) {
    newRect.m_min[index] = Min(a_rectA->m_min[index], a_rectB->m_min[index]);
    newRect.m_max[index] = Max(a_rectA->m_max[index], a_rectB->m_max[index]);
  }

  return newRect;
}

// Split a node.
// Divides the nodes branches and the extra one between two nodes.
// Old node is one of the new ones, and one really new one is created.
// Tries more than one method for choosing a partition, uses best result.
RTREE_TEMPLATE
void RTREE_QUAL::SplitNode(Node *a_node, const Branch *a_branch,
                           Node **a_newNode) {
  ASSERT(a_node);
  ASSERT(a_branch);

  // Could just use local here, but member or external is faster since it is
  // reused
  PartitionVars localVars(this);
  PartitionVars *parVars = &localVars;

  // Load all the branches into a buffer, initialize old node
  GetBranches(a_node, a_branch, parVars);

  // Find partition
  ChoosePartition(parVars, MINNODES);

  // Create a new node to hold (about) half of the branches
  *a_newNode = AllocNode();
  (*a_newNode)->m_level = a_node->m_level;

  // Put branches from buffer into 2 nodes according to the chosen partition
  a_node->m_count = 0;
  LoadNodes(a_node, *a_newNode, parVars);

  ASSERT((a_node->m_count + (*a_newNode)->m_count) == parVars->m_total);
}

// Calculate the n-dimensional volume of a rectangle
RTREE_TEMPLATE
ELEMTYPEREAL RTREE_QUAL::RectVolume(Rect *a_rect) {
  ASSERT(a_rect);

  ELEMTYPEREAL volume = (ELEMTYPEREAL)1;

  for (int index = 0; index < dims; ++index) {
    volume *= a_rect->m_max[index] - a_rect->m_min[index];
  }

  ASSERT(volume >= (ELEMTYPEREAL)0);

  return volume;
}

// The exact volume of the bounding sphere for the given Rect
RTREE_TEMPLATE
ELEMTYPEREAL RTREE_QUAL::RectSphericalVolume(Rect *a_rect) {
  ASSERT(a_rect);

  ELEMTYPEREAL sumOfSquares = (ELEMTYPEREAL)0;
  ELEMTYPEREAL radius;

  for (unsigned int index = 0; index < dims; ++index) {
    ELEMTYPEREAL halfExtent = ((ELEMTYPEREAL)a_rect->m_max[index] -
                               (ELEMTYPEREAL)a_rect->m_min[index]) *
                              (ELEMTYPEREAL)0.5;
    sumOfSquares += halfExtent * halfExtent;
  }

  radius = (ELEMTYPEREAL)sqrt(sumOfSquares);

  // Pow maybe slow, so test for common dims like 2,3 and just use x*x, x*x*x.
  if (dims == 3) {
    return (radius * radius * radius * m_unitSphereVolume);
  } else if (dims == 2) {
    return (radius * radius * m_unitSphereVolume);
  } else {
    return (ELEMTYPEREAL)(pow(radius, dims) * m_unitSphereVolume);
  }
}

// Use one of the methods to calculate retangle volume
RTREE_TEMPLATE
ELEMTYPEREAL RTREE_QUAL::CalcRectVolume(Rect *a_rect) {
#ifdef RTREE_USE_SPHERICAL_VOLUME
  return RectSphericalVolume(a_rect); // Slower but helps certain merge cases
#else                                 // RTREE_USE_SPHERICAL_VOLUME
  return RectVolume(a_rect); // Faster but can cause poor merges
#endif                                // RTREE_USE_SPHERICAL_VOLUME
}

// Load branch buffer with branches from full node plus the extra branch.
RTREE_TEMPLATE
void RTREE_QUAL::GetBranches(Node *a_node, const Branch *a_branch,
                             PartitionVars *a_parVars) {
  ASSERT(a_node);
  ASSERT(a_branch);

  ASSERT(a_node->m_count == MAXNODES);

  // Load the branch buffer
  for (int index = 0; index < MAXNODES; ++index) {
    a_parVars->m_branchBuf[index] = a_node->m_branch[index];
  }
  a_parVars->m_branchBuf[MAXNODES] = *a_branch;
  a_parVars->m_branchCount = MAXNODES + 1;

  // Calculate rect containing all in the set
  a_parVars->m_coverSplit = a_parVars->m_branchBuf[0].m_rect;
  for (int index = 1; index < MAXNODES + 1; ++index) {
    a_parVars->m_coverSplit = CombineRect(
        &a_parVars->m_coverSplit, &a_parVars->m_branchBuf[index].m_rect);
  }
  a_parVars->m_coverSplitArea = CalcRectVolume(&a_parVars->m_coverSplit);
}

// Method #0 for choosing a partition:
// As the seeds for the two groups, pick the two rects that would waste the
// most area if covered by a single rectangle, i.e. evidently the worst pair
// to have in the same group.
// Of the remaining, one at a time is chosen to be put in one of the two groups.
// The one chosen is the one with the greatest difference in area expansion
// depending on which group - the rect most strongly attracted to one group
// and repelled from the other.
// If one group gets too full (more would force other group to violate min
// fill requirement) then other group gets the rest.
// These last are the ones that can go in either group most easily.
RTREE_TEMPLATE
void RTREE_QUAL::ChoosePartition(PartitionVars *a_parVars, int a_minFill) {
  ASSERT(a_parVars);

  ELEMTYPEREAL biggestDiff;
  int group, chosen = 0, betterGroup = 0;

  InitParVars(a_parVars, a_parVars->m_branchCount, a_minFill);
  PickSeeds(a_parVars);

  while (
      ((a_parVars->m_count[0] + a_parVars->m_count[1]) < a_parVars->m_total) &&
      (a_parVars->m_count[0] < (a_parVars->m_total - a_parVars->m_minFill)) &&
      (a_parVars->m_count[1] < (a_parVars->m_total - a_parVars->m_minFill))) {
    biggestDiff = (ELEMTYPEREAL)-1;
    for (int index = 0; index < a_parVars->m_total; ++index) {
      if (PartitionVars::NOT_TAKEN == a_parVars->m_partition[index]) {
        Rect *curRect = &a_parVars->m_branchBuf[index].m_rect;
        Rect rect0 = CombineRect(curRect, &a_parVars->m_cover[0]);
        Rect rect1 = CombineRect(curRect, &a_parVars->m_cover[1]);
        ELEMTYPEREAL growth0 = CalcRectVolume(&rect0) - a_parVars->m_area[0];
        ELEMTYPEREAL growth1 = CalcRectVolume(&rect1) - a_parVars->m_area[1];
        ELEMTYPEREAL diff = growth1 - growth0;
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
        } else if ((diff == biggestDiff) && (a_parVars->m_count[group] <
                                             a_parVars->m_count[betterGroup])) {
          chosen = index;
          betterGroup = group;
        }
      }
    }
    Classify(chosen, betterGroup, a_parVars);
  }

  // If one group too full, put remaining rects in the other
  if ((a_parVars->m_count[0] + a_parVars->m_count[1]) < a_parVars->m_total) {
    if (a_parVars->m_count[0] >= a_parVars->m_total - a_parVars->m_minFill) {
      group = 1;
    } else {
      group = 0;
    }
    for (int index = 0; index < a_parVars->m_total; ++index) {
      if (PartitionVars::NOT_TAKEN == a_parVars->m_partition[index]) {
        Classify(index, group, a_parVars);
      }
    }
  }

  ASSERT((a_parVars->m_count[0] + a_parVars->m_count[1]) == a_parVars->m_total);
  ASSERT((a_parVars->m_count[0] >= a_parVars->m_minFill) &&
         (a_parVars->m_count[1] >= a_parVars->m_minFill));
}

// Copy branches from the buffer into two nodes according to the partition.
RTREE_TEMPLATE
void RTREE_QUAL::LoadNodes(Node *a_nodeA, Node *a_nodeB,
                           PartitionVars *a_parVars) {
  ASSERT(a_nodeA);
  ASSERT(a_nodeB);
  ASSERT(a_parVars);

  for (int index = 0; index < a_parVars->m_total; ++index) {
    ASSERT(a_parVars->m_partition[index] == 0 ||
           a_parVars->m_partition[index] == 1);

    int targetNodeIndex = a_parVars->m_partition[index];
    Node *targetNodes[] = {a_nodeA, a_nodeB};

    // It is assured that AddBranch here will not cause a node split.
    bool nodeWasSplit = AddBranch(&a_parVars->m_branchBuf[index],
                                  targetNodes[targetNodeIndex], NULL);
    ASSERT(!nodeWasSplit);
  }
}

// Initialize a PartitionVars structure.
RTREE_TEMPLATE
void RTREE_QUAL::InitParVars(PartitionVars *a_parVars, int a_maxRects,
                             int a_minFill) {
  ASSERT(a_parVars);

  a_parVars->m_count[0] = a_parVars->m_count[1] = 0;
  a_parVars->m_area[0] = a_parVars->m_area[1] = (ELEMTYPEREAL)0;
  a_parVars->m_total = a_maxRects;
  a_parVars->m_minFill = a_minFill;
  for (int index = 0; index < a_maxRects; ++index) {
    a_parVars->m_partition[index] = PartitionVars::NOT_TAKEN;
  }
}

RTREE_TEMPLATE
void RTREE_QUAL::PickSeeds(PartitionVars *a_parVars) {
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

  Classify(seed0, 0, a_parVars);
  Classify(seed1, 1, a_parVars);
}

// Put a branch in one of the groups.
RTREE_TEMPLATE
void RTREE_QUAL::Classify(int a_index, int a_group, PartitionVars *a_parVars) {
  ASSERT(a_parVars);
  ASSERT(PartitionVars::NOT_TAKEN == a_parVars->m_partition[a_index]);

  a_parVars->m_partition[a_index] = a_group;

  // Calculate combined rect
  if (a_parVars->m_count[a_group] == 0) {
    a_parVars->m_cover[a_group] = a_parVars->m_branchBuf[a_index].m_rect;
  } else {
    a_parVars->m_cover[a_group] = CombineRect(
        &a_parVars->m_branchBuf[a_index].m_rect, &a_parVars->m_cover[a_group]);
  }

  // Calculate volume of combined rect
  a_parVars->m_area[a_group] = CalcRectVolume(&a_parVars->m_cover[a_group]);

  ++a_parVars->m_count[a_group];
}

// Delete a data rectangle from an index structure.
// Pass in a pointer to a Rect, the tid of the record, ptr to ptr to root node.
// Returns true if something removed
// RemoveRect provides for eliminating the root.
RTREE_TEMPLATE
bool RTREE_QUAL::RemoveRect(Node **a_root, Rect *a_rect, int &a_foundCount, Callback predicate ) {
  ASSERT(a_rect && a_root);
  ASSERT(*a_root);

  ListNode *reInsertList = NULL;

  if (!RemoveRectRec(*a_root, a_rect, a_foundCount, predicate , &reInsertList)) {
    return false; // nothing removed
  }

  // removed
  while (reInsertList) {
    Node *tempNode = reInsertList->m_node;

    for (int index = 0; index < tempNode->m_count; ++index) {
      // TODO go over this code. should I use (tempNode->m_level - 1)?
      InsertRect(tempNode->m_branch[index], a_root, tempNode->m_level);
    }

    ListNode *remLNode = reInsertList;
    reInsertList = reInsertList->m_next;

    FreeNode(remLNode->m_node);
    FreeListNode(remLNode);
  }

  // Check for redundant root (not leaf, 1 child) and eliminate TODO replace
  // if with while? In case there is a whole branch of redundant roots...
  if ((*a_root)->m_count == 1 && (*a_root)->IsInternalNode()) {
    Node *tempNode = (*a_root)->m_branch[0].m_child;

    ASSERT(tempNode);
    FreeNode(*a_root);
    *a_root = tempNode;
  }
  return true;
}

// Delete a rectangle from non-root part of an index structure.
// Called by RemoveRect.  Descends tree recursively,
// merges branches on the way back up.
// Returns true if something removed
RTREE_TEMPLATE
bool RTREE_QUAL::RemoveRectRec(Node *a_node, Rect *a_rect, int &a_removedCount, Callback predicate, ListNode **a_listNode) {
  ASSERT(a_rect && a_node && a_listNode);
  ASSERT(a_node->m_level >= 0);

  if (a_node->IsInternalNode()) {
    // not a leaf node
    bool removed = false;
    for (int index = 0; index < a_node->m_count; ++index) {
      if (Overlap(a_rect, &(a_node->m_branch[index].m_rect))) {
        if (RemoveRectRec(a_node->m_branch[index].m_child, a_rect, a_removedCount,
                          predicate,
                           a_listNode)) {
          if (a_node->m_branch[index].m_child->m_count >= MINNODES) {
            // child removed, just resize parent rect
            a_node->m_branch[index].m_rect =
                NodeCover(a_node->m_branch[index].m_child);
          } else {
            // child removed, not enough entries in node, eliminate node
            ReInsert(a_node->m_branch[index].m_child, a_listNode);
            DisconnectBranch(a_node, index);
            // NB: Before remove refactor this was returning
            index--; // have to revisit same index, as now it's swapped with the last item
          }
          removed = true;
        }
      }
    }
    return removed;
  } else {
    // A leaf node
    bool removed = false;
    for (int index = 0; index < a_node->m_count; ++index) {
      if (Overlap(a_rect, &a_node->m_branch[index].m_rect)) {
        Branch& branch = a_node->m_branch[index];
        if (predicate(branch.m_data, branch.m_rect.m_min.data(), branch.m_rect.m_max.data() )) {
          removed = true;
          DisconnectBranch(a_node, index);
          // NB: Before remove refactor this was returning
          index--; // have to revisit same index, as now it's swapped with the last item
          a_removedCount++;
        }
      }
    }
    return removed;
  }
}

// Decide whether two rectangles overlap.
RTREE_TEMPLATE
bool RTREE_QUAL::Overlap(Rect *a_rectA, Rect *a_rectB) const {
  ASSERT(a_rectA && a_rectB);

  for (unsigned int index = 0; index < dims; ++index) {
    if (a_rectA->m_min[index] > a_rectB->m_max[index] ||
        a_rectB->m_min[index] > a_rectA->m_max[index]) {
      return false;
    }
  }
  return true;
}

// Add a node to the reinsertion list.  All its branches will later
// be reinserted into the index structure.
RTREE_TEMPLATE
void RTREE_QUAL::ReInsert(Node *a_node, ListNode **a_listNode) {
  ListNode *newListNode;

  newListNode = AllocListNode();
  newListNode->m_node = a_node;
  newListNode->m_next = *a_listNode;
  *a_listNode = newListNode;
}

// Search in an index tree or subtree for all data retangles that overlap the
// argument rectangle.
RTREE_TEMPLATE
bool RTREE_QUAL::Search(Node *a_node, Rect *a_rect, int &a_foundCount,
                        Callback callback) const {
  ASSERT(a_node);
  ASSERT(a_node->m_level >= 0);
  ASSERT(a_rect);

  if (a_node->IsInternalNode()) {
    // This is an internal node in the tree
    for (int index = 0; index < a_node->m_count; ++index) {
      if (Overlap(a_rect, &a_node->m_branch[index].m_rect)) {
        if (!Search(a_node->m_branch[index].m_child, a_rect, a_foundCount,
                    callback)) {
          // The callback indicated to stop searching
          return false;
        }
      }
    }
  } else {
    // This is a leaf node
    for (int index = 0; index < a_node->m_count; ++index) {
      if (Overlap(a_rect, &a_node->m_branch[index].m_rect)) {
        Branch& branch = a_node->m_branch[index];
        DATATYPE &data = branch.m_data;
        ++a_foundCount;

        if (!callback(data, branch.m_rect.m_min.data(), branch.m_rect.m_max.data())) {
          return false; // Don't continue searching
        }
      }
    }
  }

  return true; // Continue searching
}

RTREE_TEMPLATE
std::vector<typename RTREE_QUAL::Rect> RTREE_QUAL::ListTree() const {
  ASSERT(m_root);
  ASSERT(m_root->m_level >= 0);

  std::vector<Rect> treeList;

  std::vector<Node *> toVisit;
  toVisit.push_back(m_root);

  while (!toVisit.empty()) {
    Node *a_node = toVisit.back();
    toVisit.pop_back();
    if (a_node->IsInternalNode()) {
      // This is an internal node in the tree
      for (int index = 0; index < a_node->m_count; ++index) {
        treeList.push_back(a_node->m_branch[index].m_rect);
        toVisit.push_back(a_node->m_branch[index].m_child);
      }
    } else {
      // This is a leaf node
      for (int index = 0; index < a_node->m_count; ++index) {
        treeList.push_back(a_node->m_branch[index].m_rect);
      }
    }
  }

  return treeList;
}

RTREE_TEMPLATE
typename RTREE_QUAL::Rect RTREE_QUAL::Bounds() const {
  ASSERT(m_root);
  ASSERT(m_root->m_level >= 0);
  if(m_root->m_count == 0) {
    Rect bounds(this);
    return bounds;
  }

  Rect bounds(this);
  Branch &first_branch = m_root->m_branch[0];
  bounds = first_branch.m_rect; // init
  for (int branch_id = 1; branch_id < m_root->m_count; branch_id++) {
    Branch &branch = m_root->m_branch[branch_id];
    Rect &other_rect = branch.m_rect;
    for (unsigned int index = 0; index < dims; index++) {
      bounds.m_min[index] = Min(bounds.m_min[index], other_rect.m_min[index]);
      bounds.m_max[index] = Max(bounds.m_max[index], other_rect.m_max[index]);
    }
  }

  return bounds;
}

RTREE_TEMPLATE
int RTREE_QUAL::Dimensions() const { return this->dims; }

// for the copy ctor
RTREE_TEMPLATE
void RTREE_QUAL::CopyRec(Node *current, Node *other) {
  current->m_level = other->m_level;
  current->m_count = other->m_count;

  if (current->IsInternalNode()) // not a leaf node
  {
    for (int index = 0; index < current->m_count; ++index) {
      Branch *currentBranch = &current->m_branch[index];
      Branch *otherBranch = &other->m_branch[index];

      currentBranch->m_rect.m_min = otherBranch->m_rect.m_min;
      // std::copy(otherBranch->m_rect.m_min, otherBranch->m_rect.m_min + dims,
      //           currentBranch->m_rect.m_min);

      currentBranch->m_rect.m_max = otherBranch->m_rect.m_max;
      // std::copy(otherBranch->m_rect.m_max, otherBranch->m_rect.m_max + dims,
      //           currentBranch->m_rect.m_max);

      currentBranch->m_child = AllocNode();
      CopyRec(currentBranch->m_child, otherBranch->m_child);
    }
  } else // A leaf node
  {
    for (int index = 0; index < current->m_count; ++index) {
      Branch *currentBranch = &current->m_branch[index];
      Branch *otherBranch = &other->m_branch[index];

      currentBranch->m_rect.m_min = otherBranch->m_rect.m_min;
      // std::copy(otherBranch->m_rect.m_min, otherBranch->m_rect.m_min + dims,
      //           currentBranch->m_rect.m_min);

      // std::copy(otherBranch->m_rect.m_max, otherBranch->m_rect.m_max + dims,
      //           currentBranch->m_rect.m_max);
      currentBranch->m_rect.m_max = otherBranch->m_rect.m_max;

      currentBranch->m_data = otherBranch->m_data;
    }
  }
}

#undef RTREE_TEMPLATE
#undef RTREE_QUAL

#endif // RTREE_H
