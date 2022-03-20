#pragma once

#include <cstdint>
#include <optional>
#include <vector>
#include <cstddef>
#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

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

#define ASSERT assert // drtree uses ASSERT( condition )
#ifndef Min
#define Min std::min
#endif // Min
#ifndef Max
#define Max std::max
#endif // Max

#define DRTREE_USE_SPHERICAL_VOLUME // Better split classification, may be slower
                                    // on some systems

// with RECT_ONE_ARRAY it seems slower (wrt RectSphericalVolume)
// (also with -O2)
// #define RECT_ONE_ARRAY

#ifdef RECT_ONE_ARRAY
#define RECT_MIN(rect) rect.array
#define RECT_MAX(rect) rect.array + rect.dims

#define RECTP_MIN(rect) rect->array
#define RECTP_MAX(rect) rect->array + rect->dims

#define RECT_MIN_REF(rect, i) rect.array[i]
#define RECT_MAX_REF(rect, i) rect.array[rect.dims + i]

#define RECTP_MIN_REF(rect, i) rect->array[i]
#define RECTP_MAX_REF(rect, i) rect->array[i + rect->dims]

#else

#define RECT_MIN(rect) rect.m_min
#define RECT_MAX(rect) rect.m_max

#define RECTP_MIN(rect) rect.m_min
#define RECTP_MAX(rect) rect.m_max

#define RECT_MIN_REF(rect, i) rect.m_min[i]
#define RECT_MAX_REF(rect, i) rect.m_max[i]

#define RECTP_MIN_REF(rect, i) rect->m_min[i]
#define RECTP_MAX_REF(rect, i) rect->m_max[i]

#endif


class FixedAllocator {
 private:
  size_t size;
  size_t allocated = 0;
  uint8_t* block = nullptr;
 public:
  FixedAllocator() = delete;
  FixedAllocator(size_t size) : size(size) {
    // cout << "allocator ctor, size " << size << endl;
    block = new uint8_t[size];
    if(!block) {
      cout << "could not allocate" << endl;
    }
  }
  ~FixedAllocator() {
    delete[] block;
  }
  void print_stats() {
    printf("Allocator stats: Init size %zu, allocated %zu\n",
           this->size, allocated);
  }
  void* allocate(size_t size) {
    uint8_t* res = block + allocated;
    allocated += size;
    if (allocated > this->size) {
      printf("Allocator exceeded size. Init size %zu allocated %zu requested "
             "%zu\n",
             this->size, allocated, size);
    }

    // printf("allocated %zu, totally %zu out of %zu\n", size, allocated, this->size);
    // cout << "allocating " << size << endl;
    return (void*)res;
    }
};

// datatype is what is stored in the tree (for example a gui element)
// elemtype is double or float
#define DRTREE_TEMPLATE                                                         \
  template <class DATATYPE, class ELEMTYPE, class ELEMTYPEREAL,   \
            int TMAXNODES, int TMINNODES>
#define DRTREE_QUAL                                                             \
  drtree<DATATYPE, ELEMTYPE, ELEMTYPEREAL, TMAXNODES, TMINNODES>

// Fwd decl
class RTFileStream; // File I/O helper class, look below for implementation and
                    // notes.

/// \class drtree
/// Implementation of drtree, a multidimensional bounding rectangle tree.
/// Example usage: For a 3-dimensional tree use drtree<Object*, float, 3> myTree;
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
class drtree {
  static_assert(std::numeric_limits<ELEMTYPEREAL>::is_iec559,
                "'ELEMTYPEREAL' accepts floating-point types only");

  typedef std::function<bool(const DATATYPE &, const ELEMTYPE*, const ELEMTYPE*)> Callback;

protected:
  uint32_t dims; // set by the constructor
  struct Node;   // Fwd decl.  Used by other internal structs and iterator
  struct Branch;
  struct Rect; // TODO rename to _Rect (not exposed to the outside, constructed only with allocator)
  struct PublicRect; // using std::vector, to pass info to the user (deconstruction is easy)
  struct PartitionVars;
  uint32_t count = 0;

 private:
  FixedAllocator m_temps_allocator;
  Branch m_insert_branch;
  Branch m_insert_rect_rec_branch;
  Branch m_insert_rect_branch;
  Rect m_remove_rect;
  Rect m_search_rect;
  Rect m_pick_branch_rect;
  Rect m_pick_seeds_rect;
  Rect m_choose_partition_rect0;
  Rect m_choose_partition_rect1;
  PartitionVars m_split_parition_vars;
public:
  // These constant must be declared after Branch and before Node struct
  // Stuck up here for MSVC 6 compiler.  NSVC .NET 2003 is much happier.
  
  enum {
    MAXNODES = TMAXNODES, ///< Max elements in node
    MINNODES = TMINNODES, ///< Min elements in node
  };

public:
  drtree() = delete;
  drtree(int dims);
  drtree(const drtree &other);
  virtual ~drtree();

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
             Callback callback);

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
  struct PublicRect {
    std::vector<double> low;
    std::vector<double> high;

    PublicRect() = delete;
    PublicRect(const DRTREE_QUAL* tree) : PublicRect(tree->Dimensions()) {
    }
    PublicRect(int dims) {
      low.resize(dims, 0);
      high.resize(dims, 0);
    }
    
  };
  
  /// Minimal bounding rectangle (n-dimensional)
  struct Rect {
    #ifdef RECT_ONE_ARRAY
    double* array;
    #else
    double* m_min; ///< Min dimensions of bounding box
    double* m_max; ///< Max dimensions of bounding box
    #endif
    const uint8_t dims;
    Rect() = delete;
    Rect(const DRTREE_QUAL* tree)
        : Rect(tree->Dimensions())
    {
      // cout << "rect ctor " <<tree->Dimensions() << endl;
    }
    static size_t heap_size(int dims) {
      return  2 * dims * sizeof(double);
    }

    Rect(int dims, FixedAllocator* allocator = nullptr)
        :dims(dims)
    {
      ASSERT(allocator);
#ifdef RECT_ONE_ARRAY
      array = (double*) allocator->allocate(2 * dims * sizeof(double));
#else
      m_min = (double*) allocator->allocate(dims * sizeof(double));
      m_max = (double*) allocator->allocate(dims * sizeof(double));
#endif
    }

    Rect(const Rect& other) = delete;

    Rect &operator=(const Rect &other) {
      // without this when  `a_parVars->m_coverSplit = a_parVars->m_branchBuf[0].m_rect;`
      // and modifying a_parVars->m_coverSplit, it was modifying a_parVars->m_branchBuf[0].m_rect

      // so.. the default copy assignment was copying the m_min & m_max pointers?
      // cout << "Rect copy assignment" << endl;
#ifdef RECT_ONE_ARRAY
      std::copy(other.array, other.array + 2*dims, array);
#else
      std::copy(other.m_min, other.m_min + dims, m_min);
      std::copy(other.m_max, other.m_max + dims, m_max);
#endif
      return *this;
    }
  };

  /// May be data or may be another subtree
  /// The parents level determines this.
  /// If the parents level is 0, then this is data
  struct Branch {
    Rect m_rect;     ///< Bounds
    Node *m_child;   ///< Child node
    // Note: made it optional initially cause of other issues,
    // but I suppose it's useful: what if the stored datatype does not have a default constructor?
    DATATYPE m_data; ///< Data Id
    int m_dims;
    Branch() = delete;
    Branch(const DRTREE_QUAL* tree)
        : Branch(tree->Dimensions()) {
    }
    Branch(int dims, FixedAllocator* allocator = nullptr)
        : m_rect(dims, allocator)  {
      m_dims = dims;
      // cout << "Branch dims ctor" << endl;
          };
    static size_t heap_size(int dims) {
      return Rect::heap_size(dims);
    }
    Branch(const Branch &other)
        : Branch(other.m_dims){
      // cout << "Branch copy ctor" << endl;
          };

    // TODO learn a thing or two about move semantics !!
    // disabled this cause was getting problems with copy assignment having bugs?
    // TODO move rect here
    // Branch(Branch&& other)
    //     : m_rect(m_dims) {
    //       cout << "Branch move ctor" << endl;
    // };

    // copy assignment operator
    // hm.. TODO copy ctor for Rect?
  //   Branch& operator=(const Branch &rhs) {
  //     // TODO fix this one. copy data??
  //     // cout << "branch copy assignment operator" << endl;
  //     m_rect.m_min = rhs.m_rect.m_min;
  //     m_rect.m_max = rhs.m_rect.m_max;
  //     m_rect.dims = rhs.m_rect.dims;
  //     m_child = rhs.m_child;
  //     if(rhs.m_data) {
  //       cout << "copy assignment, can't move data??" << endl;
  //       m_data = rhs.m_data;
  //     };
  //     return *this;
  //   }
  };

  /// Node for each branch level
  struct Node {
    bool IsInternalNode() {
      return (m_level > 0);
    }                                        // Not a leaf, but a internal node
    bool IsLeaf() { return (m_level == 0); } // A leaf, contains data

    int m_count; ///< Count
    int m_level; ///< Leaf is zero, others positive
    Branch* m_branch; // we need MAXNODES branches
    Node() = delete;
    int dims;
    FixedAllocator allocator;
    Node(const DRTREE_QUAL *tree)
        : Node(tree->Dimensions()) { };
    Node(int dims)
        : dims(dims),
          allocator(MAXNODES * (sizeof(Branch) + Branch::heap_size(dims)))
    {
      m_branch = (Branch*)allocator.allocate(MAXNODES * sizeof(Branch));
      for(int i=0; i<MAXNODES; i++) {
        // hmm move semantics could be used instead of placement new
        // but had problems with move & copy assignment definitions
        // Branch b{dims};
        // m_branch[i] = std::move(b);
        // new((void*)(m_branch+i)) Branch(dims);
        new(m_branch+i) Branch(dims, &allocator);
      }
      // cout << "Node ctor done" <<endl;
    }
    Node(const Node &other) : Node(other.dims) {
      cout << "Node copy ctor" << endl;
    };
    ~Node() {
      // delete m_branch;
      for (int i = 0; i < MAXNODES; i++) {
        m_branch[i].~Branch();
      }
    }
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

    // important to haves this before the rects
    // or generally anything that wants to use this
    // cause of the Constructor initialization-list evaluation order
    // https://stackoverflow.com/questions/1242830/constructor-initialization-list-evaluation-order
    FixedAllocator allocator;

    Rect m_coverSplit;
    Rect m_cover[2];
    ELEMTYPEREAL m_area[2];
    Branch* m_branchBuf;
    int m_branchCount;
    ELEMTYPEREAL m_coverSplitArea;
    PartitionVars() = delete;
    PartitionVars(const DRTREE_QUAL* tree) : PartitionVars(tree->Dimensions()) {}
    PartitionVars(int dims)
        : allocator((MAXNODES+1) * (sizeof(Branch) + Branch::heap_size(dims))
                     // plus the m_coverSplit
                    + 1 * Rect::heap_size(dims)
                    // plus m_cover (Rect[2])
                    + 2 * Rect::heap_size(dims)
                    ),
          m_coverSplit{dims, &allocator },
          m_cover{{dims, &allocator}, {dims, &allocator }}
    {
      m_branchBuf = (Branch *)allocator.allocate((MAXNODES+1) * sizeof(Branch));
      for (int i = 0; i < MAXNODES+1; i++) {
        new (m_branchBuf + i) Branch(dims, &allocator);
      }
      // allocator.print_stats();
    };
    ~PartitionVars() {
      for (int i = 0; i < MAXNODES + 1; i++) {
        m_branchBuf[i].~Branch();
      }
    }
  };

  Node *AllocNode();
  void FreeNode(Node *a_node);
  void InitNode(Node *a_node);
  bool InsertRectRec(const Branch &a_branch, Node *a_node, Node **a_newNode,
                     int a_level);
  bool InsertRect(const Branch &a_branch, Node **a_root, int a_level);
  void node_cover(Rect* dst, Node *a_node);
  bool AddBranch(const Branch *a_branch, Node *a_node, Node **a_newNode);
  void DisconnectBranch(Node *a_node, int a_index);
  int PickBranch(const Rect *a_rect, Node *a_node);
  void combine_rects(const Rect *dst, const Rect *a, const Rect *b);
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
              Callback callback);
  void RemoveAllRec(Node *a_node);
  void Reset();
  void MoveChildren(Node* a_node, const ELEMTYPE *a_offset);

  void CopyRec(Node *current, Node *other);

  Node *m_root;                    ///< Root of tree
  ELEMTYPEREAL m_unitSphereVolume; ///< Unit sphere constant for required number
                                   ///< of dimensions

public:
  // return all the AABBs that form the drtree
  std::vector<Rect> ListTree() const;

  PublicRect Bounds();
  int Dimensions() const;
};


DRTREE_TEMPLATE
DRTREE_QUAL::drtree(int dims)
    : m_temps_allocator(
      3*Branch::heap_size(dims)
      + 6*Rect::heap_size(dims)
    ),
      // branches
      m_insert_branch(dims, &m_temps_allocator),
      m_insert_rect_rec_branch(dims, &m_temps_allocator),
      m_insert_rect_branch(dims, &m_temps_allocator),
      // rects
      m_remove_rect(dims, &m_temps_allocator),
      m_search_rect(dims, &m_temps_allocator),
      m_pick_branch_rect(dims, &m_temps_allocator),
      m_pick_seeds_rect(dims, &m_temps_allocator),
      m_choose_partition_rect0(dims, &m_temps_allocator),
      m_choose_partition_rect1(dims, &m_temps_allocator),
      // only PartitionVars doesn't need an allocator (has its own)
      m_split_parition_vars(dims)
{
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

// DRTREE_TEMPLATE
// DRTREE_QUAL::drtree(const drtree &other) : drtree() {
//   CopyRec(m_root, other.m_root);
// }

DRTREE_TEMPLATE
DRTREE_QUAL::~drtree() {
  Reset(); // Free, or reset node memory
}

DRTREE_TEMPLATE
void DRTREE_QUAL::Insert(const ELEMTYPE *a_min, const ELEMTYPE *a_max,
                        const DATATYPE &a_dataId) {
#ifdef _DEBUG
  for (int index = 0; index < dims; ++index) {
    ASSERT(a_min[index] <= a_max[index]);
  }
#endif //_DEBUG

  m_insert_branch.m_data = a_dataId;
  m_insert_branch.m_child = NULL;

  for (unsigned int axis = 0; axis < dims; ++axis) {
    RECT_MIN_REF(m_insert_branch.m_rect, axis) = a_min[axis];
    RECT_MAX_REF(m_insert_branch.m_rect, axis) = a_max[axis];
  }

  InsertRect(m_insert_branch, &m_root, 0);
  count++;
}

DRTREE_TEMPLATE
int DRTREE_QUAL::Remove(const ELEMTYPE *a_min, const ELEMTYPE *a_max,
                        Callback predicate) {
#ifdef _DEBUG
  for (int index = 0; index < dims; ++index) {
    ASSERT(a_min[index] <= a_max[index]);
  }
#endif //_DEBUG

  for (unsigned int axis = 0; axis < dims; ++axis) {
    RECT_MIN_REF(m_remove_rect, axis) = a_min[axis];
    RECT_MAX_REF(m_remove_rect, axis) = a_max[axis];
  }

  int removedCount = 0;
  RemoveRect(&m_root, &m_remove_rect, removedCount, predicate);
  count -= removedCount;
  return removedCount;
}

DRTREE_TEMPLATE
int DRTREE_QUAL::Search(const ELEMTYPE *a_min, const ELEMTYPE *a_max,
                       Callback callback) {
#ifdef _DEBUG
  for (int index = 0; index < dims; ++index) {
    ASSERT(a_min[index] <= a_max[index]);
  }
#endif //_DEBUG

  for (unsigned int axis = 0; axis < dims; ++axis) {
    RECT_MIN_REF(m_search_rect, axis) = a_min[axis];
    RECT_MAX_REF(m_search_rect, axis) = a_max[axis];
  }

  // NOTE: May want to return search result another way, perhaps returning the
  // number of found elements here.

  int foundCount = 0;
  Search(m_root, &m_search_rect, foundCount, callback);

  return foundCount;
}

DRTREE_TEMPLATE
int DRTREE_QUAL::Count() {
  return count;
}


DRTREE_TEMPLATE
void DRTREE_QUAL::MoveChildren(Node* a_node, const ELEMTYPE *a_offset) {
  ASSERT(a_node);
  ASSERT(a_node->m_level >= 0);
  ASSERT(a_offset);

  for (int index = 0; index < a_node->m_count; ++index) {
    Branch &branch = a_node->m_branch[index];
    Rect &rect = branch.m_rect;
    for (unsigned int i = 0; i < dims; i++) {
      RECT_MIN_REF(rect, i) += a_offset[i];
      RECT_MAX_REF(rect, i) += a_offset[i];
    }
    if (a_node->IsInternalNode()) {
      MoveChildren(branch.m_child, a_offset);
    }
  }
}

DRTREE_TEMPLATE
void DRTREE_QUAL::RemoveAll() {
  // Delete all existing nodes
  Reset();
  count = 0;

  m_root = AllocNode();
  m_root->m_level = 0;
}

DRTREE_TEMPLATE
void DRTREE_QUAL::MoveAll(const ELEMTYPE *a_offset) {
  MoveChildren(m_root, a_offset);
}

DRTREE_TEMPLATE
void DRTREE_QUAL::Reset() {
  RemoveAllRec(m_root);
}

DRTREE_TEMPLATE
void DRTREE_QUAL::RemoveAllRec(Node *a_node) {
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

DRTREE_TEMPLATE
typename DRTREE_QUAL::Node *DRTREE_QUAL::AllocNode() {
  Node *newNode;
  newNode = new Node(this);
  InitNode(newNode);
  return newNode;
}

DRTREE_TEMPLATE
void DRTREE_QUAL::FreeNode(Node *a_node) {
  ASSERT(a_node);

  delete a_node;
}

// Allocate space for a node in the list used in DeletRect to
// store Nodes that are too empty.
DRTREE_TEMPLATE
typename DRTREE_QUAL::ListNode *DRTREE_QUAL::AllocListNode() {
  return new ListNode;
}

DRTREE_TEMPLATE
void DRTREE_QUAL::FreeListNode(ListNode *a_listNode) {
  delete a_listNode;
}

DRTREE_TEMPLATE
void DRTREE_QUAL::InitNode(Node *a_node) {
  a_node->m_count = 0;
  a_node->m_level = -1;
}

// Inserts a new data rectangle into the index structure.
// Recursively descends tree, propagates splits back up.
// Returns 0 if node was not split.  Old node updated.
// If node was split, returns 1 and sets the pointer pointed to by
// new_node to point to the new node.  Old node updated to become one of two.
// The level argument specifies the number of steps up from the leaf
// level to insert; e.g. a data rectangle goes in at level = 0.
DRTREE_TEMPLATE
bool DRTREE_QUAL::InsertRectRec(const Branch &a_branch, Node *a_node,
                               Node **a_newNode, int a_level) {
  ASSERT(a_node && a_newNode);
  ASSERT(a_level >= 0);
  ASSERT(a_level <= a_node->m_level);

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
      combine_rects(&(a_node->m_branch[index].m_rect), &a_branch.m_rect, &(a_node->m_branch[index].m_rect));
      return false;
    } else {
      // Child was split. The old branches are now re-partitioned to two nodes
      // so we have to re-calculate the bounding boxes of each node
      node_cover(&a_node->m_branch[index].m_rect, a_node->m_branch[index].m_child);

      m_insert_rect_rec_branch.m_child = otherNode;
      node_cover(&m_insert_rect_rec_branch.m_rect, otherNode);

      // The old node is already a child of a_node. Now add the newly-created
      // node to a_node as well. a_node might be split because of that.
      return AddBranch(&m_insert_rect_rec_branch, a_node, a_newNode);
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
DRTREE_TEMPLATE
bool DRTREE_QUAL::InsertRect(const Branch &a_branch, Node **a_root,
                            int a_level) {
  ASSERT(a_root);
  ASSERT(a_level >= 0 && a_level <= (*a_root)->m_level);
#ifdef _DEBUG
  for (int index = 0; index < dims; ++index) {
    ASSERT(RECT_MIN_REF(a_branch.m_rect, index) <= a_branch.m_rect.m_max[index]);
  }
#endif //_DEBUG

  Node *newNode;

  if (InsertRectRec(a_branch, *a_root, &newNode, a_level)) // Root split
  {
    // Grow tree taller and new root
    Node *newRoot = AllocNode();
    newRoot->m_level = (*a_root)->m_level + 1;

    // add old root node as a child of the new root
    node_cover(&m_insert_rect_branch.m_rect, *a_root);
    m_insert_rect_branch.m_child = *a_root;
    AddBranch(&m_insert_rect_branch, newRoot, NULL);

    // add the split node as a child of the new root
    node_cover(&m_insert_rect_branch.m_rect, newNode);
    m_insert_rect_branch.m_child = newNode;
    AddBranch(&m_insert_rect_branch, newRoot, NULL);

    // set the new root as the root node
    *a_root = newRoot;

    return true;
  }

  return false;
}

// Find the smallest rectangle that includes all rectangles in branches of a
// node.
DRTREE_TEMPLATE
void DRTREE_QUAL::node_cover(Rect* dst, Node *a_node) {
  ASSERT(a_node);
  *dst = a_node->m_branch[0].m_rect;
  for (int index = 1; index < a_node->m_count; ++index) {
    combine_rects(dst, dst, &(a_node->m_branch[index].m_rect));
  }
}

// Add a branch to a node.  Split the node if necessary.
// Returns 0 if node not split.  Old node updated.
// Returns 1 if node split, sets *new_node to address of new node.
// Old node updated, becomes one of two.
DRTREE_TEMPLATE
bool DRTREE_QUAL::AddBranch(const Branch *a_branch, Node *a_node,
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
DRTREE_TEMPLATE
void DRTREE_QUAL::DisconnectBranch(Node *a_node, int a_index) {
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
DRTREE_TEMPLATE
int DRTREE_QUAL::PickBranch(const Rect *a_rect, Node *a_node) {
  ASSERT(a_rect && a_node);

  bool firstTime = true;
  ELEMTYPEREAL increase;
  ELEMTYPEREAL bestIncr = (ELEMTYPEREAL)-1;
  ELEMTYPEREAL area;
  ELEMTYPEREAL bestArea;
  int best = 0;

  for (int index = 0; index < a_node->m_count; ++index) {
    Rect *curRect = &a_node->m_branch[index].m_rect;
    area = CalcRectVolume(curRect);
    combine_rects(&m_pick_branch_rect, a_rect, curRect);
    increase = CalcRectVolume(&m_pick_branch_rect) - area;
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
DRTREE_TEMPLATE
void DRTREE_QUAL::combine_rects(const Rect *dst, const Rect *a, const Rect *b) {
  ASSERT(dst && a && b);

  for (unsigned int index = 0; index < dims; index++) {
    RECTP_MIN_REF(dst, index) = Min(RECTP_MIN_REF(a,index), RECTP_MIN_REF(b, index));
    RECTP_MAX_REF(dst, index) = Max(RECTP_MAX_REF(a,index), RECTP_MAX_REF(b, index));
  }
}



// Split a node.
// Divides the nodes branches and the extra one between two nodes.
// Old node is one of the new ones, and one really new one is created.
// Tries more than one method for choosing a partition, uses best result.
DRTREE_TEMPLATE
void DRTREE_QUAL::SplitNode(Node *a_node, const Branch *a_branch,
                           Node **a_newNode) {
  ASSERT(a_node);
  ASSERT(a_branch);

  // Load all the branches into a buffer, initialize old node
  GetBranches(a_node, a_branch, &m_split_parition_vars);

  // Find partition
  ChoosePartition(&m_split_parition_vars, MINNODES);

  // Create a new node to hold (about) half of the branches
  *a_newNode = AllocNode();
  (*a_newNode)->m_level = a_node->m_level;

  // Put branches from buffer into 2 nodes according to the chosen partition
  a_node->m_count = 0;
  LoadNodes(a_node, *a_newNode, &m_split_parition_vars);

  ASSERT((a_node->m_count + (*a_newNode)->m_count) == m_split_parition_vars.m_total);
}

// Calculate the n-dimensional volume of a rectangle
DRTREE_TEMPLATE
ELEMTYPEREAL DRTREE_QUAL::RectVolume(Rect *a_rect) {
  ASSERT(a_rect);

  ELEMTYPEREAL volume = (ELEMTYPEREAL)1;

  for (int index = 0; index < dims; ++index) {
    volume *= RECTP_MAX_REF(a_rect, index) - RECTP_MIN_REF(a_rect, index);
  }

  ASSERT(volume >= (ELEMTYPEREAL)0);

  return volume;
}

// The exact volume of the bounding sphere for the given Rect
DRTREE_TEMPLATE
ELEMTYPEREAL DRTREE_QUAL::RectSphericalVolume(Rect *a_rect) {
  ASSERT(a_rect);

  ELEMTYPEREAL sumOfSquares = (ELEMTYPEREAL)0;
  ELEMTYPEREAL radius;

  for (unsigned int index = 0; index < dims; ++index) {
    const ELEMTYPEREAL halfExtent = (RECTP_MAX_REF(a_rect, index) -
                                     RECTP_MIN_REF(a_rect, index)) *
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
DRTREE_TEMPLATE
ELEMTYPEREAL DRTREE_QUAL::CalcRectVolume(Rect *a_rect) {
#ifdef DRTREE_USE_SPHERICAL_VOLUME
  return RectSphericalVolume(a_rect); // Slower but helps certain merge cases
#else                                 // DRTREE_USE_SPHERICAL_VOLUME
  return RectVolume(a_rect); // Faster but can cause poor merges
#endif                                // DRTREE_USE_SPHERICAL_VOLUME
}

// Load branch buffer with branches from full node plus the extra branch.
DRTREE_TEMPLATE
void DRTREE_QUAL::GetBranches(Node *a_node, const Branch *a_branch,
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
  // Note: without implementing Rect copy assignment,
  // when modifying m_coverSplit, it would also modify m_branchBuf[0].m_rect
  a_parVars->m_coverSplit = a_parVars->m_branchBuf[0].m_rect;
  for (int index = 1; index < MAXNODES + 1; ++index) {
    combine_rects(&a_parVars->m_coverSplit, &a_parVars->m_coverSplit,
                  &a_parVars->m_branchBuf[index].m_rect);
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
DRTREE_TEMPLATE
void DRTREE_QUAL::ChoosePartition(PartitionVars *a_parVars, int a_minFill) {
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
        combine_rects(&m_choose_partition_rect0, curRect, &a_parVars->m_cover[0]);
        combine_rects(&m_choose_partition_rect1, curRect, &a_parVars->m_cover[1]);

        ELEMTYPEREAL growth0 = CalcRectVolume(&m_choose_partition_rect0) - a_parVars->m_area[0];
        ELEMTYPEREAL growth1 = CalcRectVolume(&m_choose_partition_rect1) - a_parVars->m_area[1];
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
DRTREE_TEMPLATE
void DRTREE_QUAL::LoadNodes(Node *a_nodeA, Node *a_nodeB,
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
DRTREE_TEMPLATE
void DRTREE_QUAL::InitParVars(PartitionVars *a_parVars, int a_maxRects,
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

DRTREE_TEMPLATE
void DRTREE_QUAL::PickSeeds(PartitionVars *a_parVars) {
  int seed0 = 0, seed1 = 0;
  ELEMTYPEREAL worst, waste;
  ELEMTYPEREAL area[MAXNODES + 1];

  for (int index = 0; index < a_parVars->m_total; ++index) {
    area[index] = CalcRectVolume(&a_parVars->m_branchBuf[index].m_rect);
  }

  worst = -a_parVars->m_coverSplitArea - 1;

  for (int indexA = 0; indexA < a_parVars->m_total - 1; ++indexA) {
    for (int indexB = indexA + 1; indexB < a_parVars->m_total; ++indexB) {
      combine_rects(&m_pick_seeds_rect, &a_parVars->m_branchBuf[indexA].m_rect,
                    &a_parVars->m_branchBuf[indexB].m_rect);
      waste = CalcRectVolume(&m_pick_seeds_rect) - area[indexA] - area[indexB];
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
DRTREE_TEMPLATE
void DRTREE_QUAL::Classify(int a_index, int a_group, PartitionVars *a_parVars) {
  ASSERT(a_parVars);
  ASSERT(PartitionVars::NOT_TAKEN == a_parVars->m_partition[a_index]);

  a_parVars->m_partition[a_index] = a_group;

  // Calculate combined rect
  if (a_parVars->m_count[a_group] == 0) {
    a_parVars->m_cover[a_group] = a_parVars->m_branchBuf[a_index].m_rect;
  } else {
    combine_rects(&a_parVars->m_cover[a_group], &a_parVars->m_cover[a_group], &a_parVars->m_branchBuf[a_index].m_rect);
  }

  // Calculate volume of combined rect
  a_parVars->m_area[a_group] = CalcRectVolume(&a_parVars->m_cover[a_group]);

  ++a_parVars->m_count[a_group];
}

// Delete a data rectangle from an index structure.
// Pass in a pointer to a Rect, the tid of the record, ptr to ptr to root node.
// Returns true if something removed
// RemoveRect provides for eliminating the root.
DRTREE_TEMPLATE
bool DRTREE_QUAL::RemoveRect(Node **a_root, Rect *a_rect, int &a_foundCount, Callback predicate ) {
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
DRTREE_TEMPLATE
bool DRTREE_QUAL::RemoveRectRec(Node *a_node, Rect *a_rect, int &a_removedCount, Callback predicate, ListNode **a_listNode) {
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
            node_cover(&a_node->m_branch[index].m_rect, a_node->m_branch[index].m_child);
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
        if (predicate(branch.m_data, RECT_MIN(branch.m_rect), RECT_MAX(branch.m_rect))) {
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
DRTREE_TEMPLATE
bool DRTREE_QUAL::Overlap(Rect *a_rectA, Rect *a_rectB) const {
  ASSERT(a_rectA && a_rectB);

  for (unsigned int index = 0; index < dims; ++index) {
    if (RECTP_MIN_REF(a_rectA, index) > RECTP_MAX_REF(a_rectB, index) ||
        RECTP_MIN_REF(a_rectB, index) > RECTP_MAX_REF(a_rectA, index)) {
      return false;
    }
  }
  return true;
}

// Add a node to the reinsertion list.  All its branches will later
// be reinserted into the index structure.
DRTREE_TEMPLATE
void DRTREE_QUAL::ReInsert(Node *a_node, ListNode **a_listNode) {
  ListNode *newListNode;

  newListNode = AllocListNode();
  newListNode->m_node = a_node;
  newListNode->m_next = *a_listNode;
  *a_listNode = newListNode;
}

// Search in an index tree or subtree for all data retangles that overlap the
// argument rectangle.
DRTREE_TEMPLATE
bool DRTREE_QUAL::Search(Node *a_node, Rect *a_rect, int &a_foundCount,
                        Callback callback)  {
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

        if (!callback(data, RECT_MIN(branch.m_rect), RECT_MAX(branch.m_rect))) {
          return false; // Don't continue searching
        }
      }
    }
  }

  return true; // Continue searching
}

DRTREE_TEMPLATE
std::vector<typename DRTREE_QUAL::Rect> DRTREE_QUAL::ListTree() const {
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

// this could be const but I'm modifying a member variable (for avoiding allocations)
DRTREE_TEMPLATE
typename DRTREE_QUAL::PublicRect DRTREE_QUAL::Bounds() {
  ASSERT(m_root);
  ASSERT(m_root->m_level >= 0);
  PublicRect rect{this};

  if(m_root->m_count == 0) {
    return rect;
  }

  Branch &first_branch = m_root->m_branch[0];
  // init
  for(unsigned int i =0; i<dims; i++) {
    rect.low[i] = RECT_MIN_REF(first_branch.m_rect, i);
    rect.high[i] = RECT_MAX_REF(first_branch.m_rect, i);
  }
  for (int branch_id = 1; branch_id < m_root->m_count; branch_id++) {
    Branch &branch = m_root->m_branch[branch_id];
    Rect &other_rect = branch.m_rect;
    for (unsigned int index = 0; index < dims; index++) {
      rect.low[index] = Min(rect.low[index], RECT_MIN_REF(other_rect, index));
      rect.high[index] = Max(rect.high[index], RECT_MAX_REF(other_rect, index));
    }
  }

  return rect;
}

DRTREE_TEMPLATE
int DRTREE_QUAL::Dimensions() const { return this->dims; }

// for the copy ctor
DRTREE_TEMPLATE
void DRTREE_QUAL::CopyRec(Node *current, Node *other) {
  current->m_level = other->m_level;
  current->m_count = other->m_count;

  if (current->IsInternalNode()) // not a leaf node
  {
    for (int index = 0; index < current->m_count; ++index) {
      Branch *currentBranch = &current->m_branch[index];
      Branch *otherBranch = &other->m_branch[index];

      currentBranch->m_rect = otherBranch->m_rect;
      currentBranch->m_child = AllocNode();
      CopyRec(currentBranch->m_child, otherBranch->m_child);
    }
  } else // A leaf node
  {
    for (int index = 0; index < current->m_count; ++index) {
      Branch *currentBranch = &current->m_branch[index];
      Branch *otherBranch = &other->m_branch[index];

      currentBranch->m_rect = otherBranch->m_rect;
      currentBranch->m_data = otherBranch->m_data;
    }
  }
}

#undef DRTREE_TEMPLATE
#undef DRTREE_QUAL

#undef RECT_MIN
#undef RECT_MAX

#undef RECTP_MIN
#undef RECTP_MAX

#undef RECT_MIN_REF
#undef RECT_MAX_REF

#undef RECTP_MIN_REF
#undef RECTP_MAX_REF