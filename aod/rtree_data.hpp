#include "./rtree_base.hpp"
namespace aod {

#define PRE template <class DATATYPE>
#define QUAL rtree_data<DATATYPE>

PRE class rtree_data : public rtree_base {
 private:
  std::vector<DATATYPE> m_data;
  void set_data(Did did, const DATATYPE &data) {
    ASSERT(did);
    m_data[did.id] = data;
  }

  const DATATYPE &get_data(Did did) {
    ASSERT(did);
    return m_data[did.id];
  };

 public:
  using rtree_base::rtree_base;

  using Predicate = std::function<bool(const DATATYPE &)>;
  using SearchCb = std::function<bool(const DATATYPE &)>;
  void insert(const Vec &low, const Vec &high, const DATATYPE &data) {
    const Did did = make_data_id();
    m_data.resize(m_data_count); // TODO
    set_data(did, data);
    rtree_base::insert(low, high, did);
  }

  std::vector<DATATYPE> search(const Vec &low, const Vec &high) {
    std::vector<DATATYPE> results;
    search(low, high, results);
    return results;
  }

  void search(const Vec &low, const Vec &high, std::vector<DATATYPE> &results) {
    std::vector<Did> intermediate_result;
    rtree_base::search(low, high, intermediate_result);
    results.clear();
    results.reserve(intermediate_result.size());
    for(Did did : intermediate_result) {
      results.push_back(get_data(did));
    }
  }
};

}
