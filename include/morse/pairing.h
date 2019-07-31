#include <vector>

// class that wraps pre-pairing of complex
// functions that compute morse pairings




// template over
//  TC - complex type
template <class TC>
class MorsePairing
{
private:
  std::vector<std::vector<bool>> ispaired;
  std::vector<std::vector<size_t>> up;
  std::vector<std::vector<size_t>> down;
  size_t _maxdim;
public:

  // morse pairing for complex up to maxdim cells
  MorsePairing(size_t maxdim) : ispaired(std::vector<std::vector<bool>>(maxdim+1)),
                                up(std::vector<std::vector<size_t>>(maxdim)),
                                down(std::vector<std::vector<size_t>>(maxdim)),
                                _maxdim(maxdim) {}

  inline size_t maxdim() const { return _maxdim; }

  void reserve(size_t dim, size_t n) {
    if ( ispaired[dim].size() < n ) {
      ispaired[dim].resize(n, false);
    }
    if (dim < maxdim())
      if (up[dim].size() < n) {
        up[dim].reserve(n);
      }
    }
    if (dim > 0) {
      if (down[dim-1].size() < n) {
        down[dim-1].reserve(n);
      }
    }
  }

  // set pair (i,j) in dimension dim
  void set_pair_unsafe(size_t dim, size_t i, size_t j) {
    up[dim].emplace_back(i);
    down[dim].emplace_back(j);
    ispaired[dim][i] = true;
    ispaired[dim+1][j] = true;
  }

  bool set_pair_safe(size_t dim, size_t i, size_t j) {
    if (ispaired.size() < dim + 1) {
      throw "Not enough dimensions in pairing!"
      return false;
    }
    // check that ispaired has enough entries
    if ( ispaired[dim].size() < i ) {
      ispaired[dim].resize(i, false);
    }
    if ( ispaired[dim+1].size() < j ) {
      ispaired[dim+1].resize(j, false);
    }
    if (ispaired[dim][i] || ispaired[dim+1][j]) {
      // was already paired
      return false;
    }
    set_pair_unsafe(dim, i, j);
    return true;
  }


};
