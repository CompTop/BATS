#pragma once
// utilities for simplices
#include <functional> // for hash
#include <vector>
#include <sstream>

#include "permutation.hpp"

namespace bats {
namespace util {

// implement has function for vectors

struct SimplexHasher
{
    std::size_t operator()(const std::vector<size_t>& k) const
    {
        using std::size_t;
        using std::hash;

        // Compute individual hash values for first,
        // second and third and combine them using XOR
        // and bit shifting:
        size_t ret = 0;
        for (auto i : k) {
            // ret = ret ^ (i);
            ret *= i;
            // ret += ((i+1) << 6);
        }
        return ret;
    }
};


class SimplexContainer {
private:
    std::vector<size_t> data;
    size_t k; // number of vertices in this simplex size
public:

    // empty constructor - not available
    // constructor with dimension specified
    SimplexContainer(size_t d) : k(d+1) {}

    // emplace_back with vector
    void emplace_back(std::vector<size_t> &s) {
        if (data.size() != k) {throw "unexpected simplex size!";}
        for (auto i : s) {
            data.emplace_back(i);
        }
    }

    // size - number of simplices
    size_t size() const {
        return data.size() / k;
    }
    // operator[] - return vector
    // reserve - reserve space for a fixed number of simplices
    // dim - dimension
    inline size_t dim() const {
        return k-1;
    }
    // function to return iterator over vertices of simplex i
};

// returns true if simplex is degenerate
// assume s is in lexicographic order
template <typename T>
bool is_degenerate(const std::vector<T> &s) {
    for (size_t i = 1; i < s.size(); i++) {
        if (s[i] == s[i-1]) { return true; }
    }
    return false;
}

// // sort simplex s into lexicographic order
// // return the sign of the sort permutation
// // returns 0 if simplex is degenerate
template <typename T>
int simplex_sign(std::vector<T> &s) {
    auto p = sortperm(s);
    apply_perm(s, p);
    return is_degenerate(s) ? 0 : perm_sign(p);
}

template <typename IO, typename T>
void write_simplex(IO &io, std::vector<T> &s) {
    for (auto x : s) {
        io << x << ", ";
    }
    io << '\n';
}

template <typename IO, typename TI>
void write_simplex(IO &io, TI&& begin, TI&& end) {
    auto it = begin;
    while (it != end) {
        io << *it++ << ",";
    }
    io << '\n';
}

template < typename T>
void read_simplex(std::string &line, std::vector<T> &s) {
    s.clear();
    // read csv line to string
    std::string token;
    std::istringstream iss(line);
    while (getline(iss, token, ',')) {
        // std::cout << token << ',' << std::endl;
        if (token.size() > 0) {
            s.emplace_back(std::stoi(token));
        }
    }
}

} // namespace util
} // namespace bats
