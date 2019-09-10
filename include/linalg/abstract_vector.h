#pragma once

#include <iostream>

/*
struct that holds index value pairs
*/
template <typename TI, typename TV>
struct nzpair {
    TI ind;
    TV val;

    nzpair() : ind(TI(0)), val(TV(0)) {};
    nzpair(const TI ind) : ind(ind), val(TV(0)) {};
    nzpair(const TI ind, const TV val) : ind(ind), val(val) {};
};

// define comparison
template <typename TI, typename TV>
bool operator<(const TI& a, const nzpair<TI, TV>& b) { return a < b.ind; }
template <typename TI, typename TV>
bool operator<(const nzpair<TI, TV>& a, const TI& b) { return a.ind < b; }
template <typename TI, typename TV>
bool operator<(const nzpair<TI, TV>& a, const nzpair<TI, TV>& b) { return a.ind < b.ind; }

// printing
template <typename TI, typename TV>
std::ostream& operator<<( std::ostream& os, const nzpair<TI, TV> &x) {
os << '(' << x.ind << " : " << x.val << ')';
return os;
}
