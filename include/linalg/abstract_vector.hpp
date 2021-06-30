#pragma once

#include <iostream>
#include <string>

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
    // nzpair(const nzpair& p) : ind(p.ind), val(p.val) {};
    nzpair(std::string &str) {
        std::string token;
		std::istringstream iss(str);
		getline(iss, token, ':');
		ind = TI(std::stoi(token));
		getline(iss, token);
		val = TV(std::stoi(token));
    }

    inline bool operator==(const nzpair &other) const { return (ind == other.ind) && (val == other.val); }
    inline bool operator!=(const nzpair &other) const { return (ind != other.ind) || (val != other.val); }
};

// define comparison
template <typename TI, typename TV>
bool operator<(const TI& a, const nzpair<TI, TV>& b) { return a < b.ind; }
template <typename TI, typename TV>
bool operator<(const nzpair<TI, TV>& a, const TI& b) { return a.ind < b; }
template <typename TI, typename TV>
bool operator<(const nzpair<TI, TV>& a, const nzpair<TI, TV>& b) { return a.ind < b.ind; }

// template <typename TI, typename TV>
// void std::swap(nzpair<TI, TV>& a, nzpair<TI, TV>& b) {
// 	TI tmp = a.ind;
// 	a.ind = b.ind;
// 	b.ind = tmp;
// 	TV tmp2 = a.val;
// 	a.val = b.val;
// 	b.val = tmp2;
// }

// printing
template <typename TI, typename TV>
std::ostream& operator<<( std::ostream& os, const nzpair<TI, TV> &x) {
os << '(' << x.ind << " : " << x.val << ')';
return os;
}
