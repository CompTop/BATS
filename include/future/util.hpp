#pragma once

namespace bats {
namespace future {

template <typename RandomAccessIterator>
ssize_t find_pivot_high(const RandomAccessIterator &v, ssize_t start, ssize_t end) {
	for (ssize_t i = start; i < end; i++) {
		if (v[i] != 0) {
			return i;
		}
	}
	return -1;
}
template <typename RandomAccessIterator>
inline ssize_t find_pivot_high(const RandomAccessIterator &v, ssize_t end) {return find_pivot_high(v, 0, end);}

}
}
