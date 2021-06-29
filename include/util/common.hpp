#pragma once
#include <cstddef>
#include <limits>


namespace bats {

// no index
const size_t NO_IND = std::numeric_limits<size_t>::max();

// get global ops
#ifdef BATS_OPCOUNT
unsigned get_global_ops() {return global_ops;}
void reset_global_ops() {global_ops = 0;}
#endif

}
