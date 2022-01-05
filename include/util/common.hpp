#pragma once
#include <cstddef>
#include <limits>


namespace bats {

// no index
const size_t NO_IND = std::numeric_limits<size_t>::max();

// get global ops
#ifdef BATS_OPCOUNT
unsigned get_global_ops() {return field_ops;}
void reset_global_ops() {field_ops = 0;}
#endif

}
