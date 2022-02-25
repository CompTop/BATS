#pragma once
/*
Single header file to include all functionality
*/
#ifdef BATS_OPCOUNT
namespace bats {
    unsigned field_ops; // global counter for field operations
    unsigned column_ops; // global counter for column operations
}
#endif

#include "linalg.hpp"
#include "complex.hpp"
#include "filtration.hpp"
#include "dgvs.hpp"
#include "chain.hpp"
#include "homology.hpp"
#include "persistence.hpp"
#include "zigzag.hpp"
#include "multigraph.hpp"
#include "quiver.hpp"
#include "topology.hpp"
// #include "future.hpp"
