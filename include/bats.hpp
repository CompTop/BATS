#pragma once
/*
Single header file to include all functionality
*/
#ifdef BATS_OPCOUNT
namespace bats {
    unsigned global_ops;
}
#endif

#include "linalg.hpp"
#include "complex.hpp"
#include "filtration.hpp"
#include "chain.hpp"
#include "homology.hpp"
#include "persistence.hpp"
#include "multigraph.hpp"
#include "quiver.hpp"
#include "topology.hpp"
#include "future.hpppp"
