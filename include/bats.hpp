#pragma once
/*
Single header file to include all functionality
*/
#ifdef BATS_OPCOUNT
namespace bats {
    unsigned global_ops;
}
#endif

#include "linalg.h"
#include "complex.h"
#include "filtration.h"
#include "chain.h"
#include "homology.h"
#include "persistence.h"
#include "multigraph.h"
#include "quiver.h"
#include "topology.h"
#include "future.hpp"
