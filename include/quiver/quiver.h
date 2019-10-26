/*
A data type to hold a the companion matrix of a quiver
*/
#include <vector>



// template over implementation
template <typename I>
class Quiver {
private:
    // storage is essentially a (directed) edge list
    // with list of matrices for each edge

    std::vector<size_t> edgelist; // edges in quiver
    std::vector<size_t> dim; // dimension of each vector space
    std::vector<I> A; // matrices on edges
public:

    // default empty constructor
    Quiver() {};

    void add_edge(size_t i, size_t j, I& Aij) {
        // see if we need to extend number of nodes in quiver

        // first see if dim[i] and dim[j] have been set
        // if so, check that dimesnions of Aij agree

        // if checks pass, we can add edge
    }

}

template <typename I>
class QuiverMorphism {
private:
    std::vector<size_t> dim_src; // dimension of vector spaces in source
    std::vector<size_t> dim_targ; // dimension of vector spaces in target
    std::vector<I> M; // block matrices for diagonal

public:
    // empty constructor
    QuiverMorphism() {};

    // quiver morphism
    QuiverMorphism(
        std::vector<size_t> dim1,
        std::vector<size_t> dim2) : dim_src(dim1), dim_targ(dim)
    {}; // TODO: initialize with empty maps

    // TODO: static member to initialize identity morphism

    // operator[] returns morphism


}
