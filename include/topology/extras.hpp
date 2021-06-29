#pragma once
/*
Some constructions for demonstrations
*/

#include <vector>
#include <set>

#include <complex/cell_complex.hpp>
#include <complex/cell_map.hpp>

namespace bats {

struct triangle {
    size_t a, b, c, ab, bc, ca;

    triangle() {}
    // node only constructor
    triangle(size_t a, size_t b, size_t c) : a(a), b(b), c(c) {}

};

// diagram up to kth level
Diagram<CellComplex, CellularMap> serpinski_diagram(
    size_t k
) {

    using Vint = SparseVector<int>;
    using MTint = ColumnMatrix<Vint>;

    Diagram<CellComplex, CellularMap> X(k, k-1);
    // create first iteration - just a triangle
    CellComplex S0;
    S0.add_vertices(3);
    for (size_t i = 0; i < 3; i++) {
        S0.add({i, (i+1) % 3}, {-1, 1}, 1);
    }
    X.set_node(0, S0);

    std::vector<triangle> old_tri(1);
    old_tri[0].a = 0;
    old_tri[0].b = 1;
    old_tri[0].c = 2;
    old_tri[0].ab = 0;
    old_tri[0].bc = 1;
    old_tri[0].ca = 2;

    // allocations for loop
    std::vector<triangle> new_tri;
    CellComplex Si;
    CellularMap Fi(2);
    std::vector<Vint> tri_bdry;
    std::vector<Vint> col0;
    std::vector<Vint> col1;
    std::vector<Vint> col2;

    std::vector<size_t> ind;
    std::vector<int> val;

    for (size_t i = 1; i < k; i++) {
        new_tri.clear();

        col0.resize(S0.ncells(0));
        col1.resize(S0.ncells(1));
        col2.resize(S0.ncells(2));

        Si.add_vertices(S0.ncells(0));

        for (auto tri : old_tri) {
            size_t a = tri.a;
            size_t b = tri.b;
            size_t c = tri.c;

            // add new vertices
            size_t d = Si.add_vertex();
            size_t e = Si.add_vertex();
            size_t f = Si.add_vertex();

            // new triangles
            auto t1 = triangle(a, e, d);
            auto t2 = triangle(e, b, f);
            auto t3 = triangle(d, f, c);

            // add new edges
            val = {-1, 1};
            t1.ab = Si.add({a, e}, val, 1);
            t1.bc = Si.add({e, d}, val, 1);
            t1.ca = Si.add({d, a}, val, 1);

            t2.ab = Si.add({e, b}, val, 1);
            t2.bc = Si.add({b, f}, val, 1);
            t2.ca = Si.add({f, e}, val, 1);

            t3.ab = Si.add({d, f}, val, 1);
            t3.bc = Si.add({f, c}, val, 1);
            t3.ca = Si.add({c, d}, val, 1);

            new_tri.emplace_back(t1);
            new_tri.emplace_back(t2);
            new_tri.emplace_back(t3);

            col1[tri.ab] = Vint({t1.ab, t2.ab}, {1, 1});
            col1[tri.bc] = Vint({t2.bc, t3.bc}, {1, 1});
            col1[tri.ca] = Vint({t2.ca, t1.ca}, {1, 1});

            tri_bdry.emplace_back(Vint({t1.bc, t2.ca, t3.ab}, {-1, -1, -1}));
        }

        for (size_t k = 0; k < S0.ncells(0); k++) {
            col0[k] = Vint(k);
        }

        Fi[0] = MTint(Si.ncells(0), S0.ncells(0), col0);
        Fi[1] = MTint(Si.ncells(1), S0.ncells(1), col1);

        // update old triangle boundaries
        for (size_t k = 0; k < S0.ncells(2); k++) {
            tri_bdry[k] = Fi[1] * tri_bdry[k];
            col2[k] = Vint(k);
        }
        for (size_t k = 0; k < tri_bdry.size(); k++) {
            ind.clear();
            val.clear();
            for (auto it = tri_bdry[k].nzbegin(); it < tri_bdry[k].nzend(); it++) {
                ind.emplace_back((*it).ind);
                val.emplace_back((*it).val);
            }
            Si.add(ind, val, 2);
        }
        Fi[2] = MTint(Si.ncells(2), S0.ncells(2), col2);

        X.set_node(i, Si);
        X.set_edge(i-1, i-1, i, Fi);

        // prepare for next iter
        S0 = Si;
        Si = CellComplex();
        std::swap(old_tri,new_tri);

    }

    return X;
}

} // namespace bats
