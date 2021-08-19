#include <vector>
#include <iostream>
#include <bats.hpp>

int main() {

    // Simplicial complexes and Homology
    {
        bats::SimplicialComplex X;
        X.add_recursive({0,1,2});
        X.add_recursive({2,3});
        X.add({1,3});

        X.print_summary();

        using F2 = ModP<int, 2>;
        auto R = bats::Reduce(X, F2());

        R.print_summary();
    }

    // Persistent homology
    {
        bats::Filtration<double, bats::SimplicialComplex> F;
        std::vector<size_t> spx;
        spx = {0,1,2}; F.add_recursive(0.0, spx);
        spx = {2,3};   F.add_recursive(1.0, spx);
        spx = {1,3};   F.add(2.0, spx);

        F.complex().print_summary();

        using F2 = ModP<int, 2>;
        auto R = bats::Reduce(F, F2());

        for (size_t k = 0; k < R.maxdim(); ++k) {
            std::cout << "\n" << k << "-dimensional barcode:\n";
            for ( auto p : R.persistence_pairs(k)) {
                std::cout << p.str() << std::endl;
            }
        }
    }

    // Maps
    {
        bats::SimplicialComplex X;
        X.add_recursive({0,1});
        X.add_recursive({1,2});
        X.add_recursive({2,3});
        X.add_recursive({0,3});
        bats::SimplicialComplex Y = X; // copy

        X.print_summary();

        std::vector<size_t> f0{2,1,0,3}; // reflection map
        auto F = bats::SimplicialMap(X, Y, f0);

        // apply the chain functor
        using F3 = ModP<int, 3>;
        auto CX = bats::Chain(X, F3());
        auto CY = bats::Chain(Y, F3());
        auto CF = bats::Chain(F, F3());

        auto RX = bats::ReducedChainComplex(CX);
        auto RY = bats::ReducedChainComplex(CY);
        RX.print_summary();
        RY.print_summary();

        for (size_t k = 0; k < 2; k++) {
            std::cout << "\nInduced map in dimension " << k << std::endl;
            auto HF = bats::induced_map(CF, RX, RY, k);
            HF.print();
        }
    }

    // Zigzag homology
    std::cout << "\nZigzag Example\n";
    {
        bats::Diagram<bats::SimplicialComplex, bats::CellularMap> D(2,1);

        bats::SimplicialComplex X;
        X.add_recursive({0,1});
        X.add_recursive({1,2});
        X.add_recursive({2,3});
        X.add_recursive({0,3});

        std::vector<size_t> f0{2,1,0,3}; // reflection map
        auto F = bats::SimplicialMap(X, X, f0);

        D.set_node(0, X);
        D.set_node(1, X);
        D.set_edge(0, 0, 1, F); // edge 0: (0,1)

        using F3 = ModP<int, 3>;
        auto CD = bats::Chain(D, F3());

        auto HD = bats::Hom(CD, (size_t) 1); // homology in dimension 1

        auto ps = bats::barcode(HD, 1);
    	for (auto p : ps) {
    		std::cout << p.str() << std::endl;
    	}


    }


    return EXIT_SUCCESS;
}
