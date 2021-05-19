#include <bats.hpp>
#include <vector>
#include <iostream>

int main() {

    bats::RightFiltration<bats::SimplicialComplex> F;


    std::vector<size_t> spx;
    // create a cycle that persists for a while
    spx = {0,1}; F.add_recursive(0.0, 10.0, spx);
    spx = {0,2}; F.add_recursive(0.0, 10.0, spx);
    spx = {1,2}; F.add_recursive(0.0, 10.0, spx);

    {
        F.complex().print_summary();

        using F2 = ModP<int, 2>;
        auto ps = bats::barcode(F, F2());

        for (auto& pk : ps) {
            for (auto p : pk) {
                if (p.length() > 0)
                    std::cout << p.str() << std::endl;
            }
        }
    }

    // now block cycle for some period of time
    std::cout <<"\nadding block:\n";
    spx = {0,1,2}; F.add(2.0, 4.0, spx);

    {
        F.complex().print_summary();

        using F2 = ModP<int, 2>;
        auto ps = bats::barcode(F, F2());

        for (auto& pk : ps) {
            for (auto p : pk) {
                if (p.length() > 0)
                    std::cout << p << std::endl;
            }
        }
    }




    return EXIT_SUCCESS;
}
