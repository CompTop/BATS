#include <bats.hpp>
#include <vector>
#include <iostream>

using FT = ModP<int, 2>;

using CpxT = bats::LightSimplicialComplex<size_t, std::unordered_map<size_t, size_t>>;
// using CpxT = bats::SimplicialComplex;

int main() {

    /*
    First we build a simplicial complex which can be used to extend filtrations
    */
    CpxT X(3,2);

    std::vector<size_t> s;

    s = {0}; X.add(s);
    s = {1}; X.add(s);
    s = {2}; X.add(s);
    s = {0,1}; X.add(s);
    s = {1,2}; X.add(s);
    s = {0,2}; X.add(s);
    s = {0,1,2}; X.add(s);

    X.print_summary();

    /*
    Now let's exend a filtration
    */
    std::vector<double> f0 = {0.0, 0.1, 0.2};
    // lower star filtration
    std::function<double(const std::vector<size_t>&)> filtfn = [&f0](
        const std::vector<size_t>& s
    ) -> double {
        return f0[*std::max_element(s.begin(), s.end(), [&f0](size_t i, size_t j) {return f0[i] < f0[j];})];
    };
    auto vals = extend_filtration(X, filtfn);
    for (auto& valsk: vals) {
        for (auto& v: valsk) {
            std::cout << v << ", ";
        }
        std::cout << "\n";
    }

    /*
    The above is wrapped by the lower_star_filtration function
    */
    vals = lower_star_filtration(X, f0);
    for (auto& valsk: vals) {
        for (auto& v: valsk) {
            std::cout << v << ", ";
        }
        std::cout << "\n";
    }

    /*
    Now let's build a filtration and reduce
    */
    auto F = bats::Filtration(X, vals);
    auto C = bats::Chain(F, FT());
    auto R = bats::Reduce(C);
    R.print_summary();

    for (auto& p: R.persistence_pairs(0)) {
        std::cout << p.str() << std::endl;
    }

    /*
    If we want to change the filtration we have a variety of options.
    First is to update the filtration on the FilteredChainComplex
    Then, we reduce the complex again
    */
    // update filtration
    f0 = {1.1, 1.0, 1.2};
    vals = lower_star_filtration(X, f0);
    C.update_filtration(vals);
    R = bats::Reduce(C);
    R.print_summary();

    for (auto& p: R.persistence_pairs(0)) {
        std::cout << p.str() << std::endl;
    }

    /*
    The second option is to update the ReducedFilteredChainComplex
    In some situations, this may be faster
    */
    // update filtration again
    f0 = {0.0, 0.1, 0.2};
    vals = lower_star_filtration(X, f0);
    R.update_filtration(vals);
    R.print_summary();

    for (auto& p: R.persistence_pairs(0)) {
        std::cout << p.str() << std::endl;
    }

    /*
    The final option would be to just call Chain() again to build
    the updated chain complex.
    */

    return EXIT_SUCCESS;
}
