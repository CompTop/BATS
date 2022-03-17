#include<iostream>
#include<vector>
#include<bats.hpp>
#define F2 ModP<int, 3>

int main(int argc, char const *argv[])
{
    using VT = SparseVector<F2, size_t>;

    auto v =  VT({1,3,5,7},{1,1,1,1});
    v.print();
    v.insert_rows({0,1,6,11,12});
    std::cout << "\nthe final answer is" << std::endl;
    v.print();

    v.insert_rows({4,6}, {2,2});
    std::cout << "\nthe final answer is" << std::endl;
    v.print();

    std::vector<size_t> perm = {5};
    bats::print_1D_vectors(perm);
    bats::print_1D_vectors(bats::perm_to_the_end(perm, 6));
    return 0;
}
