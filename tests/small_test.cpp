#include<iostream>
#include<vector>
#include<bats.hpp>
#define F2 ModP<int, 2>

int main(int argc, char const *argv[])
{
    using VT = SparseVector<F2, size_t>;

    auto v =  VT({1,3,5,7},{1,1,1,1});
    v.insert_rows({0,1,6,11,12});
    std::cout << "\nthe final answer is" << std::endl;
    v.print();
    return 0;
}
