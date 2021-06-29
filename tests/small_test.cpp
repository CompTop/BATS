#include<iostream>
#include<vector>
#include<yuan/update_information.hpp>
int main(int argc, char const *argv[])
{
    std::vector<std::vector<size_t>> nbrs = {{1,2},{7,8}};
    std::cout << nbrs[0][1] << std::endl;

    return 0;
}
