#include <iostream>
#include <vector>
#include "linkded_list.hpp"
using namespace std;

int main() {
    LinkedList<int> list;
    std::cout << "Creating List\n";
    std::shared_ptr<size_t> pt1 = std::make_shared<size_t>(0);
    list.insert(0, pt1);

    std::shared_ptr<size_t> pt2 = std::make_shared<size_t>(2);
    list.insert(0, pt2);

    std::shared_ptr<size_t> pt3 = std::make_shared<size_t>(4);
    list.insert(3, pt3);

    std::shared_ptr<size_t> pt6 = std::make_shared<size_t>(5);
    list.insert(0, pt6);

    std::shared_ptr<size_t> pt4 = std::make_shared<size_t>(6);
    list.insert(5, pt4);

    std::shared_ptr<size_t> pt5 = std::make_shared<size_t>(7);
    list.insert(0, pt5);

    std::cout << "Linked List 1 data:\n";
    list.display();

    std::cout << "Clear zeros of Linked List 1:\n";
    list.clear_zeros();
    list.display();

    // LinkedList<int> list2(list);
    // std::cout << "Linked List 2 data:\n";
    // list2.display();

    // std::cout << "Clearing list 1\n";
    // list.clear();
    // std::cout << "Linked List1 data:\n";
    // list.display();

    // std::cout << "Linked List 2 data:\n";
    // list2.display();

    // std::cout << "Call the move " << std::endl;

    // // list3 = std::move(list2);
    // LinkedList<int> list3(std::move(list2));
    // list3.display();
    // list2.display();
    return 0;
}