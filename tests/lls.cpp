#include <iostream>
#include <vector>
#include "linkded_list.hpp"
using namespace std;

int main() {
    LinkedList<int> list;
    std::cout << "Creating List\n";
    std::shared_ptr<size_t> pt1 = std::make_shared<size_t>(0);
    list.insert(99, pt1);

    std::shared_ptr<size_t> pt2 = std::make_shared<size_t>(2);
    list.insert(2, pt2);

    std::shared_ptr<size_t> pt3 = std::make_shared<size_t>(4);
    list.insert(3, pt3);
    std::cout << "Linked List 1 data:\n";
    list.display();
    LinkedList<int> list2(list);
    std::cout << "Linked List 2 data:\n";
    list2.display();

    std::cout << "Clearing list 1\n";
    list.clear();
    std::cout << "Linked List1 data:\n";
    list.display();

    std::cout << "Linked List 2 data:\n";
    list2.display();

    std::cout << "Call the move " << std::endl;
    // LinkedList<int> list3;
    // std::shared_ptr<size_t> pt4 = std::make_shared<size_t>(1);
    // list3.insert(1, pt4);
    // std::shared_ptr<size_t> pt5 = std::make_shared<size_t>(2);
    // list3.insert(2, pt5);

    // list3 = std::move(list2);
    LinkedList<int> list3(std::move(list2));
    list3.display();
    list2.display();
    return 0;
}