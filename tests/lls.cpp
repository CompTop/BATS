#include <iostream>
#include <vector>
#include "linkded_list.hpp"
using namespace std;

int main() {
    LinkedList list;
    cout << "Creating List\n";
    size_t ind{0};
    size_t* pt1{&ind};
    list.insert(99, pt1);

    ind = 2;
    size_t* pt2{&ind};
    list.insert(2, pt2);

    ind = 4;
    size_t* pt3{&ind};
    list.insert(3, pt3);
    cout << "Linked List 1 data:\n";
    list.display();
    LinkedList list2(list);
    cout << "Linked List 2 data:\n";
    list2.display();

    cout << "Clearing list 1\n";
    list.clear();
    cout << "Linked List1 data:\n";
    list.display();

    cout << "Linked List 2 data:\n";
    list2.display();
    return 0;
}