#pragma once
#include <iostream>
using namespace std;

// This is special linked-list used for column of Vineyard matrix, where
// row index is shared by all lists
class Node { 
public:
    int data;
    size_t* index_ptr; // row index
    Node* next;
    // Defalt constructor
    Node() : data(0), index_ptr(nullptr), next(nullptr) {
        // cout << "Constructed default node\n";
    }
    // copy constructor 
    Node(const Node &node2) : data(node2.data){
        // cout << "Copy constructor for data" << data <<  "\n";
        index_ptr = node2.index_ptr;
        next = node2.next;
    }


    Node(int data, size_t* ptr) : data(data), index_ptr(ptr), next(nullptr) {
        // cout << "Constructed node: " << data << "\n";
    }
};

class LinkedList{
public:
    LinkedList() { // constructor
        head = nullptr;
    }

    // iterator constructor by
    // passing into sparse vector and row index pointer
    template <typename IT1, typename IT2>
    LinkedList(IT1 begin_it, IT1 end_it, IT2 row_it) {
        head = nullptr;
        for (auto it = begin_it; it < end_it; ++it) {
            // std::cout << "(*it).val = "<< (*it).val << ", row_it+(*it).ind " << row_it+(*it).ind << std::endl;
            insert((*it).val, *(row_it+(*it).ind));
        }
	}

    // copy constructor used for deepcopy
    LinkedList(const LinkedList &ll2) {
        // std::cout << "copy constructor is called." << std::endl;
        if (ll2.head == nullptr) {
            head = nullptr;
        }else{
            Node* temp2 = ll2.head;
            Node* next2 = temp2->next;
            head = new Node(temp2->data, temp2->index_ptr);
            Node* temp = head; // for new linkedlist
            while (temp2->next != nullptr){
                Node* new_node = new Node(next2->data, next2->index_ptr);
                temp->next = new_node;
                temp = temp->next;
                temp2 = next2;
                next2 = temp2->next;
            }
        }
    }

    // move constructor (constructing new object, where head now is null)
    LinkedList(LinkedList&& ll2) {
        // std::cout << "move constructor is called." << std::endl;
        if (ll2.head == nullptr) 
            head = nullptr;
        else{
            head = ll2.head; // steal temp obj
            ll2.head = nullptr; // nullize temp obj
        }
    }
    // move assignment operator (replace head by other object, where head might not be null)
    LinkedList& operator=(LinkedList&& other){
        // std::cout << "move assignment operator is called." << std::endl;
        if (this != &other){
            delete head; // delete current object first 
            if (other.head == nullptr) 
                head = nullptr;
            else{
                head = other.head; // steal temp obj
                other.head = nullptr; // nullize temp obj
            }
        }
       
        return *this;
    }

    // LinkedList& operator=(const LinkedList& other)
    // {
    //     s = other.s;
    //     std::cout << "copy assigned\n";
    //     return *this;
    // }

    // destructor
    ~LinkedList() {
        // std::cout << "call the destructor of LinkedList" << std::endl;
        clear();
    }
    void insert(int val, size_t* ptr);
    void display();
    void clear();
    int getval(size_t i) const;

    // TODO: find y <- ax + y
    // void axpy(int a, LinkedList x){
    // }
private:
    Node* head;
};

// Add node to a list
void LinkedList::insert(int val, size_t* ptr) {
    Node* newnode = new Node(val, ptr);
    if (head == nullptr) {
        head = newnode;
    }
    else {
        Node* temp = head; // head is not nullptr
        while (temp->next != nullptr) { 
            temp = temp->next; // go to end of list
        }
        temp->next = newnode; // linking to newnode
    }
}

// Delete the entire list
void LinkedList::clear() {
    Node* next; // next node 
    Node* temp = head; // current node
    while (temp != nullptr) {
        next = temp->next;
        delete temp;
        temp = next;
    }
    head = nullptr;
}

// get value give a row index
int LinkedList::getval (const size_t i) const{
    Node* temp = head;
    // iterate over all nodes, might be inefficient
    while (temp != nullptr) { 
        if (*(temp->index_ptr) == i){
            return temp->data;
        }
        temp = temp->next;
    }
    return 0;
}

// function to display the entire list
void LinkedList::display() {
    if (head == nullptr) {
        cout << "List is empty!" << endl;
    }
    else {
        Node* temp = head;
        while (temp != nullptr) {
            cout << "("<< temp->data << ", " << *(temp->index_ptr)<< ") ";
            temp = temp->next;
        }
        cout << endl;
    }
}