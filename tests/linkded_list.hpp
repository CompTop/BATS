#pragma once
#include <iostream>
using namespace std;

// This is special linked-list used for column of Vineyard matrix, where
// row index is shared by all lists
// TODO: current constructor implementation assumes that row index values are in ascending order, 
// please check it
template <typename TV>
class Node { 
public:
    TV data;
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


    Node(TV data, size_t* ptr) : data(data), index_ptr(ptr), next(nullptr) {
        // cout << "Constructed node: " << data << "\n";
    }
};

template <typename TV> // type of node
class LinkedList{
private:
    Node<TV>* head;
public:
    using node_type = Node<TV>;
    LinkedList() { // constructor
        head = nullptr;
    }

    // constructor from a sparse vector 
    // begin_it: begin iterator of sparse vector
    // end_it: end iterator of sparse vector.
    // row_it: begin iterator of the row indices vector
    template <typename IT1, typename IT2>
    LinkedList(IT1 begin_it, IT1 end_it, IT2 row_it) {
        head = nullptr;
        for (auto it = begin_it; it < end_it; ++it) {
            // std::cout << "(*it).val = "<< (*it).val << ", row_it+(*it).ind " << row_it+(*it).ind << std::endl;
            auto r_ptr = *(row_it+ it->ind);
            insert(it->val, r_ptr);
        }
	}

    // Constructor from a column of CSC matrix
    // rptr: row begin iterator 
    // vptr: value begin iterator 
    // row_it: begin iterator of the row indices vector
    // col_size: the number of (elements) nodes needed in the column
    template <typename IT1, typename IT2, typename IT3>
    LinkedList(IT1 rptr, IT2 vptr, IT3 row_it, size_t col_size) {
        head = nullptr;
        auto itr = rptr;
        auto itv = vptr;
        while(itr < rptr + col_size){
            if (TV(*itv) != 0){ // zero node is not needed
                insert(*itv, *(row_it + (*itr)));
            }
            itv++;
            itr++;
        }
	}

    // copy constructor used for deepcopy
    LinkedList(const LinkedList &ll2) {
        // std::cout << "copy constructor is called." << std::endl;
        if (ll2.head == nullptr) {
            head = nullptr;
        }else{
            node_type* temp2 = ll2.head;
            node_type* next2 = temp2->next;
            head = new node_type(temp2->data, temp2->index_ptr);
            node_type* temp = head; // for new linkedlist
            while (temp2->next != nullptr){
                node_type* new_node = new node_type(next2->data, next2->index_ptr);
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
            // TODO: check it delete head is enough or should we delete all nodes like clear() (memory leak) 
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
        clear();
    }
    void insert(TV val, size_t* ptr){
        node_type* newnode = new node_type(val, ptr);
        if (head == nullptr) {
            head = newnode;
        }
        else {
            node_type* temp = head; // head is not nullptr
            while (temp->next != nullptr) { 
                temp = temp->next; // go to end of list
            }
            temp->next = newnode; // linking to newnode
        }
    }
    void display(){
        if (head == nullptr) {
            cout << "List is empty!" << endl;
        }
        else {
            node_type* temp = head;
            while (temp != nullptr) {
                cout << "("<< temp->data << ", " << *(temp->index_ptr)<< ") ";
                temp = temp->next;
            }
            cout << endl;
        }
    }
    void clear(){
        node_type* next; // next node 
        node_type* temp = head; // current node
        while (temp != nullptr) {
            next = temp->next;
            delete temp;
            temp = next;
        }
        head = nullptr;
    };
    TV getval(size_t i) const{
        node_type* temp = head;
        // iterate over all nodes, might be inefficient
        while (temp != nullptr) { 
            if (*(temp->index_ptr) == i){
                return temp->data;
            }
            temp = temp->next;
        }
        return TV(0);
    };

    // TODO: find y <- ax + y
    // void axpy(int a, LinkedList x){
    // }

};
