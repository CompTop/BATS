#pragma once
#include <iostream>
using namespace std;

// This is special linked-list used for column of Vineyard matrix, where
// row index is shared by all lists
template <typename TV>
class Node { 
public:
    TV data;
    const size_t* & index_ptr; // row index
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


    Node(TV data, const size_t* & ptr) : data(data), index_ptr(ptr), next(nullptr) {
        // cout << "Constructed node: " << data << "\n";
    }
};

template <typename TV> // type of node
class LinkedList{
private:
    Node<TV>* head = nullptr;
public:
    using node_type = Node<TV>;
    using rptr_type = size_t *;
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
            const rptr_type & r_ptr = *(row_it+ it->ind);
            // const auto& r_ptr = *(row_it+ it->ind); // const reference to rvalue
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
                const auto& r_ptr = *(row_it + (*itr)); // const reference to rvalue
                insert(*itv, r_ptr);
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
    // template<typename TPtr>
    // void insert(TV val, TPtr* ptr){
    void insert(const TV &val, const size_t* & ptr){
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

    // TODO: find y <- ax + y still onging!
    // assume x shares the same vector of row index pointers
    template <typename TV_new>
    void axpy(TV_new _a, LinkedList<TV> x){
        LinkedList<TV> temp_ll;
        TV a = TV(_a);
        auto yptr = this->head;
        auto xptr = x.head;
        while (xptr != nullptr and yptr != nullptr){
            // compare memory addresses of row pointers instead of their values
            if (&(xptr->index_ptr) < &(yptr->index_ptr)){
                temp_ll.insert(a * xptr->data, xptr->index_ptr);
                xptr = xptr->next; 
            }else if (&(xptr->index_ptr) > &(yptr->index_ptr)){ 
                temp_ll.insert(yptr->data, yptr->index_ptr);
                yptr = yptr->next;
            }else{
                temp_ll.insert(yptr->data + a *  xptr->data, yptr->index_ptr);
                yptr = yptr->next;
                xptr = xptr->next; 
            }
        }
        this->clear();
        this->head = temp_ll.head;
    }

};
