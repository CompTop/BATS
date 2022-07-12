#pragma once
#include <memory>
#include <iostream>
#include <bats.hpp>
using namespace std;

// This is special linked-list used for column of Vineyard matrix, where
// row index is shared by all lists
template <typename TV>
class Node { 
public:
    TV data;
    const std::shared_ptr<size_t> * index_ptr; // pointer to the address of pointers of row index
    Node* next;
    // Defalt constructor
    Node() : data(0), index_ptr(nullptr), next(nullptr) {
        // cout << "Constructed default node\n";
    }
    // copy constructor 
    Node(const Node &node2) : data(node2.data){
        // cout << "Copy constructor for data " << data <<  "\n";
        index_ptr = node2.index_ptr;
        next = node2.next;
    }

    Node(TV data, const std::shared_ptr<size_t>& sptr) : data(data), index_ptr(&sptr), 
    next(nullptr) {// pass into pointer and store its address
        // cout << "Constructed node: " << data << "\n";
    }
    
    // find row index value of current node 
    inline size_t get_ind_val(){
        return *(get_ind_ptr());
    }

    // find the address of row index pointer of current node 
    // i.e., the row position in row_inds_ptr
    inline auto get_ind_ptr (){
        return index_ptr->get();
    }

};

template <typename TV> // type of node
class LinkedList{
private:
    Node<TV>* head = nullptr;
public:
    using Node_type = Node<TV>;
    using tmp_type = LinkedList<TV>;
    LinkedList() { // constructor
        head = nullptr;
    }

    // constructor from a sparse vector 
    // begin_it: begin iterator of sparse vector
    // end_it: end iterator of sparse vector.
    // row_it: begin iterator of the row indices vector
    template <typename IT1, typename IT2>
    LinkedList(IT1 begin_it, IT1 end_it, IT2 row_it) {
        if (begin_it == end_it){// check empty sparse vector 
            head = nullptr;
        }else{
            auto it = begin_it;
            head = new Node_type(it->val, *(row_it + it->ind));
            Node_type* temp = head;
            it++;
            while (it < end_it) {
                temp->next = new Node_type(it->val, *(row_it + it->ind));
                temp = temp->next;
                it++;
            }
        }
        this->clear_zeros();
	}

    // Constructor from a column of CSC matrix
    // (assume the vector of row index pointer is the default sorted one)
    // rptr: row begin iterator 
    // vptr: value begin iterator 
    // row_it: begin iterator of the row indices vector
    // col_size: the number of (elements) nodes needed in the column
    template <typename IT1, typename IT2, typename IT3>
    LinkedList(IT1 rptr, IT2 vptr, IT3 row_it, size_t col_size) {
        head = nullptr;
        auto itr = rptr;
        auto itv = vptr;
        // deal with head first
        while(itr < rptr + col_size){
            if (TV(*itv) != 0){ // zero node is not needed
                // insert(*itv, *(row_it + (*itr)));
                head = new Node_type(*itv, *(row_it + (*itr)));
                itv++;
                itr++;
                break;
            }
            itv++;
            itr++;
        }
        if (head == nullptr){return;} // no head, i.e., zero column
        Node_type* temp = head;
        // for the rest
        while(itr < rptr + col_size){
            if (TV(*itv) != 0){
                Node_type* new_node = new Node_type(*itv, *(row_it + (*itr)));
                temp->next = new_node;
                temp = temp->next;
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
            Node_type* temp2 = ll2.head;
            Node_type* next2 = temp2->next;
            head = new Node_type(temp2->data, *(temp2->index_ptr));
            // head = new Node_type(temp2);
            Node_type* temp = head; // for new linkedlist
            while (temp2->next != nullptr){
                Node_type* new_node = new Node_type(next2->data, *(next2->index_ptr));
                // Node_type* new_node = new Node_type(next2);
                temp->next = new_node;
                temp = temp->next;
                temp2 = next2;
                next2 = temp2->next;
            }
        }
    }

    // copy assigenment operator 
    LinkedList& operator=(const LinkedList& ll2)
    {
        // std::cout << "Call copy assignment operator \n";
        if (this != &ll2){
            // Free the existing resource.
            this->~LinkedList();
            // same as copy constructor as above
            if (ll2.head == nullptr) {
                head = nullptr;
            }else{
                Node_type* temp2 = ll2.head;
                Node_type* next2 = temp2->next;
                head = new Node_type(temp2->data, *(temp2->index_ptr));
                // head = new Node_type(temp2);
                Node_type* temp = head; // for new linkedlist
                while (temp2->next != nullptr){
                    Node_type* new_node = new Node_type(next2->data, *(next2->index_ptr));
                    // Node_type* new_node = new Node_type(next2);
                    temp->next = new_node;
                    temp = temp->next;
                    temp2 = next2;
                    next2 = temp2->next;
                }
            }
        }
        return *this;
    }

    // move constructor (constructing new object, where head now is null)
    LinkedList(LinkedList&& ll2) noexcept{
        std::cout << "move constructor is called." << std::endl;
        if (ll2.head == nullptr) 
            return;
        else{
            // head = ll2.head; // steal temp obj
            // ll2.head = nullptr; // nullize temp obj
            head = std::move(ll2.head);
            ll2.head = nullptr;
        }
    }
    
    // move assignment operator (replace head by other object, where head might not be null)
    LinkedList& operator=(LinkedList&& other) noexcept{
        // std::cout << "move assignment operator is called." << std::endl;
        if (this != &other){
            this->~LinkedList(); // delete current object first
            head = other.head; // steal temp obj
            other.head = nullptr; // nullize temp obj
        }
        return *this;
    }

    // destructor
    ~LinkedList() {
        // std::cout << "Call the destructor of LinkedList" << std::endl;
        clear();
    }

    // Return tail node 
    Node_type* tail_node(){
        if (head == nullptr) {
            return head;
        }else {
            Node_type* temp = head; // head is not nullptr
            while (temp->next != nullptr) { 
                temp = temp->next; // go to end of list
            }
            return temp;
        }
    }

	// number of nonzeros
	size_t nnz() const {
		size_t ct = 0;
		if (head == nullptr) {
            return ct;
        }else{
            ct++;
            Node_type* temp = head; // head is not nullptr
            while (temp->next != nullptr) { 
                ct++;
                temp = temp->next; // go to end of list
            }
        }
		return ct;
	}

    // Return the last non-zero node (tail node)
    nzpair<size_t, TV> lastnz(){
        if (head == nullptr) {
            return nzpair<size_t, TV>();
        }else {
            TV last_data = head->data;
            size_t last_ind = head->get_ind_val();
            Node_type* temp = head; // head is not nullptr
            while (temp->next != nullptr) {
                temp = temp->next;
                if (temp->get_ind_val() > last_ind){
                    last_ind = temp->get_ind_val();
                    last_data = temp->data;
                }
            }
            return nzpair<size_t, TV>(last_ind, last_data);
        }
    }

    // insert node at the end of linked list 
    void insert(const TV &val, const std::shared_ptr<size_t> & sptr){
        if(val == TV(0)){return;}

        Node_type* newnode = new Node_type(val, sptr);
        if (head == nullptr) {
            head = newnode;
        }
        else {
            Node_type* temp = head; // head is not nullptr
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
            Node_type* temp = head;
            std::cout << "(index, value): ";
            while (temp != nullptr) {
                cout << "("<< temp->get_ind_val() << ", " << temp->data << ") ";
                temp = temp->next;
            }
            cout << endl;
        }
    }

    void print(){
        if (head == nullptr) {
            cout << "List is empty!" << endl;
        }
        else {
            Node_type* temp = head;
            std::cout << "(index, value): ";
            while (temp != nullptr) {
                cout << "("<< temp->get_ind_val() << ", " << temp->data << ") ";
                temp = temp->next;
            }
            cout << endl;
        }
    }

    // clear zeros inside
    void clear_zeros(){
        if (head == nullptr){return;}

        Node_type* prev_nd; // previous node
        Node_type* next_nd; // next node 
        Node_type* temp_nd = head; // current node
        // find the first non-zero as its new head 
        while(temp_nd != nullptr){
            if (temp_nd->data == TV(0)){ 
                next_nd = temp_nd->next;
                delete temp_nd;
                head = next_nd;
                temp_nd = next_nd;
            }else{
                prev_nd = temp_nd;
                temp_nd = temp_nd->next;
                break;
            }
        }

        while (temp_nd != nullptr) {
            // if temp is zero, no need to modify previous node
            if (temp_nd->data == TV(0)){ 
                next_nd = temp_nd->next;
                prev_nd->next = next_nd;
                delete temp_nd;
                temp_nd = next_nd;
            }else{
                prev_nd = temp_nd;
                temp_nd = temp_nd->next;
            }
        }
    };

    void clear(){
        Node_type* next; // next node 
        Node_type* temp = head; // current node
        while (temp != nullptr) {
            next = temp->next;
            delete temp;
            temp = next;
        }
    };

    TV getval(size_t i) const{
        Node_type* temp = head;
        // iterate over all nodes, might be inefficient
        while (temp != nullptr) { 
            if (temp->get_ind_val() == i){
                return temp->data;
            }
            temp = temp->next;
        }
        return TV(0);
    };

    // find y <- ax + y by modifying y in-place
    // assume x shares the same vector of row index pointers
    template <typename TV_new>
    void axpy_inplace(TV_new _a, LinkedList<TV> x){
        TV a = TV(_a);
        auto yptr = this->head;
        auto xptr = x.head;

        if (a == TV(0) or xptr == nullptr){return;}
        
        if (yptr == nullptr){
            // call the copy assignment operator
            *this = x;
            auto temp = this->head;
            while (temp != nullptr){
                temp->data = a * (temp->data); 
                temp = temp->next;
            }
        }else{ // none of x and y are nullptr
            // Deal with head first to make sure yptr always has a node above it
            auto x_address = xptr->get_ind_ptr(); // address in row indices vector
            auto y_address = yptr->get_ind_ptr();
            // compare their addresses of row index pointers (since they initialized in vector)
            if (x_address > y_address){
                yptr = yptr->next;
            }else if (x_address < y_address){
                head = new Node_type(a * xptr->data, *(xptr->index_ptr));
                head->next = yptr;
                xptr = xptr->next;
            }else{
                yptr->data += a * xptr->data; 
                yptr = yptr->next;
                xptr = xptr->next;
            }
            auto yptr_pre = head; // previous node of yptr
            // Strat look at the rest nodes 
            while (xptr != nullptr and yptr != nullptr){ 
                x_address = xptr->get_ind_ptr();
                y_address = yptr->get_ind_ptr();
                if (x_address > y_address){
                    // keep y the same and only increment yptr
                    yptr_pre = yptr;
                    yptr = yptr->next;
                }else if (x_address < y_address){
                    // insert new node to y and only increment xptr 
                    Node_type* new_node = new Node_type(a * xptr->data, *(xptr->index_ptr));
                    yptr_pre->next = new_node;
                    new_node->next = yptr;
                    xptr = xptr->next;
                }else{
                    // add a*x to y and increment both
                    yptr->data += a * xptr->data; 
                    yptr_pre = yptr;
                    yptr = yptr->next;
                    xptr = xptr->next;
                }
            } 

            // Now at least one of the x y pointers are nullptr
            while (xptr != nullptr){ // yptr == nullptr is true 
                // which means in the last iteration of the above while loop, 
                // yptr = yptr->next; has been called (if case 1 and 3).
                // For case 1: yptr_pre < xptr, only need to copy the rest of x to y
                // For case 3: yptr_pre = xptr_pre, in other words xptr > yptr_pre, back to the case above
                assert((xptr->index_ptr)->get() > (yptr_pre->index_ptr)->get());
                Node_type* new_node = new Node_type(a * xptr->data, *(xptr->index_ptr));
                yptr_pre->next = new_node;
                yptr_pre = yptr_pre->next;
                xptr = xptr->next;
            }
            // if xptr == nullptr is true, no need to merge
        }
        this->clear_zeros();
    }

    // find y <- ax + y by using temporary linked list
    // assume x shares the same vector of row index pointers
    template <typename TV_new>
    void axpy(TV_new _a, LinkedList<TV> x){
        TV a = TV(_a);
        auto yptr = this->head;
        auto xptr = x.head;

        if (a == TV(0) or xptr == nullptr){return;}
        
        if (yptr == nullptr){
            // call the copy assignment operator
            *this = x;
            auto temp = this->head;
            while (temp != nullptr){
                temp->data = a * (temp->data); 
                temp = temp->next;
            }
            this->clear_zeros();
            return;
        }
        
        // Now, none of x and y are nullptr
        LinkedList<TV> templl;
        // find head first
        auto x_address = xptr->get_ind_ptr(); // address in the shared row indices vector
        auto y_address = yptr->get_ind_ptr();
        if (x_address > y_address){ 
            templl.head = new Node_type(yptr->data, *(yptr->index_ptr));
            yptr = yptr->next;
        }else if (x_address < y_address){
            templl.head = new Node_type(a * xptr->data, *(xptr->index_ptr));
            xptr = xptr->next;
        }else{
            templl.head = new Node_type(yptr->data + a * xptr->data, *(yptr->index_ptr));
            xptr = xptr->next;
            yptr = yptr->next;
        }
        
        
        Node_type* temp_ptr = templl.head;
        while (xptr != nullptr and yptr != nullptr){ 
            x_address = xptr->get_ind_ptr(); // address in row indices vector
            y_address = yptr->get_ind_ptr();

            if (x_address > y_address){ 
                // only deal with y
                temp_ptr->next = new Node_type(yptr->data, *(yptr->index_ptr));
                yptr = yptr->next;
            }else if (x_address < y_address){
                temp_ptr->next = new Node_type(a * xptr->data, *(xptr->index_ptr));
                xptr = xptr->next;
            }else{
                temp_ptr->next = new Node_type(yptr->data + a * xptr->data, *(yptr->index_ptr));
                xptr = xptr->next;
                yptr = yptr->next;
            }
            temp_ptr = temp_ptr->next;
        }
        // Now at least one of the x y pointers are nullptr
        while (xptr != nullptr){ // yptr == nullptr is true 
            temp_ptr->next = new Node_type(a * xptr->data, *(xptr->index_ptr));
            xptr = xptr->next;
            temp_ptr = temp_ptr->next;
        }

        while (yptr != nullptr){ // xptr == nullptr is true 
            temp_ptr->next = new Node_type(yptr->data, *(yptr->index_ptr));
            yptr = yptr->next;
            temp_ptr = temp_ptr->next;
        }
        
        // std::cout << "\nshould call move assignment operator now" << std::endl;
        *this = std::move(templl);

        this->clear_zeros();
    }

    // find y <- ax + y by passing temporary linked list
    // assume x shares the same vector of row index pointers
    template <typename TV_new>
    void axpy(TV_new _a, LinkedList<TV> x, LinkedList<TV> templl){
        TV a = TV(_a);
        auto yptr = this->head;
        auto xptr = x.head;

        if (a == TV(0) or xptr == nullptr){return;}
        
        if (yptr == nullptr){
            // call the copy assignment operator
            *this = x;
            auto temp = this->head;
            while (temp != nullptr){
                temp->data = a * (temp->data); 
                temp = temp->next;
            }
            this->clear_zeros();
            return;
        }
        
        // Now, none of x and y are nullptr
        templl.clear();
        // find head first
        auto x_address = xptr->get_ind_ptr(); // address in the shared row indices vector
        auto y_address = yptr->get_ind_ptr();
        if (x_address > y_address){ 
            templl.head = new Node_type(yptr->data, *(yptr->index_ptr));
            yptr = yptr->next;
        }else if (x_address < y_address){
            templl.head = new Node_type(a * xptr->data, *(xptr->index_ptr));
            xptr = xptr->next;
        }else{
            templl.head = new Node_type(yptr->data + a * xptr->data, *(yptr->index_ptr));
            xptr = xptr->next;
            yptr = yptr->next;
        }
        
        
        Node_type* temp_ptr = templl.head;
        while (xptr != nullptr and yptr != nullptr){ 
            x_address = xptr->get_ind_ptr(); // address in row indices vector
            y_address = yptr->get_ind_ptr();

            if (x_address > y_address){ 
                // only deal with y
                temp_ptr->next = new Node_type(yptr->data, *(yptr->index_ptr));
                yptr = yptr->next;
            }else if (x_address < y_address){
                temp_ptr->next = new Node_type(a * xptr->data, *(xptr->index_ptr));
                xptr = xptr->next;
            }else{
                temp_ptr->next = new Node_type(yptr->data + a * xptr->data, *(yptr->index_ptr));
                xptr = xptr->next;
                yptr = yptr->next;
            }
            temp_ptr = temp_ptr->next;
        }
        // Now at least one of the x y pointers are nullptr
        while (xptr != nullptr){ // yptr == nullptr is true 
            temp_ptr->next = new Node_type(a * xptr->data, *(xptr->index_ptr));
            xptr = xptr->next;
            temp_ptr = temp_ptr->next;
        }

        while (yptr != nullptr){ // xptr == nullptr is true 
            temp_ptr->next = new Node_type(yptr->data, *(yptr->index_ptr));
            yptr = yptr->next;
            temp_ptr = temp_ptr->next;
        }
        
        // std::cout << "\nshould call move assignment operator now" << std::endl;
        *this = std::move(templl);

        this->clear_zeros();
    }
};
