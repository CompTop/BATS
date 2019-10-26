#include <vector>
#include <iostream>

// node type templated over data
// store incoming and outgoing edges

// edge type templated over data
// store source, target, and data

template <typename T>
struct Node {
private:
    T* data;
    std::vector< void* > inputs;
    std::vector< void* > outputs;

public:
    Node(T &x) : data(&x) {};

    //inline void add_input(ET in) { inputs.emplace_back(&in); }
    template <typename ET>
    inline void add_input(ET* in) {
        std::cout << "adding input" << std::endl;
        std::cout << inputs.size();
        inputs.emplace_back((void*) in);
        std::cout << inputs.size() << std::endl;
    }
    //inline void add_output(ET out) { outputs.emplace_back(&out); }
    template <typename ET>
    inline void add_output(ET* out) { outputs.emplace_back((void*) out); }

    inline T* get_data() { return data; }

    void print() {
        std::cout << '[' << this << "] Node : " << std::endl;
        std::cout << "  [" << data << "] data : " << *data << std::endl;
        std::cout << "  " << inputs.size() << " inputs" << std::endl;
        std::cout << "  " << outputs.size() << " outputs" << std::endl;

    }
};

// template over T - data type
template <typename T>
struct Edge {
private:
    void* src;
    void* targ;
    T* data;

public:

    template <typename NT>
    Edge(Node<NT> &sin, Node<NT> &tin, T& x) : src((void*) &sin), targ((void*) &tin), data(&x) {
        // add edge to src outputs
        sin.add_output(this);
        // add edge to targ inputs
        tin.add_input(this);
    };

    inline T* get_data() { return data; }



    void print() {
        std::cout << '[' << this << "] Edge : " << std::endl;
        std::cout << "  [" << src << "] source" << std::endl;
        std::cout << "  [" << targ << "] target" << std::endl;
        std::cout << "  [" << data << "] data " << *data << std::endl;
    }

};
