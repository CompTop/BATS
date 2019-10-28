#include <vector>
#include <iostream>

// Multigraph type
// template over node data, edge data
// inner types for nodes and edges
// multigraph is not responsible for data - always pointer to data

template <typename NT, typename ET>
class MultiGraph {
public:
    struct Node;
    struct Edge;

    struct Node {
        NT* data;
        std::vector<Edge*> input;
        std::vector<Edge*> output;

        Node(NT& x) : data(&x) {};
        Node(NT* x) : data( x) {};

        //inline void add_input(ET in) { inputs.emplace_back(&in); }
        inline void add_input(Edge* in) { input.emplace_back(in); }
        //inline void add_output(ET out) { outputs.emplace_back(&out); }
        inline void add_output(Edge* out) { output.emplace_back(out); }

        inline NT* get_data() { return data; }

        void print() {
            std::cout << '[' << this << "] Node : " << std::endl;
            std::cout << "  [" << data << "] data : " << *data << std::endl;
            std::cout << "  " << input.size() << " inputs" << std::endl;
            std::cout << "  " << output.size() << " outputs" << std::endl;
        }
    };

    struct Edge {
        ET* data;
        Node* src;
        Node* targ;

        Edge(Node& sin, Node& tin, ET& x) : src(&sin), targ(&tin), data(&x) {
            // add edge to src outputs
            sin.add_output(this);
            // add edge to targ inputs
            tin.add_input(this);
        }

        inline ET* get_data() { return data; }

        void print() {
            std::cout << '[' << this << "] Edge : " << std::endl;
            std::cout << "  [" << src << "] source" << std::endl;
            std::cout << "  [" << targ << "] target" << std::endl;
            std::cout << "  [" << data << "] data " << *data << std::endl;
        }
    };

private:
    std::vector<Node> node;
    std::vector<Edge> edge;

public:

    MultiGraph() {}

    Node& add_node(NT &a) {
        node.emplace_back(Node(&a));
        return node.back();
    }

    Edge& add_edge(Node &a, Node &b, ET &data) {
        edge.emplace_back(Edge(a, b, data));
        return edge.back();
    }
};
