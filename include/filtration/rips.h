// Struct for holding edge with a value
template <typename TF, typename TI>
struct tedge {
  TF v; // value
  TI s; // source
  TI t; // target
};

template <typename TF, typename TI>
tedge<TF, TI> make_edge(TI s, TI t, TF v) {
  return tedge<TF, TI> {v, s, t};
}

template <typename TF, typename TI>
inline bool operator<(const tedge<TF, TI> &a, const tedge<TF, TI> &b) {
  return a.v < b.v;
}

template <typename TF, typename TI>
std::ostream& operator<<( std::ostream& os, tedge<TF, TI> &x) {
  os << "(" << x.s << ", " << x.t << "): " << x.v;
  return os;
}

// construct edges for rips complex
template <typename T>
void rips_edges(std::vector<T> &x, std::vector<size_t> &edges, std::vector<T> &t) {
    edges.clear();
    size_t nedges = x.size() * (x.size() - 1) / 2;
    edges.reserve(2*nedges);
    t.clear();
    t.reserve(nedges);
    for (size_t i = 0; i < x.size(); i++) {
        for (size_t j = 0; j < i; j++) {
            edges.push_back(j);
            edges.push_back(i);
            t.push_back(std::abs(x[i] - x[j]));
        }
    }
    return;
}

// template <typename TF, typename TI>
// std::vector<tedge<TF, TI>> adj_edges_vec(std::vector<TF> D, TF thresh) {
//
// }
