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


// template <typename TF, typename TI>
// std::vector<tedge<TF, TI>> adj_edges_vec(std::vector<TF> D, TF thresh) {
//
// }
