#include <unordered_map>
#include <vector>

// there is a boost::trie implementation, but it is based on std::map, not std::unordered_map
template <typename A, typename T>
class SparseTrie {
private:
public:
    typedef SparseTrie<A, T> child_type;
    // todo: move this to template type
    typedef std::unordered_map<A, child_type*> child_container;
    T val;
    child_container *children = NULL;

    SparseTrie(T v) : val(v), children(NULL) {}



    SparseTrie() : children(NULL) {};

    ~SparseTrie() {
//         std::cout << "trie is being deleted at " << this << std::endl;
        if (children != NULL) {
//             std::cout << "children at " << children << std::endl;
            // delete each child
            for (auto kvpair : *children) {
//                 std::cout << "deleting child at " << kvpair.second << std::endl;
//                 std::cout << "key " << kvpair.first << std::endl;
//                 std::cout << "value " << kvpair.second->val << std::endl;
                delete kvpair.second;
//                 std::cout << "deleted." << std::endl;
            }
//             std::cout << "deleting map at " << children << std::endl;
            delete children;
        }
    }


    void insert(A* stptr, const A* endptr, const T &v) {
        SparseTrie* current = this;
        while (stptr < endptr) {
            if (current->children == NULL) { current->children = new child_container; }
            if (current->children->count(*stptr) == 0) { current->children->emplace(*stptr, new child_type); }
            current = current->children->at(*stptr++);
        }
        current->val = v;
        return;
    }


    inline void emplace(std::vector<A> &k, T &&v) {
        return insert(k.data(), k.data() + k.size(), v);
    }

    inline void emplace(std::vector<A> &k, T &v) {
        return insert(k.data(), k.data() + k.size(), v);
    }


    T get(A* stptr, A* endptr) {
        SparseTrie* current = this;
        while (stptr < endptr) {
            current = current->children->at(*stptr++);
            //stptr++;
        }
        return current->val;
    }

    inline T operator[](std::vector<T> &k) {
        return get(k.data(), k.data() + k.size());
    }

    size_t count(A* stptr, A* endptr) {
        SparseTrie* current = this;
        while (stptr < endptr) {
            if (current->children == NULL || current->children->count(*stptr) == 0) {
                return 0;
            }
            current = current->children->at(*stptr++);
            //stptr++;
        }
        return 1;
    }

    inline size_t count(std::vector<T> &k) {
        return count(k.data(), k.data() + k.size());
    }

};
