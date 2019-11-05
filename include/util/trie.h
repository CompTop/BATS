#include <unordered_map>
#include <vector>

// there is a boost::trie implementation, but it is based on std::map, not std::unordered_map
// template over alphabet type A, storage type T
template <typename A, typename T>
class SparseTrie {
private:
public:
    typedef SparseTrie<A, T> child_type;
    // todo: move this to template type
    typedef std::unordered_map<A, child_type*> child_container;
    T val;
    child_container *children = NULL;

    SparseTrie(const SparseTrie &t) {
        // std::cout << "trie is being copied to " << this << std::endl;
        val = t.val;
        if (t.children != NULL) {
            children = new child_container;
            for (auto& [k, v] : *(t.children)) {
                children->emplace(k, new child_type(*v));
            }
        }
    }

    SparseTrie(T v) : val(v), children(NULL) {}

    SparseTrie() : children(NULL) {};

    ~SparseTrie() {
        // std::cout << "trie is being deleted at " << this << std::endl;
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

    template <typename ITT>
    T& get(ITT stptr, ITT endptr) {
        SparseTrie* current = this;
        while (stptr < endptr) {
            current = current->children->at(*stptr++);
            //stptr++;
        }
        return current->val;
    }

    inline T& operator[](std::vector<T> &k) {
        return get(k.data(), k.data() + k.size());
    }

    // get value with default return
    template <typename ITT>
    T get(ITT stptr, ITT endptr, const T &def_ret) {
        SparseTrie* current = this;
        while (stptr < endptr) {
            if (current->children == NULL || current->children->count(*stptr) == 0) {
                // return the default value
                return def_ret;
            }
            current = current->children->at(*stptr++);
            //stptr++;
        }
        return current->val;
    }

    inline T get(const std::vector<T> &k, const T &def_ret) {
        return get(k.cbegin(), k.cend(), def_ret);
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
