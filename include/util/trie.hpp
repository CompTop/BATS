#pragma once

// #include <unordered_map>
#include <vector>
// #include "tsl/hopscotch_map.h"
// #include "tsl/robin_map.h"

namespace bats {

// there is a boost::trie implementation, but it is based on std::map, not std::unordered_map
// template over alphabet type A, storage type T
template <typename A, typename T>
class SparseTrie {
private:
public:
    typedef SparseTrie<A, T> child_type;
    // todo: move this to template type
    typedef std::unordered_map<A, child_type*> child_container;
	// typedef tsl::hopscotch_map<A, child_type*> child_container;
    T val;
    child_container *children = nullptr;

	// Copy constructor - deep copy
    SparseTrie(const SparseTrie &t) {
        //std::cout << "trie is being copied to " << this << std::endl;
        val = t.val;
        if (t.children != nullptr) {
            children = new child_container;
            for (auto& [k, v] : *(t.children)) {
                children->emplace(k, new child_type(*v));
            }
        }
    }

    // Move constructor
	// Transfer ownership of children dictionary
	SparseTrie(SparseTrie&& t)
		: val(t.val), children(t.children)
	{
		t.children = nullptr; //
	}

    SparseTrie(T v) : val(v), children(nullptr) {}

    SparseTrie() : val(T(0)), children(nullptr) {};

    ~SparseTrie() {
        if (children != nullptr) {
			// delete each child
            for (auto kvpair : *children) {
                delete kvpair.second;
            }
            delete children;
        }
    }

	// Copy assignment
	// performs deep copy
	SparseTrie& operator=(const SparseTrie& t) {

		if (&t == this) {return *this;}

		// release children
		if (children != nullptr) { delete children; }

		val = t.val;
        if (t.children != nullptr) {
            children = new child_container;
            for (auto& [k, v] : *(t.children)) {
                children->emplace(k, new child_type(*v));
            }
        }

		return *this;
	}

	// Move assignment
	// Transfer ownership of children
	SparseTrie& operator=(SparseTrie&& t)
	{
		// Self-assignment detection
		if (&t == this)
			return *this;

		// release children
		if (children != nullptr) { delete children; }

		val = t.val;
		// Transfer ownership of children
		children = t.children;
		t.children = nullptr;

		return *this;
	}

    template <typename ITT>
    void insert(ITT stptr, const ITT endptr, const T &v) {
        SparseTrie* current = this;
        while (stptr < endptr) {
            if (current->children == nullptr) { current->children = new child_container; }
            if (current->children->count(*stptr) == 0) { current->children->emplace(*stptr, new child_type); }
            current = current->children->at(*stptr++);
        }
        current->val = v;
        return;
    }


    inline void emplace(const std::vector<A> &k, T &&v) {
        return insert(k.cbegin(), k.cend(), v);
    }

    inline void emplace(const std::vector<A> &k, const T &v) {
        return insert(k.cbegin(), k.cend(), v);
    }

    template <typename ITT>
    T& get(ITT stptr, const ITT endptr) {
        SparseTrie* current = this;
        while (stptr < endptr) {
            current = current->children->at(*stptr++);
            //stptr++;
        }
        return current->val;
    }

    inline T& operator[](const std::vector<T> &k) {
        return get(k.cbegin(), k.cend());
    }

    // get value with default return
    template <typename ITT>
    T get(ITT stptr, const ITT endptr, const T &def_ret) {
        SparseTrie* current = this;
        while (stptr < endptr) {
            if (current->children == nullptr || current->children->count(*stptr) == 0) {
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

    // get value with default return
    // maintains const correctness
    template <typename ITT>
    T get(ITT stptr, const ITT endptr, const T &def_ret) const {
        if (stptr < endptr) {
            if (children == nullptr || children->count(*stptr) == 0) {
                return def_ret;
            }
            return children->at(*stptr)->get(++stptr, endptr, def_ret);
        }
        return val;
    }

    inline T get(const std::vector<T> &k, const T &def_ret) const {
        return get(k.cbegin(), k.cend(), def_ret);
    }


    template <typename ITT>
    size_t count(ITT stptr, const ITT endptr) {
        SparseTrie* current = this;
        while (stptr < endptr) {
            if (current->children == nullptr || current->children->count(*stptr) == 0) {
                return 0;
            }
            current = current->children->at(*stptr++);
            //stptr++;
        }
        return 1;
    }

    inline size_t count(const std::vector<T> &k) {
        return count(k.cbegin(), k.cend());
    }

};

} // namespace bats
