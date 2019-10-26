#include <assert.h>
#include <unordered_map>
#include <vector>
#include <util/common.h>
#include <iostream>
#include <algorithm>

// dense vector with entries T
template<typename T>
struct ColumnView{
    T* start;
    T* end;

    ColumnView(T* start, T* end) : start(start), end(end) {};

    // self += a * x
    void axpy(const T a, const ColumnView x){
        T* ptr = start;
        T* xptr = x.start;
        while (ptr < end) {
            *ptr += a * (*xptr);
            ++ptr;
            ++xptr;
        }
        return;
    }

    inline size_t size() const { return end - start; }

    inline T& operator[](size_t i) {return *(start + i); }
}


template<typename F>
struct Dense{};

template<typename F>
struct A<Dense<F>>{
    using DI = Dense<F>;
    F* mat;
    size_t m,n;

    A<DI>(size_t mm, size_t nn, F* mat) : m(mm), n(nn), mat(mat) {}

    A<DI>(size_t mm, size_t nn) : m(mm), n(nn) {
        mat = new F[m*n];
        // std::fill
        std::fill(mat, mat + m*n, F(0));
    }

    A<DI> copy() {
        F* newmat = new F[m*n];
        std::copy(mat, mat + m*n, newmat);

        return A<DI>(m, n, newmat);
    }

    inline size_t nrow() const { return m; }
    inline size_t ncol() const { return n; }

    void print(){
        for( size_t i=0; i<m; i++){
            for( size_t j=0; j<n; j++){
                std::cout<<(mat[j*m+i])<<" ";
            }
            std::cout<<"\n";
        }
    }

    // return a column view of column j
    inline ColumnView<F> operator[](size_t j) {
        return ColumnView<F>(mat + m*j, mat + m*(j+1));
    }

    inline F& operator()(int i, int j) {
        return mat[j*m+i];
    }
    void free(){
        delete mat;
    }

    void add_col_to(size_t i, size_t j, F a);

	void add_row_to(size_t i, size_t j, F a);

    void swap_cols(size_t i, size_t j);
};

template<typename F>
struct T<Dense<F>>:A<Dense<F>>{
      T<Dense<F>>(size_t mm, size_t nn) : A<Dense<F>> (mm,nn) {}
      T<Dense<F>>(size_t mm, size_t nn, F* mat) : A<Dense<F>> (mm,nn,mat) {}
};

template<typename F>
struct L<Dense<F>>:T<Dense<F>> {
      L<Dense<F>>(size_t mm, size_t nn) : T<Dense<F>> (mm,nn) {}
      L<Dense<F>>(size_t mm, size_t nn, F* mat) : T<Dense<F>> (mm,nn,mat) {}
};

template<typename F>
struct EL<Dense<F>> :L<Dense<F>>{
      EL<Dense<F>>(size_t mm, size_t nn) : L<Dense<F>> (mm,nn) {}
      EL<Dense<F>>(size_t mm, size_t nn, F* mat) : L<Dense<F>> (mm,nn,mat) {}

};

template< typename F>
A<Dense<F>> matmul(A<Dense<F>> m1, A<Dense<F>> m2){
    assert( m1.n==m2.m );
    A<Dense<F>> prod(m1.m,m2.n);
    for(size_t i=0;i<m1.m;i++)
        for(size_t j=0;j<m2.n;j++)
            for(size_t k=0;k<m1.n;k++){
                prod(i,j) += m1(i,k)*m2(k,j);
            }
    return prod;
}

/*
Lret * ELmat = ELmat *Lmat
*/
template <typename F>
L<Dense<F>> el_commute(EL<Dense<F>> ELmat, L<Dense<F>> Lmat) {
    size_t m = ELmat.nrow();
    size_t n = ELmat.ncol();
    assert(n == Lmat.ncol());
    assert(Lmat.ncol() == Lmat.nrow());
    L<Dense<F>> Lret = L<Dense<F>>(m, m);

    // step 1 is to make an index map from ELmat
    std::vector<size_t> idx_map(n);
    size_t i = 0; // pointer to current row
    size_t j = 0; // pointer to column
    // loop over columns of ELmat
    for (j = 0; j < n && i < m; j++) {
        // increment i until we find a non-zero, or run out of column
        while (ELmat(i,j) == 0 && i < m) {
            i++;
        }
        // column is all zero
        if (i >= m) {
            idx_map[j] = bats::NO_IND;
            continue;
        }
        // i is the non-zer index
        idx_map[j] = i;
    }

    // step 2 is to form Lret
    for (j = 0; j < m; j++) {
        Lret(j,j) = F(1); // default unit-diagonal
    }
    // perform actual commutation
    for (j = 0; j < n; j++) {
        if (idx_map[j] == bats::NO_IND) { break; }
        for (i = j; i < n; i++) {
            if (idx_map[i] == bats::NO_IND) { break; }
            Lret(idx_map[i], idx_map[j]) = Lmat(i,j);
        }
    }

    return Lret;

}

// solve Lmat x = y in-place
template <typename F>
void l_solve(L<Dense<F>> Lmat, ColumnView<F> y) {
    size_t n = y.size();
    for (size_t k = 0; k < n; k++) {
        y[k] /= Lmat(k,k);
        // apply forward-looking update
        for (size_t j = k+1; j < n; j++) {
            y[j] -= Lmat(j, k) * y[k];
        }
    }
    return;
}

/*
solve Aret = Lmat \ Amat
*/
template <typename F>
A<Dense<F>> l_solve(L<Dense<F>> Lmat, A<Dense<F>> Amat) {

    size_t m = Amat.nrow();
    size_t n = Amat.ncol();

    A<Dense<F>> Aret = Amat.copy();

    // apply lower triangular solve to each column
    for (size_t j = 0; j < n; j++) {
        l_solve(Lmat, Aret[j]);
    }

    return Aret;

}
