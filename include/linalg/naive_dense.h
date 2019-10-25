#include <assert>
#include <unordered_map>
#include <vector>
#include <util/common.h>
#include<iostream>

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
        for( size_t i=0; i<m; i++)
            for( size_t j=0; j<n; j++){
                mat[i*m+j]=0;
            }
    }


    inline size_t nrow() const { return m; }
    inline size_t ncol() const { return n; }

    inline F operator()(int i, int j) {
        return mat[i*m+j];
    }
    void free(){
        delete mat;
    }
    
    void print(){
        for( size_t i=0; i<m; i++){
            for( size_t j=0; j<n; j++){
                std::cout<<(mat[i*m+j])<<" ";
            }
            std::cout<<"\n";
        }
    }
};

template< typename F>
A<Dense<F>> matmul(A<Dense<F>> m1, A<Dense<F>> m2){
    assert( m1.n==m2.m );
    A<Dense<F>> prod(m1.m,m2.n);
    for(size_t i=0;i<m1.m;i++)
        for(size_t j=0;j<m2.n;j++)
            for(size_t k=0;k<m1.n;k++){
                prod.mat[i*prod.m+j] += m1.mat[i*m1.m+k]*m2.mat[k*m2.m+j];
            }
    return prod;
}

/*
L' E_L = E_L L
*/
template <typename F>
L<Dense<F>> el_commute(EL<Dense<F>> ELmat, L<Dense<F>> Lmat) {
    size_t m = ELmat.nrow();
    size_t n = ELmat.ncol();
    assert(n == Lmat.ncol());
    assert(Lmat.ncol() == Lmat.nrow());
    L<Dense<F>> Lret = L<Dense<F>>(m, m);

    // step 1 is to make an index map from ELmat
    std::vector<int> idx_map(n);
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
        for (i = j; i < n; i++) {
            Lret(idx_map[i], idx_map[j]) = Lmat(i,j);
        }
    }

    return Lret;

}
