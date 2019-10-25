#include <assert.h>
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
                mat[j*m+i]=0;
            }
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
    
    inline F& operator()(int i, int j) {
        return mat[j*m+i];
    }
    void free(){
        delete mat;
    }
    
    void add_col_to(int i, int j){
        
    }
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
                prod.mat[i*prod.m+j] += m1.mat[i*m1.m+k]*m2.mat[k*m2.m+j];
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
        for (i = j; i < n; i++) {
            Lret(idx_map[i], idx_map[j]) = Lmat(i,j);
        }
    }

    return Lret;

}

/*
solve Aret = Lmat \ Amat
*/
template <typename F>
A<Dense<F>> l_solve(L<Dense<F>> Lmat, A<Dense<F>> Amat) {

    size_t m = Amat.nrow();
    size_t n = Amat.ncol();

    A<Dense<F>> Aret(m,n);
    return Aret;

}
