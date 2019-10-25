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

