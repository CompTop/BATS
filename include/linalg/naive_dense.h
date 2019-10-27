#include <assert.h>
#include <unordered_map>
#include <vector>
#include <util/common.h>
#include <iostream>
#include <algorithm>
#include "matrix_interface.h"

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

	void operator=(const ColumnView x){
        T* ptr = start;
        T* xptr = x.start;
        while (ptr < end) {
            *ptr = (*xptr);
            ++ptr;
            ++xptr;
        }
        return;
    }

    inline size_t size() const { return end - start; }

    inline T& operator[](size_t i) {return *(start + i); }
};


template<typename T>
struct RowView{
    T* start;
    T* end;
	size_t stride;

    RowView(T* start, T* end, size_t s) : start(start), end(end), stride(s) {};

    // self += a * x
    void axpy(const T a, const RowView x){
        T* ptr = start;
        T* xptr = x.start;
        while (ptr < end) {
            *ptr += a * (*xptr);
            ptr+=stride;
            xptr+=x.stride;
        }
        return;
    }

	void operator=(const RowView x){
        T* ptr = start;
        T* xptr = x.start;
        while (ptr < end) {
            *ptr =  (*xptr);
            ptr+=stride;
            xptr+=x.stride;
        }
        return;
    }

    //inline size_t size() const { return (end - start)/stride; }

    inline T& operator[](size_t i) {return *(start + stride*i); }
};



template<typename F>
struct Dense{};

template<typename F>
struct A<Dense<F>>{
    using DI = Dense<F>;

    size_t m,n;
    F* mat;

	// null constructor
	A<DI>(){ m=0;n=0;}

	//"copy" constructor
	A<DI>(const A<DI> &m2) {m = m2.m; n = m2.n; mat = m2.mat; }

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
                std::cout<<((*this)(i,j))<<" ";
            }
            std::cout<<"\n";
        }
    }

	void print_arr(){
		std::cout<<"{\n";
        for( size_t i=0; i<m; i++){
            for( size_t j=0; j<n; j++){
                std::cout<<((*this)(j,i))<<",";
            }
            std::cout<<"\n";
        }
		std::cout<<"};";
    }


    // return a column view of column j
    inline ColumnView<F> operator[](size_t j) {
        return ColumnView<F>(mat + m*j, mat + m*(j+1));
    }

	// return a row view of row i
    inline RowView<F> r(size_t i) {
        return RowView<F>(mat + i, mat + m*n + i, m );
    }

    inline F& operator()(int i, int j) {
        return mat[j*m+i];
    }
    void free(){
        delete mat;
    }

	//equality check
	template<typename M>
	bool operator==(M&& other) {
		if( m!=other.m || n!=other.n)
			return false;

        for(size_t i=0;i<m;i++){
    		for(size_t j=0;j<n;j++)
				if( (*this)(i,j) != other(i,j) )
					return false;
		}

		return true;
    }

    void add_col_to(size_t i, size_t j, F a);

	void add_row_to(size_t i, size_t j, F a);

    void swap_cols(size_t i, size_t j);

	void swap_rows(size_t i, size_t j);
};


// Inherit all the constructors
#define INHERIT(T1,T2) \
template<typename F> \
struct T1<Dense<F>>:T2<Dense<F>>{ \
	using Base = T2<Dense<F>>; \
	using Base::Base ; \
}; \

INHERIT(T,A)
INHERIT(L,T)
INHERIT(U,T)
INHERIT(E,A)
INHERIT(EL,E)
INHERIT(EU,E)
INHERIT(P,E)


// naive matmul
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

    //size_t m = Amat.nrow();
    size_t n = Amat.ncol();

    A<Dense<F>> Aret = Amat.copy();

    // apply lower triangular solve to each column
    for (size_t j = 0; j < n; j++) {
        l_solve(Lmat, Aret[j]);
    }

    return Aret;

}


// implementations of naive matrix ops

template<typename F>
void A<Dense<F>>::swap_cols(size_t i, size_t j){
    for(size_t k=0;k<m;k++){
        std::swap((*this)(k,i),(*this)(k,j));
    }
}

template<typename F>
void A<Dense<F>>::swap_rows(size_t i, size_t j){
    for(size_t k=0;k<n;k++){
        std::swap((*this)(i,k),(*this)(j,k));
    }
}

template<typename F>
void A<Dense<F>>::add_col_to(size_t i, size_t j, F a){
    for(size_t k=0;k<m;k++){
       (*this)(k,i) = (*this)(k,i) + a*(*this)(k,j);
    }
}

template<typename F>
void A<Dense<F>>::add_row_to(size_t i, size_t j, F a){
    for(size_t k=0;k<n;k++){
       (*this)(i,k) = (*this)(i,k) + a*(*this)(j,k);
    }
}


// mat creation utils

template<typename M>
void make_diag_ones(M mat){
    for(size_t i=0;i<mat.m;i++)
        mat(i,i)=1;
}
template<typename M>
void fill_zeros(M mat){
    for(size_t i=0;i<mat.m;i++)
        for(size_t j=0;j<mat.n;j++)
            mat(i,j)=0;
}
template<typename M>
void fill_rand(M mat){
    for(size_t i=0;i<mat.m;i++)
        for(size_t j=0;j<mat.n;j++)
            mat(i,j)=rand();
}


// Find pivot in lower right block of pr,pc
template<typename M>
auto find_pivot(M mat,size_t pr, size_t pc){
    size_t prow = pr-1; // start 1 less as its immediately incremented
    size_t pcol = pc;

    // outer loop searches rows
    do{
        prow ++; // try new row
        pcol = pc;
        //search col, give up if u reach last col
        while( (mat(prow,pcol)==0) & (pcol < mat.n)){
            pcol++;
        }
    }while((pcol == mat.n) && !(prow == mat.m-1));

    // if col is outside limit, whole sublock is 0
    if(pcol == mat.n){
        prow=mat.m;
    }
    return std::make_tuple(prow,pcol);
}


// Find the L EL U P factorization
template<typename F>
auto LEUP_fact(A<Dense<F>>& mat_arg){
    //create copy
    auto mat = mat_arg.copy();
    size_t m = mat.nrow();
    size_t n = mat.ncol();
    size_t p_col;
    size_t pr,pc;
    //return values
    L<Dense<F>> Lmat(m,m);
    EL<Dense<F>> ELmat(m,n);
    P<Dense<F>> Pmat(n,n);
    U<Dense<F>> Umat(n,n);

    std::vector<std::pair<size_t,size_t>> pivots;

    make_diag_ones(Lmat);
    make_diag_ones(Pmat);
    make_diag_ones(Umat);

    for(pr=0,pc=0;pc<n;pc++,pr++){
        //if pr++ exceeds range
        if(pr==m)
            break;
        // find pivot in lower right region of pr,pc
        std::tie(pr,p_col)=find_pivot(mat,pr,pc);
        // if full zero sublock then done
        if(pr==m)
            break;
        // apply and record permutation
        if(pc != p_col){
            mat.swap_cols(pc,p_col);
            Pmat.swap_rows(pc,p_col); // !inefficient!
        }

        //zero out column pc below pr - apply schur complement
        for(size_t i=pr+1;i<m;i++){
            if( !(mat(i,pc)==0) ){
                auto coef = mat(i,pc)/mat(pr,pc);
                mat.add_row_to(i,pr, -coef);
                //Lmat.add_col_to(pr,i, coef); //inefficient
                Lmat(i,pr) = coef; // record L
            }
        }
        //record pivot
        pivots.push_back(std::make_pair(pr,pc));
        ELmat(pr,pc)=1;

    }
    //extract U from upper echelon mat
    for(auto pair:pivots){
        std::tie(pr,pc) = pair;
        Umat.r(pc) = mat.r(pr);
    }

	mat.free();
    return std::make_tuple( Lmat, ELmat, Umat, Pmat );
}
