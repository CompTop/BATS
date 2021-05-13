#pragma once

#include <assert.h>
#include <unordered_map>
#include <vector>
#include <util/common.hpp>
#include <iostream>
#include <algorithm>
#include "matrix_interface.hpp"
#include <string>

//required for fill_rand<Rational<int>>()
#include <linalg/field.hpp>

/// @cond
// this is all obsolete and should not be in documentation

// storage type
template<typename F,typename Acc>
struct Dense{};

// memory access types
struct RowMaj{};
struct ColMaj{};
struct RevRowMaj{};
struct RevColMaj{};

// represents either a row view or column view
template<typename T>
struct VectorView{
    T* first;
    T* last;
	int stride;

    VectorView(T* start, T* end, size_t s) : first(start), last(end), stride(s) {};

	class Iterator {
	private:
		T* ptr;
		int stride;
	public:
		Iterator(T* p, const int s) noexcept : ptr(p), stride(s) { }

		// Prefix ++ overload
        Iterator& operator++() {
            ptr+=stride;
            return *this;
        }

		// Postfix ++ overload
		Iterator operator++(int) {
			Iterator iterator = *this;
			++*this;
			return iterator;
		}

		inline Iterator operator+(const size_t offset) const {
			return Iterator(ptr + stride*offset, stride);
		}

		inline bool operator!=(const Iterator& other) { return ptr != other.ptr;}
		inline size_t operator-(const Iterator& other) const {
			return (ptr - other.ptr)/stride;
		}

		inline T operator*() { return *ptr;}

	};

	Iterator begin() { return Iterator(first, stride); }
	Iterator end() { return Iterator(last, stride); }

    // self += a * x
    void axpy(const T a, const VectorView<T> x){
        T* ptr = first;
        T* xptr = x.first;
        while (ptr != last) {
            *ptr += a * (*xptr);
            ptr+=stride;
            xptr+=x.stride;
        }
        return;
    }

    VectorView& operator=(const VectorView<T> x){
        T* ptr = first;
        T* xptr = x.first;
        while (ptr != last) {
            *ptr = (*xptr);
            ptr+=stride;
            xptr+=x.stride;
        }
        return *this;
    }

    inline size_t size() const {
		return (last - first)/stride ;
	}

    inline T& operator[](size_t i) {return *(first + stride*i); }
	inline const T& operator[](size_t i) const {return *(first + stride*i); }
	inline T& operator()(size_t i) {return *(first + stride*i); }
	inline const T& operator()(size_t i) const {return *(first + stride*i); }

	T norm() const {
		// 2-norm of vector
		T n = T(0);
		T* ptr = first;
		while (ptr != last) {
			n += (*ptr)*(*ptr);
			ptr += stride;
		}
		return std::sqrt(n);
	}

	// basic arithmetic
	VectorView& operator/=(const T c) {
		T* ptr = first;
		while (ptr != last) {
			*ptr /= c;
			ptr += stride;
		}
		return *this;
	}
};

template<typename T>
inline T norm(const VectorView<T> &v) { return v.norm(); }

template<typename F, typename Acc>
struct MemAcc{
    size_t m,n;
    F* mat;

    MemAcc(size_t mm, size_t nn, F* mat) : m(mm), n(nn), mat(mat) {}
	// null constructor
	MemAcc(){ m=0;n=0;mat=NULL;}

	void init(size_t mm, size_t nn, F* mmat){
		m=mm;
		n=nn;
		mat=mmat;
	}

	// single element access
	inline F& operator()(const size_t k) { return mat[k]; }
	inline const F& operator()(const size_t k) const { return mat[k]; }

    inline F& operator()(const int i, const int j) {
        if constexpr(std::is_same<Acc,ColMaj>::value)
            return mat[j*m+i];
        else if constexpr(std::is_same<Acc,RowMaj>::value)
            return mat[i*n+j];
        else if constexpr(std::is_same<Acc,RevColMaj>::value)
            return mat[(n-j-1)*m+(m-i-1)];
        else if constexpr(std::is_same<Acc,RevRowMaj>::value)
            return mat[(m-i-1)*n+(n-j-1)];
    }

	inline const F& operator()(const int i, const int j) const {
		if constexpr(std::is_same<Acc,ColMaj>::value)
			return mat[j*m+i];
		else if constexpr(std::is_same<Acc,RowMaj>::value)
			return mat[i*n+j];
		else if constexpr(std::is_same<Acc,RevColMaj>::value)
			return mat[(n-j-1)*m+(m-i-1)];
		else if constexpr(std::is_same<Acc,RevRowMaj>::value)
			return mat[(m-i-1)*n+(n-j-1)];
	}

    // return a column view of column j
    inline VectorView<F> operator[](const size_t j) {
        if constexpr(std::is_same<Acc,ColMaj>::value)
            return VectorView<F>(mat + m*j, mat + m*(j+1),1);
        else if constexpr(std::is_same<Acc,RowMaj>::value)
            return VectorView<F>(mat + j, mat + m*n + j, n );
        else if constexpr(std::is_same<Acc,RevColMaj>::value)
            return VectorView<F>( mat + m*(n-j)-1,mat + m*(n-1-j)-1,-1);
        else if constexpr(std::is_same<Acc,RevRowMaj>::value)
            return VectorView<F>(mat + m*n -j-1, mat -1-j, -n );
    }

	// return a column view of column j
	inline const VectorView<F> operator[](const size_t j) const {
		if constexpr(std::is_same<Acc,ColMaj>::value)
			return VectorView<F>(mat + m*j, mat + m*(j+1),1);
		else if constexpr(std::is_same<Acc,RowMaj>::value)
			return VectorView<F>(mat + j, mat + m*n + j, n );
		else if constexpr(std::is_same<Acc,RevColMaj>::value)
			return VectorView<F>( mat + m*(n-j)-1,mat + m*(n-1-j)-1,-1);
		else if constexpr(std::is_same<Acc,RevRowMaj>::value)
			return VectorView<F>(mat + m*n -j-1, mat -1-j, -n );
	}

	// return a view of row i
    inline VectorView<F> r(const size_t i) {
        if constexpr(std::is_same<Acc,ColMaj>::value)
            return VectorView<F>(mat + i, mat + m*n + i, m );
        else if constexpr(std::is_same<Acc,RowMaj>::value)
            return VectorView<F>(mat + n*i, mat + n*(i+1),1);
        else if constexpr(std::is_same<Acc,RevColMaj>::value)
            return VectorView<F>(mat + m*n -i-1, mat -1-i, -m );
        else if constexpr(std::is_same<Acc,RevRowMaj>::value)
            return VectorView<F>( mat + n*(m-i)-1,mat + n*(m-1-i)-1,-1);
    }

	// return a view of row i
	inline const VectorView<F> r(const size_t i) const {
		if constexpr(std::is_same<Acc,ColMaj>::value)
			return VectorView<F>(mat + i, mat + m*n + i, m );
		else if constexpr(std::is_same<Acc,RowMaj>::value)
			return VectorView<F>(mat + n*i, mat + n*(i+1),1);
		else if constexpr(std::is_same<Acc,RevColMaj>::value)
			return VectorView<F>(mat + m*n -i-1, mat -1-i, -m );
		else if constexpr(std::is_same<Acc,RevRowMaj>::value)
			return VectorView<F>( mat + n*(m-i)-1,mat + n*(m-1-i)-1,-1);
	}

    void print(){
        for( size_t i=0; i<m; i++){
            for( size_t j=0; j<n; j++){
                std::cout<<((*this)(i,j))<<" ";
            }
            std::cout<<"\n";
        }
    }

	std::string str(){
		std::ostringstream oss;
		for( size_t i=0; i<m; i++){
            for( size_t j=0; j<n; j++){
                oss<<((*this)(i,j))<<" ";
            }
            oss<<"\n";
        }
	    return oss.str();
	}
};



// define the type changes for trp, Jconj versions
template<typename T>
struct Trp_T{};

template<typename T>
struct JConj_T{};

template<typename T>
struct TJConj_T{};

#define RULE(K,A,B) \
template<> \
struct K<A>{ \
    using type = B; \
};

//memory accessor rules

RULE( Trp_T, RowMaj, ColMaj )
RULE( Trp_T, ColMaj, RowMaj )
RULE( Trp_T, RevRowMaj, RevColMaj )
RULE( Trp_T, RevColMaj, RevRowMaj )

RULE( JConj_T, RowMaj, RevRowMaj )
RULE( JConj_T, ColMaj, RevColMaj )
RULE( JConj_T, RevRowMaj, RowMaj )
RULE( JConj_T, RevColMaj, ColMaj )

RULE( TJConj_T, RowMaj, RevColMaj )
RULE( TJConj_T, ColMaj, RevRowMaj )
RULE( TJConj_T, RevRowMaj, ColMaj )
RULE( TJConj_T, RevColMaj, RowMaj )

//matrix shape rules

#define SRULE(K,A,B) \
template<typename F, typename Acc> \
struct K<A<Dense<F,Acc>>>{ \
    using type = B<Dense<F,typename K<Acc>::type>>; \
};

SRULE( Trp_T, A, A )
SRULE( Trp_T, P, P )
SRULE( Trp_T, L, U )
SRULE( Trp_T, U, L )
SRULE( Trp_T, EL, EU )
SRULE( Trp_T, EU, EL )
SRULE( Trp_T, ELH, EUH )
SRULE( Trp_T, EUH, ELH )

SRULE( JConj_T, A, A )
SRULE( JConj_T, P, P )
SRULE( JConj_T, L, U )
SRULE( JConj_T, U, L )
SRULE( JConj_T, EL, EUH )
SRULE( JConj_T, EU, ELH )
SRULE( JConj_T, ELH, EU )
SRULE( JConj_T, EUH, EL )

SRULE( TJConj_T, A, A )
SRULE( TJConj_T, P, P )
SRULE( TJConj_T, L, L )
SRULE( TJConj_T, U, U )
SRULE( TJConj_T, EL, ELH )
SRULE( TJConj_T, EU, EUH )
SRULE( TJConj_T, ELH, EL )
SRULE( TJConj_T, EUH, EU )

SRULE( Trp_T, T, T )
SRULE( JConj_T, T, T )
SRULE( TJConj_T, T, T )

#define IMPLEMENT_TRP(M) \
auto Trp(){ \
		typename Trp_T<M>::type rmat(Base::n,Base::m,Base::mat); \
		return rmat; \
} \

#define IMPLEMENT_JCONJ(M) \
auto JConj(){ \
		typename JConj_T<M>::type rmat(Base::m,Base::n,Base::mat); \
		return rmat; \
} \

#define IMPLEMENT_TJCONJ(M) \
auto TJConj(){ \
		typename TJConj_T<M>::type rmat(Base::n,Base::m,Base::mat); \
		return rmat; \
} \


template<typename F,typename Acc>
struct A<Dense<F,Acc>>{
    using DI = Dense<F,Acc>;
	// to make things consistent with all the derived classes, Makes the macro work
	using Base = A<DI>;
	using FieldType = F;

    size_t m,n;
    F* mat = NULL;
	MemAcc<F,Acc> macc;

	// null constructor
	A<DI>(){ m=0;n=0;}

	//"copy" constructor
	A<DI>(const A<DI> &m2) {m = m2.m; n = m2.n; mat = m2.mat; macc=m2.macc; }

    A<DI>(size_t mm, size_t nn, F* mat) : m(mm), n(nn), mat(mat) {
		macc.init(mm,nn,mat);
	}

    A<DI>(size_t mm, size_t nn) : m(mm), n(nn) {
        mat = new F[m*n];
        // std::fill
        std::fill(mat, mat + m*n, F(0));
		macc.init(m,n,mat);
    }

    A<DI> copy() {
        F* newmat = new F[m*n];
        std::copy(mat, mat + m*n, newmat);

        return A<DI>(m, n, newmat);
    }

    inline size_t nrow() const { return m; }
    inline size_t ncol() const { return n; }
	inline size_t size() const { return m * n; }

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
        for( size_t i=0; i<n; i++){
            for( size_t j=0; j<m; j++){
                std::cout<<((*this)(j,i))<<",";
            }
            std::cout<<"\n";
        }
		std::cout<<"};";
    }

	std::string str(){
		std::ostringstream oss;
		for( size_t i=0; i<m; i++){
            for( size_t j=0; j<n; j++){
                oss<<((*this)(i,j))<<" ";
            }
            oss<<"\n";
        }
	    return oss.str();
	}

    // return a column view of column j
    inline VectorView<F> operator[](size_t j) {
        return macc[j];
    }
	inline const VectorView<F> operator[](size_t j) const {
		return macc[j];
	}

	// return a row view of row i
    inline VectorView<F> r(size_t i) {
        return macc.r(i);
    }
	inline const VectorView<F> r(size_t i) const {
        return macc.r(i);
    }

	inline F& operator()(const size_t k) { return macc(k); }
	inline const F& operator()(const size_t k) const { return macc(k); }
    inline F& operator()(int i, int j) { return macc(i,j); }
	inline const F& operator()(int i, int j) const { return macc(i,j); }

    void free(){
        delete [] mat;
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

	static A<DI> identity(size_t n) {
		A<DI> mat = A<DI>(n,n);
		for (size_t j = 0; j < n; j++) {
			mat(j,j) = F(1);
		}
		return mat;
	}


	IMPLEMENT_TRP(A<DI>)
	IMPLEMENT_JCONJ(A<DI>)
	IMPLEMENT_TJCONJ(A<DI>)

};



// Inherit all the constructors and implement mem acc transforms
#define INHERIT(T1,T2) \
template<typename F,typename Acc> \
struct T1<Dense<F,Acc>>:T2<Dense<F,Acc>>{ \
	using Base = T2<Dense<F,Acc>>; \
	using Self = T1<Dense<F,Acc>>; \
	using Base::Base ; \
	IMPLEMENT_TRP(Self)\
	IMPLEMENT_JCONJ(Self)\
	IMPLEMENT_TJCONJ(Self)\
}; \

INHERIT(T,A)
INHERIT(L,T)
INHERIT(U,T)
INHERIT(E,A)
INHERIT(EL,E)
INHERIT(EU,E)
INHERIT(ELH,E)
INHERIT(EUH,E)
INHERIT(P,E)


// naive matmul
template< typename F,typename Acc1, typename Acc2>
auto matmul(A<Dense<F,Acc1>> m1, A<Dense<F,Acc2>> m2){
    assert( m1.n==m2.m );
    A<Dense<F,Acc1>> prod(m1.m,m2.n);
    for(size_t i=0;i<m1.m;i++)
        for(size_t j=0;j<m2.n;j++)
            for(size_t k=0;k<m1.n;k++){
                prod(i,j) += m1(i,k)*m2(k,j);
            }
    return prod;
}

//inplace
template< typename F,typename Acc1>
void matmul(A<Dense<F,Acc1>> m1, VectorView<F> y){
    F temp[m1.m];
    A<Dense<F,Acc1>> prod(m1.m,1,temp);
    for(size_t i=0;i<m1.m;i++)
            for(size_t k=0;k<m1.n;k++){
                prod(i,0) += m1(i,k)*y[k];
            }
    y=prod[0];
}

template< typename F,typename Acc1>
void matmul(VectorView<F> y, A<Dense<F,Acc1>> m1 ){
    F temp[m1.n];
    A<Dense<F,Acc1>> prod(1,m1.n,temp);
    for(size_t j=0;j<m1.n;j++)
            for(size_t k=0;k<m1.n;k++){
                prod(0,j) += y[k]*m1(k,j);
            }
    y=prod.r(0);
}

template <typename F,typename Acc1,typename Acc2>
A<Dense<F,Acc2>> apply_matmul_on_left(A<Dense<F,Acc1>> Amat1, A<Dense<F,Acc2>> Amat2) {
    size_t n = Amat2.ncol();
    for (size_t j = 0; j < n; j++) {
        matmul(Amat1, Amat2[j]);
    }
    return Amat2;
}

template <typename F,typename Acc1,typename Acc2>
A<Dense<F,Acc2>> apply_matmul_on_right(A<Dense<F,Acc1>> Amat1, A<Dense<F,Acc2>> Amat2) {
    size_t m = Amat1.nrow();
    for (size_t i = 0; i < m; i++) {
        matmul(Amat1.r(i), Amat2);
    }
    return Amat2;
}


/*
Lret * ELmat = ELmat *Lmat
*/
template <typename F,typename Acc1,typename Acc2>
L<Dense<F,Acc2>> commute(EL<Dense<F,Acc1>> ELmat, L<Dense<F,Acc2>> Lmat) {
    size_t m = ELmat.nrow();
    size_t n = ELmat.ncol();
    assert(n == Lmat.ncol());
    assert(Lmat.ncol() == Lmat.nrow());
    L<Dense<F,Acc2>> Lret = L<Dense<F,Acc2>>(m, m);

    // step 1 is to make an index map from ELmat
    std::vector<size_t> idx_map(n);
    size_t i = 0; // pointer to current row
    size_t j = 0; // pointer to column
    // loop over columns of ELmat
    for (j = 0; j < n && i < m; j++) {
        // increment i until we find a non-zero, or run out of column
        while (ELmat(i,j) == F(0) && i < m) {
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


template <typename F,typename Acc1,typename Acc2>
auto commute(L<Dense<F,Acc2>> Lmat, ELH<Dense<F,Acc1>> ELHmat ) {
    auto ELm = ELHmat.TJConj();
    auto Lm = Lmat.TJConj();
    auto Lp = commute(ELm,Lm);
    return Lp.TJConj();
}

template <typename F,typename Acc1,typename Acc2>
auto commute(EUH<Dense<F,Acc2>> EUHmat, U<Dense<F,Acc1>> Umat ) {
    auto ELm = EUHmat.JConj();
    auto Lm = Umat.JConj();
    auto Lp = commute(ELm,Lm);
    return Lp.JConj();
}

template <typename F,typename Acc1,typename Acc2>
auto commute(U<Dense<F,Acc1>> Umat, EU<Dense<F,Acc2>> EUmat ) {
    auto ELm = EUmat.Trp();
    auto Lm = Umat.Trp();
    auto Lp = commute(ELm,Lm);
    return Lp.Trp();
}

// solve Lmat x = y in-place
template <typename F,typename Acc>
void l_solve(L<Dense<F,Acc>> Lmat, VectorView<F> y) {
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
template <typename F,typename Acc1,typename Acc2>
A<Dense<F,Acc2>> apply_inverse_on_left(L<Dense<F,Acc1>> Lmat, A<Dense<F,Acc2>> Amat) {

    //size_t m = Amat.nrow();
    size_t n = Amat.ncol();

    //A<Dense<F,Acc2>> Aret = Amat.copy();

    // apply lower triangular solve to each column
    for (size_t j = 0; j < n; j++) {
        l_solve(Lmat, Amat[j]);
    }

    return Amat;

}


// implementations of naive matrix ops

template<typename F,typename Acc>
void A<Dense<F,Acc>>::swap_cols(size_t i, size_t j){
    for(size_t k=0;k<m;k++){
        std::swap((*this)(k,i),(*this)(k,j));
    }
}

template<typename F,typename Acc>
void A<Dense<F,Acc>>::swap_rows(size_t i, size_t j){
    for(size_t k=0;k<n;k++){
        std::swap((*this)(i,k),(*this)(j,k));
    }
}

template<typename F,typename Acc>
void A<Dense<F,Acc>>::add_col_to(size_t i, size_t j, F a){
    for(size_t k=0;k<m;k++){
       (*this)(k,i) = (*this)(k,i) + a*(*this)(k,j);
    }
}

template<typename F,typename Acc>
void A<Dense<F,Acc>>::add_row_to(size_t i, size_t j, F a){
    for(size_t k=0;k<n;k++){
       (*this)(i,k) = (*this)(i,k) + a*(*this)(j,k);
    }
}


// mat creation utils

// template<typename M>
// void make_diag_ones(M mat){
//     for(size_t i=0;i<mat.m;i++)
//         mat(i,i)=1;
// }

template<typename F,typename M>
void make_diag_ones(M mat){
    for(size_t i=0;i<mat.m;i++)
        mat(i,i)=F(1);
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

template<typename F,typename M>
void fill_rand(M mat, int vmax=10){
	// fill put random numbers in matrix
    for(size_t i=0;i<mat.m;i++)
        for(size_t j=0;j<mat.n;j++)
			mat(i,j)=F(rand()%vmax);
}



// Find pivot in lower right block of pr,pc
template<typename M>
auto find_pivot(M mat,size_t pr, size_t pc){
    size_t prow = pr-1; // start 1 less as its immediately incremented
    size_t pcol = pc;
	using F = typename M::FieldType;

    // outer loop searches rows
    do{
        prow ++; // try new row
        pcol = pc;
        //search col, give up if u reach last col
        while( (mat(prow,pcol)==F(0)) & (pcol < mat.n)){
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
template<typename F,typename Acc>
auto LEUP_fact(A<Dense<F,Acc>>& mat_arg){
    //create copy
    auto mat = mat_arg.copy();
    size_t m = mat.nrow();
    size_t n = mat.ncol();
    size_t p_col;
    size_t pr,pc;
    //return values
    L<Dense<F,Acc>> Lmat(m,m);
    EL<Dense<F,Acc>> ELmat(m,n);
    P<Dense<F,Acc>> Pmat(n,n);
    U<Dense<F,Acc>> Umat(n,n);

    std::vector<std::pair<size_t,size_t>> pivots;

    make_diag_ones<F>(Lmat);
    make_diag_ones<F>(Pmat);
    make_diag_ones<F>(Umat);

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
            if( !(mat(i,pc)==F(0)) ){
                auto coef = mat(i,pc)/mat(pr,pc);
                mat.add_row_to(i,pr, -coef);
                //Lmat.add_col_to(pr,i, coef); //inefficient
                Lmat(i,pr) = coef; // record L
            }
        }
        //record pivot
        pivots.push_back(std::make_pair(pr,pc));
        ELmat(pr,pc)=F(1);

    }
    //extract U from upper echelon mat
    for(auto pair:pivots){
        std::tie(pr,pc) = pair;
        Umat.r(pc) = mat.r(pr);
    }

	mat.free();
    return std::make_tuple( Lmat, ELmat, Umat, Pmat );
}


template<typename F,typename Acc>
auto PLEU_fact(A<Dense<F,Acc>>& mat_arg){
    auto  mat_trp = mat_arg.Trp();
    auto [Lmat, ELmat, Umat, Pmat] = LEUP_fact(mat_trp);
    auto Pret = Pmat.Trp();
    auto Lret = Umat.Trp();
    auto EUret = ELmat.Trp();
    auto Uret = Lmat.Trp();
    return std::make_tuple( Pret,Lret ,EUret,Uret );
}

template<typename F,typename Acc>
auto UELP_fact(A<Dense<F,Acc>>& mat_arg){
    auto  mat_jconj = mat_arg.JConj();
    auto [Lmat, ELmat, Umat, Pmat] = LEUP_fact(mat_jconj);
    auto Uret = Lmat.JConj();
    auto EUHret = ELmat.JConj();
    auto Lret = Umat.JConj();
    auto Pret = Pmat.JConj();
    return std::make_tuple( Uret,EUHret ,Lret,Pret );
}

template<typename F,typename Acc>
auto PUEL_fact(A<Dense<F,Acc>>& mat_arg){
    auto  mat_tjconj = mat_arg.TJConj();
    auto [Lmat, ELmat, Umat, Pmat] = LEUP_fact(mat_tjconj);
    auto Pret = Pmat.TJConj();
    auto Uret = Umat.TJConj();
    auto ELHret = ELmat.TJConj();
    auto Lret = Lmat.TJConj();
    return std::make_tuple( Pret,Uret ,ELHret,Lret );
}


template <typename F,typename Acc1,typename Acc2>
auto apply_inverse_on_left(U<Dense<F,Acc1>> Umat, A<Dense<F,Acc2>> Amat) {
    auto Amat_jconj = Amat.JConj();
    auto Lm = Umat.JConj();
    auto ret = apply_inverse_on_left(Lm,Amat_jconj);
    auto ret2 = ret.JConj();
    return ret2;
}

template <typename F,typename Acc1,typename Acc2>
auto apply_inverse_on_right(A<Dense<F,Acc2>> Amat,U<Dense<F,Acc1>> Umat) {
    auto Amat_trp = Amat.Trp();
    auto Lm = Umat.Trp();
    auto ret = apply_inverse_on_left(Lm,Amat_trp);
    auto ret2 = ret.Trp();
    return ret2;
}

template <typename F,typename Acc1,typename Acc2>
auto apply_inverse_on_right(A<Dense<F,Acc2>> Amat,L<Dense<F,Acc1>> Lmat) {
    auto Amat_tjconj = Amat.TJConj();
    auto Lm = Lmat.TJConj();
    auto ret = apply_inverse_on_left(Lm,Amat_tjconj);
    auto ret2 = ret.TJConj();
    return ret2;
}

/// @endcond
