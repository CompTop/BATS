#include <bats.hpp>
#include <iostream>

#define FT ModP<int, 3>

int main() {

	bats::future::Matrix A(5,5, FT(1));
	A(1,0) = 0;
	A(2,1) = 2;
	A(3,2) = 2;
	A.print();

	// auto H1 = HessenbergTransform(A);
	// H1.print();
	// H1.prod().print();
	//
	// std::cout << "hessenberg_to_companion" << std::endl;
	//
	// hessenberg_to_companion(H1);
	// H1.print();
	// H1.prod().print();

	auto F = bats::future::LU(A);
	F.print();
	std::cout << "product:" << std::endl;
	auto P = F.prod();
	P.print();

	size_t test[2] = {1,2};
	std::cout << test[0] << test[1] << std::endl;

	P = bats::future::Matrix<FT>(P, bats::future::ColumnMajor());
	P.append_column();
	P.print();

	auto r = bats::future::range(0, 4);
	std::vector<size_t> c = {0,2,3} ;
	// auto view = MatrixView(P, r, c);
	auto view = P.view(r,c);
	// auto view = MatrixView(P, range(0, 4), range(0, 3));

	view.print();

	auto S = bats::future::Span(2, FT());
	S.Pt.print();
	S.L.print();

	std::vector<FT> v = {1,1};
	std::cout << S.add(v) << std::endl;

	S.Pt.print();
	S.L.print();

	std::cout << S.add(v) << std::endl;

	S.Pt.print();
	S.L.print();

	v = {1,0};

	std::cout << S.add(v) << std::endl;

	S.Pt.print();
	S.L.print();
	// std::cout << S.L.column(0)[1] << std::endl;

	std::cout << S.contains(v) << std::endl;

	v = {0,1};

	std::cout << S.contains(v) << std::endl;

	return 0;
}
