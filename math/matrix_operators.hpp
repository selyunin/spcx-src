/*
 * matrix_operators.hpp
 *
 *  Created on: Oct 27, 2010
 *      Author: frehse
 */

#ifndef MATRIX_OPERATORS_HPP_
#define MATRIX_OPERATORS_HPP_

#include "matrix_operators.h"

#include "math/matrix.h"
#include "math/matrix_utility.h"
#include "math/vector.h"
#include "math/unique_vector_to_value_store.h"
#include "utility/logger_stopwatch.h"
#include "extern/ublasJama-1.0.2.2/EigenvalueDecomposition.hpp"
#include "math/numeric/container_comp.h" // for lex_comp_less

#include "math/numeric/invalid_number_exception.h"

namespace math {

template<typename scalar_type> bool operator==(const matrix<scalar_type>& m1,
		const matrix<scalar_type>& m2) {
	if (m1.size1() != m2.size1() || m1.size2() != m2.size2()) {
		return false;
	}
	for (unsigned int i = 0; i < m1.size1(); i++) {
		for (unsigned int j = 0; j < m1.size2(); j++) {
			if (m1(i, j) != m2(i, j)) {
				return false;
			}
		}
	}
	return true;
}

template<typename scalar_type> bool operator!=(const matrix<scalar_type>& m1,
		const matrix<scalar_type>& m2) {
	return !(m1 == m2);
}

template<typename scalar_type> matrix<scalar_type> operator*(
		const scalar_type& s, const matrix<scalar_type>& M) {
	matrix<scalar_type> res = M;
	res.multiply_scalar(s);
	return res;
}

template<typename scalar_type> matrix<scalar_type> operator*(
		const matrix<scalar_type>& M, const scalar_type& s) {
	matrix<scalar_type> res = M;
	res.multiply_scalar(s);
	return res;
}

template<typename scalar_type> vector<scalar_type> operator*(const matrix<
		scalar_type>& M, const vector<scalar_type>& v) {
	assert(M.size2()==v.size());
	return M.multiply_vector(v);
}

template<typename scalar_type> vector<scalar_type> operator*(const vector<
		scalar_type>& v, const matrix<scalar_type>& M) {
	assert(M.size1()==v.size());
	return vector<scalar_type> (prod(v.get_vector_impl(), M.get_matrix_impl()));
}

template<typename scalar_type> matrix<scalar_type> operator*(const matrix<
		scalar_type>& M1, const matrix<scalar_type>& M2) {
	assert(M1.size2()==M2.size1());
	return matrix<scalar_type> (
			prod(M1.get_matrix_impl(), M2.get_matrix_impl()));
}
template<typename scalar_type> matrix<scalar_type> operator+(const matrix<
		scalar_type>& M1, const matrix<scalar_type>& M2 ) {
	assert(M1.size1() == M2.size1() && M1.size2() == M2.size2());
	return matrix<scalar_type>(M1.get_matrix_impl() + M2.get_matrix_impl());
}
template<typename scalar_type> matrix<scalar_type> operator-(const matrix<
		scalar_type>& M1, const matrix<scalar_type>& M2 ) {
	assert(M1.size1() == M2.size1() && M1.size2() == M2.size2());
	return matrix<scalar_type>(M1.get_matrix_impl() - M2.get_matrix_impl());
}

template<typename scalar_type> matrix<scalar_type> operator-(
		const matrix<scalar_type>& M) {
	return matrix<scalar_type>(-M.get_matrix_impl());
}

template<typename scalar_type> scalar_type matrix_determinant(
		const matrix<scalar_type>& A) {
	return boost::numeric::ublas::matrix_determinant(A.get_matrix_impl());
}

template<typename scalar_type> bool check_consistency(const matrix<scalar_type>& A) {
	for (typename matrix<scalar_type>::size_type i=0;i<A.size1();++i) {
		for (typename matrix<scalar_type>::size_type j=0;j<A.size2();++j) {
			if (!is_finite(A(i,j))) {
				return false;
			}
		}
	}
	return true;
}

template<typename scalar_type> matrix<scalar_type> matrix_exponential(
		const matrix<scalar_type>& A) {
	LOGGERSWOC(DEBUG4,"matrix_exponential","Computing matrix exponential");

	// Cache matrices
	//typedef typename matrix<scalar_type>::matrix_impl_type::array_type cache_type;
	typedef float cache_scalar_type;
	static unique_vector_to_value_store<cache_scalar_type,matrix,matrix<scalar_type> > exp_store;
	matrix<cache_scalar_type> A_small = A.template convert_to<cache_scalar_type>();
	typename unique_vector_to_value_store<cache_scalar_type,matrix,matrix<scalar_type> >::const_iterator it = exp_store.find(A_small);

	if (it == exp_store.end()) {
		typedef boost::numeric::ublas::compressed_matrix<scalar_type> compressed_type;

//		compressed_type A_compr=get_compressed(A);
//				matrix<scalar_type> res(
//						boost::numeric::ublas::expm_pad(A_compr));
		matrix<scalar_type> res;
		{
			LOGGERSWOC(DEBUG6,"matrix_exponential","Calculating matrix exponential not in cache");
			res = matrix<scalar_type>(boost::numeric::ublas::expm_pad(A.get_matrix_impl()));
		}
		if (res.size1()!=A.size1() || res.size2()!=A.size2()) {
			throw std::runtime_error("matrix exponential gives matrix of wrong size ("+to_string(res.size1())+"x"+to_string(res.size2())+") instead of "+to_string(A.size1()));
		}
		if (!check_consistency(res))  {
			throw invalid_number_exception("matrix exponential gives inf or nan");
		}
		exp_store.insert_missing(A_small, res);
		return res;
	} else {
		const matrix<scalar_type>& res = it->second;
//		numeric::lex_comp_less<cache_scalar_type, matrix> comparator;
//		std::cout << comparator(A_small,it->first) << comparator(it->first,A_small) << std::endl;
		if (res.size1()!=A.size1() || res.size2()!=A.size2()) {
			throw std::runtime_error("matrix exponential gives matrix of wrong size ("+to_string(res.size1())+"x"+to_string(res.size2())+") instead of "+to_string(A.size1()));
		}
		return res;
	}
}

template<typename scalar_type> boost::numeric::ublas::compressed_matrix<scalar_type> get_compressed(
		const matrix<scalar_type>& A) {
	size_t count = 0;
	for (size_t i = 0; i < A.size1(); ++i) {
		for (size_t j = 0; j < A.size2(); ++j) {
			if (!numeric::is_MEQ(A(i,j),scalar_type(0))) {
				++count;
			}
		}
	}
	boost::numeric::ublas::compressed_matrix<scalar_type> M(A.size1(),A.size2(),count);
	for (size_t i = 0; i < A.size1(); ++i) {
		for (size_t j = 0; j < A.size2(); ++j) {
			if (!numeric::is_MEQ(A(i,j),scalar_type(0))) {
				M.push_back(i,j,A(i,j));
			}
		}
	}
	return M;
}

template<> inline matrix<Rational> matrix_exponential(const matrix<Rational>& A) {
	matrix<double> m = A.convert_to<double> ();
	m = matrix_exponential(m);
	return m.convert_to<Rational> ();
}

template<typename scalar_type> matrix<scalar_type> absolute(const matrix<scalar_type>& A) {
	matrix<scalar_type> a(A.size1(),A.size2());
	for(unsigned int i=0;i<a.size1();i++)
		for(unsigned int j=0;j<a.size2();j++)
			a(i,j) = abs(A(i,j));
	return a;
}

template<typename scalar_type> void get_special_matrices(const matrix<scalar_type>& A, const scalar_type& delta,
		matrix<scalar_type>& special_matrix1,
		matrix<scalar_type>& special_matrix2,
		matrix<scalar_type>& special_matrix3){
	if(A.size1() != A.size2())
		throw std::runtime_error("math::get_special_matrices: Dont know what to do with non-square matrices\n");
	namespace b = boost::numeric::ublas ;
	std::size_t d = A.size1();
	math::matrix<scalar_type> M(d*3,d*3);

	//M.submatrix_assign(absolute(A),0,d,0,d);
	M.submatrix_assign(A,0,d,0,d);
	M.submatrix_assign(matrix<scalar_type>(b::identity_matrix<scalar_type>(d)),0,d,d,2*d);
	M.submatrix_assign(matrix<scalar_type>(b::identity_matrix<scalar_type>(d)),d,2*d,2*d,3*d);

	M.multiply_scalar(delta);
	M = matrix_exponential(M);

	special_matrix1 = M.project_submatrix(0,d,0,d);
	special_matrix2 = M.project_submatrix(0,d,d,2*d);
	special_matrix3 = M.project_submatrix(0,d,2*d,3*d);
}

template<typename scalar_type> matrix<scalar_type> augment_cols(const matrix<scalar_type>& A, const matrix<scalar_type>& B){
	if(A.size1() != B.size1()){
		throw std::runtime_error("Matrix: Augment: Matrices not compatible for column augmentation");
	}
	unsigned int m,n;
	m = A.size1();
	n = A.size2() + B.size2();
	matrix<scalar_type> aug_mat(m,n);
	for(unsigned int i=0;i<m;i++){
		for(unsigned int j=0;j<A.size2();j++){
			aug_mat(i,j) = A(i,j);
		}
	}
	for(unsigned int i=0;i<m;i++){
		for(unsigned int j=A.size2(),k=0;j<n;j++,k++){
			aug_mat(i,j) = B(i,k);
		}
	}

	return aug_mat;

}

template<typename scalar_type> matrix<scalar_type> augment_rows(const matrix<scalar_type>& A, const matrix<scalar_type>& B){
	if(A.size2() != B.size2()){
		throw std::runtime_error("Matrix: Augment_rows: Matrices not compatible for row augmentation");
	}
	unsigned int m,n;
	m = A.size1() + B.size2();
	n = A.size2();
	matrix<scalar_type> aug_mat(m,n);
	for(unsigned int i=0;i<A.size1();i++){
		for(unsigned int j=0;j<n;j++){
			aug_mat(i,j) = A(i,j);
		}
	}
	for(unsigned int i=B.size1(),k=0;i<m;i++,k++){
		for(unsigned int j=0;j<n;j++){
			aug_mat(i,j) = B(k,j);
		}
	}

	return aug_mat;

}
template<typename scalar_type> void compute_eigenvalues(const matrix<scalar_type>&A,
													    vector<scalar_type>& D_real,
														vector<scalar_type>& D_imag,
														matrix<scalar_type>& V){
	using namespace boost::numeric::ublas;

	if (A.size1()==0 || A.size2()==0) {
		D_real = vector<scalar_type>();
		D_imag = vector<scalar_type>();
		V = matrix<scalar_type>();
		return;
	}

	try {
		EigenvalueDecomposition ed(A.get_matrix_impl());

		D_real = vector<scalar_type>(ed.getRealEigenvalues());
		D_imag = vector<scalar_type>(ed.getImagEigenvalues());
		V = matrix<scalar_type>(ed.getV());
	} catch (std::runtime_error& e) {
		std::stringstream ss;
		ss << "Could not compute Eigenvalues of matrix: " << A;
		throw math_exception(ss.str(), e);
	}

}
template<> void compute_eigenvalues(const matrix<Rational>&A,
													    vector<Rational>& D_real,
														vector<Rational>& D_imag,
														matrix<Rational>& V){
	using namespace boost::numeric::ublas;
	matrix<double> A_d = A.convert_to<double>();
	EigenvalueDecomposition ed(A_d.get_matrix_impl());

	vector<double> d_real_d(ed.getRealEigenvalues());
	vector<double> d_imag_d(ed.getImagEigenvalues());
	matrix<double> v_d(ed.getV());

	D_real = d_real_d.convert_to<Rational>();
	D_imag = d_imag_d.convert_to<Rational>();
	V = v_d.convert_to<Rational>();
}

template<typename scalar_type> void minimize_norm(const matrix<scalar_type>&A,
		matrix<scalar_type>&M, matrix<scalar_type>&M_inv,
		matrix<scalar_type>& A_min,
		const scalar_type eps){

	assert(A.size1() == A.size2());
	assert(M_inv.size1() == M_inv.size2());
	assert(M.size1() == M.size2());

	assert(M_inv.size1() == A.size1());
	assert(M.size1() == A.size1());

	unsigned int n = A.size1();
	scalar_type r_sum = scalar_type(0), c_sum = scalar_type(0);
	math::vector<scalar_type> alpha_vec(n);

	scalar_type old_norm,new_norm,alpha;
	matrix<scalar_type> A_t(A);

	do{
		old_norm = A_t.infinity_norm();
//		std::cout << "old norm:" << old_norm << std::endl;

		for(unsigned int i=0;i<n;i++){
			r_sum = c_sum = scalar_type(0);

			for(unsigned int j=0;j<n;j++){
				r_sum = r_sum + std::abs(A_t(i,j));
				c_sum = c_sum + std::abs(A_t(j,i));

				assert(r_sum != scalar_type(0) && c_sum != scalar_type(0));

			}

			alpha = std::sqrt(c_sum/r_sum);
			alpha_vec[i] = alpha;
			A_t.multiply_row_scalar(alpha_vec[i],i);
			A_t.multiply_col_scalar(scalar_type(1)/alpha_vec[i],i);
			A_t(i,i) = A(i,i); // since the diagonal entries dont change,
							   // we put the original values to avoid floating point errors.

//			std::cout << "alpha:" << i << "=" <<alpha << std::endl;

		}
		new_norm = A_t.infinity_norm();
//		std::cout << "new norm:" << new_norm << std::endl;

	}while(old_norm - new_norm > eps);

	A_min = A_t;
	for(unsigned int i=0;i<n;i++){
		M_inv(i,i) = alpha_vec[i];
		M(i,i) = scalar_type(1)/alpha_vec[i];
	}
}
template<> inline void minimize_norm(const matrix<Rational>&A,
		matrix<Rational>&M, matrix<Rational>&M_inv,
		matrix<Rational>&A_min,
		const Rational eps){
	assert(A.size1() == A.size2());
	unsigned int n = A.size1();
	matrix<double> A_d = A.convert_to<double>();
	matrix<double> M_d(n,n);
	matrix<double> M_inv_d(n,n);
	matrix<double> A_min_d(n,n);
	minimize_norm<double>(A_d,M_d,M_inv_d,A_min_d,eps.get_double());
	M = M_d.convert_to<Rational>();
	M_inv = M_inv_d.convert_to<Rational>();
	A_min = A_min_d.convert_to<Rational>();
}

template<typename scalar_type> void reduced_row_echelon_form(matrix<scalar_type>& M) {
	using namespace numeric;
	size_t lead = 0;
    size_t rowCount = M.size1();
    size_t columnCount = M.size2();
    vector<scalar_type> row_vec(columnCount);
    for (size_t r = 0; r < rowCount && lead < columnCount; ++r) {
    	// Search for nonzero element on column lead
        size_t i = r;
        while (is_MEQ(M(i, lead),scalar_type(0))) {
            i = i + 1;
            if (i == rowCount) {
            	// whole column is zero, move lead to the right
                i = r;
                lead = lead + 1;
                if (lead == columnCount) {
                    return;
                }
            }
        }
        //Swap rows i and r
        M.row_swap(i,r);
        if (!is_MEQ(M(r, lead),scalar_type(0))) {
        	// divide row r by M[r, lead]
        	M.multiply_row_scalar(scalar_type(1/M(r, lead)),r);
        }
        for (i=0; i < rowCount; ++i) {
            if (i != r) {
                //Subtract M[i, lead] multiplied by row r from row i
            	row_vec = M.vector_from_row(i);
            	row_vec -= M(i,lead)*M.vector_from_row(r);
            	M.assign_row(i,row_vec);
            }
        }
        lead = lead + 1;
    }

}

}

#endif /* MATRIX_OPERATORS_HPP_ */
