#ifndef SPARSE_LU_H_
#define SPARSE_LU_H_

#include "sparse_matrix.h"
#include "symbolic_sparse.h"

template<typename T>
class SparseLU {

public:
    inline void SolverLU(const SparseMatrix<T>& spmat, const SymbolicSpMat& sym_spmat, std::vector<T>& vec_b) {
        /*LeftLookingLU(spmat, sym_spmat);
        SparseLTriSolver(sym_spmat, vec_b);
        SparseUTriSolver(sym_spmat, vec_b);*/

        LeftLookingLU(spmat.num_cols_, spmat.values_, spmat.col_compressed_slice_, spmat.row_indices_, sym_spmat.l_col_compressed_slice_, sym_spmat.l_row_indice_, sym_spmat.u_col_compressed_slice_, sym_spmat.u_row_indice_);
        SparseLTriSolver(vec_b.size(), sym_spmat.l_col_compressed_slice_, sym_spmat.l_row_indice_, vec_b);
        SparseUTriSolver(vec_b.size(), sym_spmat.u_col_compressed_slice_, sym_spmat.u_row_indice_, vec_b);
    }
    inline std::vector<T> SolverLU(const SparseMatrix<T>& spmat, const SymbolicSpMat& sym_spmat, const std::vector<T>& vec_b) {
        std::vector<T> x(vec_b);
        LeftLookingLU(spmat, sym_spmat);
        SparseLTriSolver(sym_spmat, x);
        SparseUTriSolver(sym_spmat, x);

        return x;
    }

    inline std::vector<T> SolverLU(const SparseMatrix<std::complex<double> >& spmat, const SymbolicSpMat& sym_spmat, std::vector<double>& vec_b) {
        std::vector<T> x(vec_b.size());
        for (int i = 0; i < vec_b.size(); ++i) {
            x[i] = vec_b[i];
        }
        LeftLookingLU(spmat, sym_spmat);
        SparseLTriSolver(sym_spmat, x);
        SparseUTriSolver(sym_spmat, x);

        return x;
    }

	inline std::vector<T> SolverLU(SparseMatrix<std::complex<double> >&& spmat, const SymbolicSpMat& sym_spmat, std::vector<double>& vec_b) {
		std::vector<T> x(vec_b.size());
		for (int i = 0; i < vec_b.size(); ++i) {
			x[i] = vec_b[i];
		}
		LeftLookingLU(std::forward<SparseMatrix<std::complex<double> > >(spmat), sym_spmat);
		SparseLTriSolver(sym_spmat, x);
		SparseUTriSolver(sym_spmat, x);

		return x;
	}

private:
    std::vector<T> l_value_vec_;
    std::vector<T> u_value_vec_;

    inline void LeftLookingLU(const SparseMatrix<T>& spmat, const SymbolicSpMat& sym_spmat) {
        const int mat_size = spmat.num_cols_;
        const std::vector<T>& value_vec = spmat.values_;
        const std::vector<int>& col_compressed_slice = spmat.col_compressed_slice_;
		const std::vector<int>& row_indices = spmat.row_indices_;

        const std::vector<int>& l_col_compressed_slice = sym_spmat.l_col_compressed_slice_;
        const std::vector<int>& l_row_indices = sym_spmat.l_row_indice_;
        const std::vector<int>& u_col_compressed_slice = sym_spmat.u_col_compressed_slice_;
        const std::vector<int>& u_row_indices = sym_spmat.u_row_indice_;


        l_value_vec_.resize(l_col_compressed_slice[mat_size]);
        u_value_vec_.resize(u_col_compressed_slice[mat_size]);

		std::vector<T> col_vec(mat_size);
        for (int i = 0; i < mat_size; ++i) {	
			/*col_vec.resize(mat_size);
			for (int j = col_compressed_slice[i]; j < col_compressed_slice[i + 1]; ++j) {
				col_vec[row_indices[j]] = value_vec[j];
			}*/
			for (int k = u_col_compressed_slice[i]; k < u_col_compressed_slice[i + 1] - 1; ++k) {
				col_vec[u_row_indices[k]] = 0;
			}
			for (int k = l_col_compressed_slice[i] + 1; k < l_col_compressed_slice[i + 1]; ++k) {
				col_vec[l_row_indices[k]] = 0;
			}
			for (int j = col_compressed_slice[i]; j < col_compressed_slice[i + 1]; ++j) {
				col_vec[row_indices[j]] = value_vec[j];
			}

            for (int k = u_col_compressed_slice[i]; k < u_col_compressed_slice[i + 1] - 1; ++k) {
                int col_index = u_row_indices[k];
                for (int m = l_col_compressed_slice[col_index] + 1; m < l_col_compressed_slice[col_index + 1]; ++m) {
                    int row_index = l_row_indices[m];
                    col_vec[row_index] -= l_value_vec_[m] * col_vec[col_index];
                }
            }

            for (int k = u_col_compressed_slice[i]; k < u_col_compressed_slice[i + 1]; ++k) {
                u_value_vec_[k] = col_vec[u_row_indices[k]];
            }
            
            for (int k = l_col_compressed_slice[i]; k < l_col_compressed_slice[i + 1]; ++k) {
                l_value_vec_[k] = col_vec[l_row_indices[k]] / col_vec[i];
            } 
        }
    }

    inline void SparseLTriSolver(const SymbolicSpMat& sym_spmat, std::vector<T>& vec_b) {
        const std::vector<int>& l_col_compressed_slice = sym_spmat.l_col_compressed_slice_;
        const std::vector<int>& l_row_indices = sym_spmat.l_row_indice_;

        int mat_size = l_col_compressed_slice.size() - 1;

        for (int i = 0; i < mat_size; ++i) {

            for (int j = l_col_compressed_slice[i] + 1; j < l_col_compressed_slice[i + 1]; ++j) {
                vec_b[l_row_indices[j]] -= l_value_vec_[j] * vec_b[i];
            }
        }
    }
    
    inline void SparseUTriSolver(const SymbolicSpMat& sym_spmat, std::vector<T>& vec_b) {
        const std::vector<int>& u_col_compressed_slice = sym_spmat.u_col_compressed_slice_;
        const std::vector<int>& u_row_indices = sym_spmat.u_row_indice_;

        int mat_size = u_col_compressed_slice.size() - 1;

        for (int i = mat_size - 1; i >= 0; --i) {
            vec_b[i] /= u_value_vec_[u_col_compressed_slice[i + 1] - 1];

            for (int j = u_col_compressed_slice[i]; j < u_col_compressed_slice[i + 1] - 1; ++j) {
                vec_b[u_row_indices[j]] -= u_value_vec_[j] * vec_b[i];
            }
        }
    }

    inline void LeftLookingLU(const int mat_size, const std::vector<T>& value_vec, const std::vector<int>& col_compressed_slice, const std::vector<int>& row_indices, const std::vector<int>& l_col_compressed_slice, const std::vector<int>& l_row_indices,
        const std::vector<int>& u_col_compressed_slice, const std::vector<int>& u_row_indices) {
        /*const int mat_size = spmat.num_cols_;
        const std::vector<T>& value_vec = spmat.values_;
        const std::vector<int>& col_compressed_slice = spmat.col_compressed_slice_;
        const std::vector<int>& row_indices = spmat.row_indices_;

        const std::vector<int>& l_col_compressed_slice = sym_spmat.l_col_compressed_slice_;
        const std::vector<int>& l_row_indices = sym_spmat.l_row_indice_;
        const std::vector<int>& u_col_compressed_slice = sym_spmat.u_col_compressed_slice_;
        const std::vector<int>& u_row_indices = sym_spmat.u_row_indice_;*/


        l_value_vec_.resize(l_col_compressed_slice[mat_size]);
        u_value_vec_.resize(u_col_compressed_slice[mat_size]);

        std::vector<T> col_vec(mat_size);
        for (int i = 0; i < mat_size; ++i) {
            for (int k = u_col_compressed_slice[i]; k < u_col_compressed_slice[i + 1] - 1; ++k) {
                col_vec[u_row_indices[k]] = 0;
            }
            for (int k = l_col_compressed_slice[i] + 1; k < l_col_compressed_slice[i + 1]; ++k) {
                col_vec[l_row_indices[k]] = 0;
            }
            for (int j = col_compressed_slice[i]; j < col_compressed_slice[i + 1]; ++j) {
                col_vec[row_indices[j]] = value_vec[j];
            }

            for (int k = u_col_compressed_slice[i]; k < u_col_compressed_slice[i + 1] - 1; ++k) {
                int col_index = u_row_indices[k];
                for (int m = l_col_compressed_slice[col_index] + 1; m < l_col_compressed_slice[col_index + 1]; ++m) {
                    int row_index = l_row_indices[m];
                    col_vec[row_index] -= l_value_vec_[m] * col_vec[col_index];
                }
            }

            for (int k = u_col_compressed_slice[i]; k < u_col_compressed_slice[i + 1]; ++k) {
                u_value_vec_[k] = col_vec[u_row_indices[k]];
            }

            for (int k = l_col_compressed_slice[i]; k < l_col_compressed_slice[i + 1]; ++k) {
                l_value_vec_[k] = col_vec[l_row_indices[k]] / col_vec[i];
            }
        }
    }

    inline void SparseLTriSolver(const int mat_size, const std::vector<int>& l_col_compressed_slice, const std::vector<int>& l_row_indices, std::vector<T>& vec_b) {
        /*const std::vector<int>& l_col_compressed_slice = sym_spmat.l_col_compressed_slice_;
        const std::vector<int>& l_row_indices = sym_spmat.l_row_indice_;

        int mat_size = l_col_compressed_slice.size() - 1;*/

        for (int i = 0; i < mat_size; ++i) {

            for (int j = l_col_compressed_slice[i] + 1; j < l_col_compressed_slice[i + 1]; ++j) {
                vec_b[l_row_indices[j]] -= l_value_vec_[j] * vec_b[i];
            }
        }
    }

    inline void SparseUTriSolver(const int mat_size, const std::vector<int>& u_col_compressed_slice, const std::vector<int>& u_row_indices, std::vector<T>& vec_b) {
        /* const std::vector<int>& u_col_compressed_slice = sym_spmat.u_col_compressed_slice_;
        const std::vector<int>& u_row_indices = sym_spmat.u_row_indice_;

        int mat_size = u_col_compressed_slice.size() - 1;*/

        for (int i = mat_size - 1; i >= 0; --i) {
            vec_b[i] /= u_value_vec_[u_col_compressed_slice[i + 1] - 1];

            for (int j = u_col_compressed_slice[i]; j < u_col_compressed_slice[i + 1] - 1; ++j) {
                vec_b[u_row_indices[j]] -= u_value_vec_[j] * vec_b[i];
            }
        }
    }
	
	// for rvalue reference
	inline void LeftLookingLU(SparseMatrix<T>&& spmat, const SymbolicSpMat& sym_spmat) {
		int&& mat_size = std::forward<int>(spmat.num_cols_);
		std::vector<T>&& value_vec = std::forward< std::vector<T> > (spmat.values_);
		std::vector<int>&& col_compressed_slice = std::forward< std::vector<int> > (spmat.col_compressed_slice_);
		std::vector<int>&& row_indices = std::forward< std::vector<int> >(spmat.row_indices_);

		const std::vector<int>& l_col_compressed_slice = sym_spmat.l_col_compressed_slice_;
		const std::vector<int>& l_row_indices = sym_spmat.l_row_indice_;
		const std::vector<int>& u_col_compressed_slice = sym_spmat.u_col_compressed_slice_;
		const std::vector<int>& u_row_indices = sym_spmat.u_row_indice_;


		l_value_vec_.resize(l_col_compressed_slice[mat_size]);
		u_value_vec_.resize(u_col_compressed_slice[mat_size]);

		std::vector<T> col_vec(mat_size);
		for (int i = 0; i < mat_size; ++i) {
			/*col_vec.resize(mat_size);
			for (int j = col_compressed_slice[i]; j < col_compressed_slice[i + 1]; ++j) {
				col_vec[row_indices[j]] = value_vec[j];
			}*/
			for (int k = u_col_compressed_slice[i]; k < u_col_compressed_slice[i + 1] - 1; ++k) {
				col_vec[u_row_indices[k]] = 0;
			}
			for (int k = l_col_compressed_slice[i] + 1; k < l_col_compressed_slice[i + 1]; ++k) {
				col_vec[l_row_indices[k]] = 0;
			}
			for (int j = col_compressed_slice[i]; j < col_compressed_slice[i + 1]; ++j) {
				col_vec[row_indices[j]] = value_vec[j];
			}

			for (int k = u_col_compressed_slice[i]; k < u_col_compressed_slice[i + 1] - 1; ++k) {
				int col_index = u_row_indices[k];
				for (int m = l_col_compressed_slice[col_index] + 1; m < l_col_compressed_slice[col_index + 1]; ++m) {
					int row_index = l_row_indices[m];
					col_vec[row_index] -= l_value_vec_[m] * col_vec[col_index];
				}
			}

			for (int k = u_col_compressed_slice[i]; k < u_col_compressed_slice[i + 1]; ++k) {
				u_value_vec_[k] = col_vec[u_row_indices[k]];
			}

			for (int k = l_col_compressed_slice[i]; k < l_col_compressed_slice[i + 1]; ++k) {
				l_value_vec_[k] = col_vec[l_row_indices[k]] / col_vec[i];
			}
		}
	}
};

#endif