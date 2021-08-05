#ifndef C_SPARSE_H_
#define C_SPARSE_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <complex>
#include <iterator>
#include <string>
#include <fstream>
#include <map>

#include "udcmp.h"
#include "modec_global_variables.h"

template<typename T>
class SparseMatrix {
public:
    SparseMatrix() : num_rows_(MODEC::kDim), num_cols_(MODEC::kDim), num_elements_(0) {
        col_compressed_slice_.resize(MODEC::kDim + 1);

        col_indices_.reserve(MODEC::kMaxEle);
        row_indices_.reserve(MODEC::kMaxEle);
        values_.reserve(MODEC::kMaxEle);

        diagonal_index_.resize(MODEC::kDim);

        row_index_support_.resize(MODEC::kDim);
    }
    SparseMatrix(const int mat_dim) : num_rows_(mat_dim), num_cols_(mat_dim), num_elements_(0) {
        col_compressed_slice_.resize(mat_dim + 1);
        
        col_indices_.reserve(MODEC::kMaxEle);
        row_indices_.reserve(MODEC::kMaxEle);
        values_.reserve(MODEC::kMaxEle);
        
        diagonal_index_.resize(mat_dim);

        row_index_support_.resize(mat_dim);
    }
    
    SparseMatrix(const SparseMatrix& sp_mat) {
        num_cols_ = sp_mat.num_cols_;
        num_rows_ = sp_mat.num_rows_;
        col_indices_ = sp_mat.col_indices_;
        row_indices_ = sp_mat.row_indices_;
        col_compressed_slice_ = sp_mat.col_compressed_slice_;
        values_ = sp_mat.values_;
        num_elements_ = sp_mat.num_elements_;

        diagonal_index_ = sp_mat.diagonal_index_;
        index_mapper_ = sp_mat.index_mapper_;
        is_duplicated_ = sp_mat.is_duplicated_;
    }
    SparseMatrix(SparseMatrix&& sp_mat) {
        num_cols_ = sp_mat.num_cols_;
        num_rows_ = sp_mat.num_rows_;
        col_indices_ = sp_mat.col_indices_;
        row_indices_ = sp_mat.row_indices_;
        col_compressed_slice_ = sp_mat.col_compressed_slice_;
        values_ = sp_mat.values_;
        num_elements_ = sp_mat.num_elements_;
        diagonal_index_ = sp_mat.diagonal_index_;
        index_mapper_ = sp_mat.index_mapper_;
        is_duplicated_ = sp_mat.is_duplicated_;
    }

    SparseMatrix& operator=(const SparseMatrix& sp_mat) {
        if (this != &sp_mat) {
            num_cols_ = sp_mat.num_cols_;
            num_rows_ = sp_mat.num_rows_;
            col_indices_ = sp_mat.col_indices_;
            row_indices_ = sp_mat.row_indices_;
            col_compressed_slice_ = sp_mat.col_compressed_slice_;
            values_ = sp_mat.values_;
            num_elements_ = sp_mat.num_elements_;

            diagonal_index_ = sp_mat.diagonal_index_;
            index_mapper_ = sp_mat.index_mapper_;
            is_duplicated_ = sp_mat.is_duplicated_;
        }
        return *this;
    }

    SparseMatrix& operator=(SparseMatrix&& sp_mat) {
        if (this != &sp_mat) {
            num_cols_ = sp_mat.num_cols_;
            num_rows_ = sp_mat.num_rows_;
            col_indices_ = sp_mat.col_indices_;
            row_indices_ = sp_mat.row_indices_;
            col_compressed_slice_ = sp_mat.col_compressed_slice_;
            values_ = sp_mat.values_;
            num_elements_ = sp_mat.num_elements_;

            diagonal_index_ = sp_mat.diagonal_index_;
            index_mapper_ = sp_mat.index_mapper_;
            is_duplicated_ = sp_mat.is_duplicated_;
        }
        return *this;
    }

    inline void AddElement(const int i, const int j, const T value) {
        if (i < 0 || j < 0) return;
        row_indices_.emplace_back(i);
        col_indices_.emplace_back(j);
        values_.emplace_back(value);
        //row_index_support_[i].emplace_back(num_elements_);
        ++num_elements_;
    }

    inline void AddElementForCouple(const int i, const int j, const T value) {
        if (i < 0 || j < 0) return;
        row_indices_.emplace_back(i);
        col_indices_.emplace_back(j);
        values_.emplace_back(value);
        row_index_support_[i].emplace_back(num_elements_++);
        //++num_elements_;
    }
    inline void AddElementCompressed(const int i, const int j, const T value) {
        for (int k = j + 1; k <= num_rows_; ++k)
            ++col_compressed_slice_[k];
        
        row_indices_.emplace_back(i);
        values_.emplace_back(value);
    }

    inline int get_matrix_dimension() {
        return num_cols_;
    }

    inline void Resize(int new_size) {
        num_cols_ = new_size;
        num_rows_ = new_size;

        int nnz = col_compressed_slice_.back();
        col_compressed_slice_.resize(new_size + 1, nnz);

        row_index_support_.resize(new_size);

        diagonal_index_.resize(new_size);
    }

    inline T operator()(int row_index, int col_index) {
        if (num_elements_ == -1) {
            for (int i = col_compressed_slice_[col_index]; i < col_compressed_slice_[col_index + 1]; ++i) {
                if (row_indices_[i] == row_index) return values_[i];
            }
        }
        else {
            for (int i = 0; i < num_elements_; ++i) {
                if (row_indices_[i] == row_index && col_indices_[i] == col_index) return values_[i];
            }
        }
        return static_cast<T>(0.0);
    }

    // inline bool IsExisted(int row_index, int col_index, T value) {
    //     if (value == this->operator()(row_index, col_index)) return true;
    //     else return false;
    // }
    inline bool IsExisted(int row_index, int col_index, T value) {
        for(auto ele : row_index_support_[row_index]) {
            if (col_indices_[ele] == col_index && values_[ele] == value) {
                return true;
            }
        }
        return false;
    }

    inline void CompressSpMat() {
		if (num_elements_ == -1) { 
			// std::cout << "The matrix has been compressed.\n";
			return; 
		}
        std::vector<int> col_counts;
        col_counts.resize(num_cols_);
        for (int i = 0; i < num_elements_; ++i) {
            ++col_counts[col_indices_[i]];
        }
        int num_slice = 0;
        for (int i = 0; i < num_cols_; ++i) {
            col_compressed_slice_[i] = num_slice;
            num_slice += col_counts[i];
        }
        col_compressed_slice_[num_cols_] = num_slice;

        std::vector<int> tmp_vec(col_compressed_slice_.begin(), col_compressed_slice_.end() - 1);
        std::vector<int> tmp_row_vec(num_elements_);
        std::vector<T> tmp_value_vec(num_elements_);
        int tmp;
        for (int i = 0; i < num_elements_; ++i) {
            tmp_row_vec[tmp = tmp_vec[col_indices_[i]]++] = row_indices_[i];
            tmp_value_vec[tmp] = values_[i];
        }

        row_indices_ = tmp_row_vec;
        values_ = tmp_value_vec;
        num_elements_ = -1; 
    }

    inline void MappingCompressedMatrix() {
		if (num_elements_ != -1) { 
			// std::cout << "The matrix has not been compressed.\n";
			return; 
		}
        for (int i = 0; i < num_cols_; ++i) {
            std::map<int, int> row_mapper;
            int j = col_compressed_slice_[i];
            for (; j < col_compressed_slice_[i + 1]; ++j) {
                int row_index = row_indices_[j];
                row_mapper.insert(std::make_pair(row_index, j));
            }
            index_mapper_.insert(std::make_pair(i, row_mapper));
        }
    }

    inline void UpdateCompressedMatrix(const int row, const int col, const T value, const std::string& update_mode = "replace") {
        if (row < 0 || row >= num_cols_ || col < 0 || col >= num_cols_) {
            //// std::cout << "The index has exceeded the range of row or col numbers.\n";
            //// std::cout << "Row index: " << row << "; Col index: " << col << ".\n";
            return;
        }
        auto iter_row_mapper = index_mapper_[col].find(row);
        if (iter_row_mapper != index_mapper_[col].end()) {
            if (update_mode == "replace")
                values_[iter_row_mapper->second] = value;
            else if (update_mode == "increment")
                values_[iter_row_mapper->second] += value;
            else if (update_mode == "multiplication")
                values_[iter_row_mapper->second] *= value;
            else
                std::cerr << "No update_mode " << update_mode << std::endl;
        }
        else {
            // std::cout << "There is no entry at (" << row <<", " << col << ") \n";
            return; 
        }
    }

    inline void DuplicateElements() {
		if (is_duplicated_) { 
			// std::cout << "The matrix has been duplicated.\n";
			return; 
		}
        std::vector<int> tmp_vec;
        tmp_vec.resize(num_rows_, -1);
        
        int nz = 0;
        int q;
        for (int i = 0; i < num_cols_; ++i) {
            q = nz;
            for (int j = col_compressed_slice_[i]; j < col_compressed_slice_[i + 1]; ++j) {
                int tmp = row_indices_[j];

                if (tmp_vec[tmp] >= q) {
                    values_[tmp_vec[tmp]] += values_[j];
                }
                else {
                    tmp_vec[tmp] = nz;
                    row_indices_[nz] = tmp;
                    values_[nz++] = values_[j];
                }
            }
            col_compressed_slice_[i] = q;
        }
        col_compressed_slice_[num_cols_] = nz;
        values_.erase(values_.begin() + nz, values_.end());
        row_indices_.erase(row_indices_.begin() + nz, row_indices_.end());
        is_duplicated_ = true;
    }

    inline void EliminateZeroEntries() {
        for (int i = 0; i < num_cols_; ++i) {
            for (int j = col_compressed_slice_[i]; j < col_compressed_slice_[i + 1]; ++j) {
                if (values_[j] == 0.0) { //&& row_indices_[j] != i
                    values_.erase(values_.begin() + j);
                    row_indices_.erase(row_indices_.begin() + j);
                    for (int k = i + 1; k <= num_cols_; ++k) {
                        col_compressed_slice_[k] -= 1;
                    }
                }
            }
        }
    }


    inline void SortColumn() {
        /*for (int i = 0; i < num_cols_; ++i) {
            std::vector<bool> cmpres;

            std::sort(row_indices_.begin() + col_compressed_slice_[i], row_indices_.begin() + col_compressed_slice_[i + 1], udcmp<int>(cmpres));

            std::vector<bool>::const_iterator cmpit = cmpres.begin();

            std::sort(values_.begin() + col_compressed_slice_[i], values_.begin() + col_compressed_slice_[i + 1], 
            [&cmpit] (const T&, const T&) {
                return *cmpit++;
            });
        }*/

        // test for another method
        // the performance is much better than the above method (15~20x faster)
        typedef sort_helper::value_iterator_t<T, int> IndexIt;
        for (int i = 0; i < num_cols_; ++i) {
            std::sort(IndexIt(values_.begin() + col_compressed_slice_[i], row_indices_.begin() + col_compressed_slice_[i]),
                IndexIt(values_.begin() + col_compressed_slice_[i + 1], row_indices_.begin() + col_compressed_slice_[i + 1]));
        }
    }

    inline void DiagonalIndex() {
        for (int i = 0; i < num_cols_; ++i) {
            auto index = std::lower_bound(row_indices_.begin() + col_compressed_slice_[i], row_indices_.begin() + col_compressed_slice_[i + 1], i);

            diagonal_index_[i] = index - (row_indices_.begin() + col_compressed_slice_[i]);
        }
    }

    // inline void MarixLoad(std::string filename) {
    //     FILE* input_file = fopen(filename.c_str(), "r");

    //     int row, col;
    //     T value;

    //     while(fscanf(input_file, "%d %d %lg\n", &row, &col, &value) == 3) {
    //         AddElement(row, col, value);
    //     }

    //     fclose(input_file);
    // }

    // inline void PrintMatrix(std::string filename) {
    //     FILE* output_file = fopen(filename.c_str(), "w");
    //     fprintf(output_file, "%d-by-%d, tot-elements: %d\n", num_rows_, num_cols_, col_compressed_slice_[num_cols_]);
    //     for (int i = 0; i < num_cols_; ++i) {
    //         fprintf(output_file, "    col %d : locations %d to %d\n", i, col_compressed_slice_[i], col_compressed_slice_[i + 1] - 1);
    //         for (int j = col_compressed_slice_[i]; j < col_compressed_slice_[i + 1]; ++j) {
    //             fprintf(output_file, "        %d : %g\n", row_indices_[j], values_[j]);
    //         }
    //     }

    //     fclose(output_file);        
    // }

	inline void PrintMatrixStd(std::string filename) {
        std::ofstream output_file(filename);
		for (int i = 0; i < num_cols_; ++i) {
			for (int j = col_compressed_slice_[i]; j < col_compressed_slice_[i + 1]; ++j) {
                output_file << row_indices_[j] << " " << i << " " << values_[j] << "\n";
			}
		}

        output_file.close();
	}

    inline void PrintMatrixCpp(std::string filename) {
        std::ofstream output_file(filename);
        output_file << num_rows_ << "-by-" << num_cols_ << ", tot-elements: " << col_compressed_slice_[num_cols_] << '\n';

        for (int i = 0; i < num_cols_; ++i) {
            output_file << "    col " << i << " : locations " << col_compressed_slice_[i] << " to " << col_compressed_slice_[i + 1] - 1<< '\n';
            for (int j = col_compressed_slice_[i]; j < col_compressed_slice_[i + 1]; ++j) {
                //fprintf(output_file, "        %d : %g\n", row_indices_[j], values_[j]);
                output_file << "        " << row_indices_[j] << " : " << values_[j] << '\n';
            }
        }
        output_file.close();
    }

    // friend SparseMatrix operator* (const SparseMatrix& A, const T& s);

    // friend SparseMatrix operator* (const SparseMatrix& A, const SparseMatrix& B);

    // friend SparseMatrix operator+ (const SparseMatrix& A, const SparseMatrix& B);

    // friend SparseMatrix<std::complex<double> >  operator+ (const SparseMatrix& A, const std::complex<double>& s);

    inline void operator+=(const SparseMatrix& B) {
        int tmp_size = col_compressed_slice_.back() + B.col_compressed_slice_.back();
        //C.values_.resize(tmp_size);
        //C.row_indices_.resize(tmp_size);

        std::vector<int> col_compressed_slice(num_cols_ + 1);
        std::vector<int> row_indices(tmp_size);
        std::vector<T> values(tmp_size);

        std::vector<int> w(num_rows_);
        std::vector<T> x(num_rows_);

        int nz = 0;
        for (int i = 0; i < B.num_cols_; ++i) {
            col_compressed_slice[i] = nz;
            int mark = i + 1;
            for (int j = col_compressed_slice_[i]; j < col_compressed_slice_[i + 1]; ++j) {
                int tmp = row_indices_[j];
                //if (w[tmp] < mark) {
                    w[tmp] = mark;
                    row_indices[nz++] = tmp;
                    x[tmp] = values_[j];
                    //++nz;
                // }
                // else {
                //     x[tmp] += A.values_[j];
                // }
            }
            for (int j = B.col_compressed_slice_[i]; j < B.col_compressed_slice_[i + 1]; ++j) {
                int tmp = B.row_indices_[j];
                if (w[tmp] < mark) {
                    w[tmp] = mark;
                    row_indices[nz++] = tmp;
                    x[tmp] = B.values_[j];
                    //++nz;
                }
                else {
                    x[tmp] += B.values_[j];
                }
            }
            for (int j = col_compressed_slice[i]; j < nz; ++j) {
                values[j] = x[row_indices[j]];
            }
        }
        col_compressed_slice[num_cols_] = nz;
            
        values.erase(values.begin() + nz, values.end());
        row_indices.erase(row_indices.begin() + nz, row_indices.end());

        col_compressed_slice_ = std::move(col_compressed_slice);
        values_ = std::move(values);
        row_indices_ = std::move(row_indices);

        // col_compressed_slice_.swap(col_compressed_slice);
        // values_.swap(values);
        // row_indices_.swap(row_indices);
    }
    
    inline void operator+=(SparseMatrix&& B) {
        int tmp_size = col_compressed_slice_.back() + B.col_compressed_slice_.back();
        //C.values_.resize(tmp_size);
        //C.row_indices_.resize(tmp_size);

        std::vector<int> col_compressed_slice(num_cols_ + 1);
        std::vector<int> row_indices(tmp_size);
        std::vector<T> values(tmp_size);

        std::vector<int> w(num_rows_);
        std::vector<T> x(num_rows_);

        int nz = 0;
        for (int i = 0; i < B.num_cols_; ++i) {
            col_compressed_slice[i] = nz;
            int mark = i + 1;
            for (int j = col_compressed_slice_[i]; j < col_compressed_slice_[i + 1]; ++j) {
                int tmp = row_indices_[j];
                //if (w[tmp] < mark) {
                    w[tmp] = mark;
                    row_indices[nz++] = tmp;
                    x[tmp] = values_[j];
                    //++nz;
                // }
                // else {
                //     x[tmp] += A.values_[j];
                // }
            }
            for (int j = B.col_compressed_slice_[i]; j < B.col_compressed_slice_[i + 1]; ++j) {
                int tmp = B.row_indices_[j];
                if (w[tmp] < mark) {
                    w[tmp] = mark;
                    row_indices[nz++] = tmp;
                    x[tmp] = B.values_[j];
                    //++nz;
                }
                else {
                    x[tmp] += B.values_[j];
                }
            }
            for (int j = col_compressed_slice[i]; j < nz; ++j) {
                values[j] = x[row_indices[j]];
            }
        }
        col_compressed_slice[num_cols_] = nz;
            
        values.erase(values.begin() + nz, values.end());
        row_indices.erase(row_indices.begin() + nz, row_indices.end());

        col_compressed_slice_ = std::move(col_compressed_slice);
        values_ = std::move(values);
        row_indices_ = std::move(row_indices);

        // col_compressed_slice_.swap(col_compressed_slice);
        // values_.swap(values);
        // row_indices_.swap(row_indices);
    }

    inline void operator*= (const T& s) {
        for (T& ele : values_) {
            ele *= s;
        }
    }

    int num_rows_;    // number of rows
    int num_cols_;    // number of cols
    std::vector<int> col_indices_;   // col indices (size nzmax)
    std::vector<int> col_compressed_slice_;   // column pointers (sizes n + 1)
    std::vector<int> row_indices_;   // row indices, size nzmax
    std::vector<T> values_; //numerical values, size nzmax

    std::vector<std::vector<int> > row_index_support_;

    std::vector<int> diagonal_index_; 

    std::map<int, std::map<int, int> > index_mapper_;

    bool is_duplicated_ = false;

    int num_elements_; // # of entries in triplet matix, -1 for compressed-col

};

template<typename T>
inline SparseMatrix<T> operator* (const SparseMatrix<T>& A, const T& s) {
    SparseMatrix<T> B(A);

    for (int i = 0; i < B.values_.size(); ++i) {
        B.values_[i] *= s;
    }

    return B;
};

template<typename T>
inline SparseMatrix<T> operator* (const SparseMatrix<T>& A, const SparseMatrix<T>& B){
    SparseMatrix<T> C(A.num_cols_);

    int tmp_size = A.col_compressed_slice_.back() + B.col_compressed_slice_.back();
    // std::vector<T> C_values(tmp_size);
    // std::vector<int> C_row(tmp_size);

    C.values_.resize(tmp_size);
    C.row_indices_.resize(tmp_size);

    std::vector<int> w(A.num_rows_);
    std::vector<T> x(A.num_rows_);
    
    int nz = 0;
    for (int i = 0; i < B.num_cols_; ++i) {
        if (nz + A.num_rows_ > tmp_size) {
            tmp_size = 2 * tmp_size + A.num_rows_; 
            C.values_.resize(tmp_size);
            C.row_indices_.resize(tmp_size);
        }

        C.col_compressed_slice_[i] = nz;

        for (int j = B.col_compressed_slice_[i]; j < B.col_compressed_slice_[i + 1]; ++j) {
            int row_index = B.row_indices_[j];
            int mark = i + 1;
            for(int k = A.col_compressed_slice_[row_index]; k < A.col_compressed_slice_[row_index + 1]; ++k) {
                int tmp = A.row_indices_[k];
                if (w[tmp] < mark) {
                    w[tmp] = mark;
                    C.row_indices_[nz++] = tmp;
                    x[tmp] = B.values_[j] * A.values_[k];
                    //++nz;
                }
                else {
                    x[tmp] += B.values_[j] * A.values_[k];
                }
            }
        }
        
        for (int j = C.col_compressed_slice_[i]; j < nz; ++j) {
            C.values_[j] = x[C.row_indices_[j]];
        }
    }
    C.col_compressed_slice_[A.num_rows_] = nz;

    //C.values_.insert(C.values_.end(), std::make_move_iterator(C_values.begin()), std::make_move_iterator(C_values.begin() + nz)); 
    //C.row_indices_.insert(C.row_indices_.end(), std::make_move_iterator(C_row.begin()), std::make_move_iterator(C_row.begin() + nz)); 

    C.values_.erase(C.values_.begin() + nz, C.values_.end());
    C.row_indices_.erase(C.row_indices_.begin() + nz, C.row_indices_.end());

    return C;
};

template<typename T>
inline SparseMatrix<T> operator+ (const SparseMatrix<T>& A, const SparseMatrix<T>& B) {
    SparseMatrix<T> C(A.num_cols_);
    int tmp_size = A.col_compressed_slice_.back() + B.col_compressed_slice_.back();
    C.values_.resize(tmp_size);
    C.row_indices_.resize(tmp_size);

    std::vector<int> w(A.num_rows_);
    std::vector<T> x(A.num_rows_);

    int nz = 0;
    for (int i = 0; i < B.num_cols_; ++i) {
        C.col_compressed_slice_[i] = nz;
        int mark = i + 1;
        for (int j = A.col_compressed_slice_[i]; j < A.col_compressed_slice_[i + 1]; ++j) {
            int tmp = A.row_indices_[j];
            //if (w[tmp] < mark) {
                w[tmp] = mark;
                C.row_indices_[nz++] = tmp;
                x[tmp] = A.values_[j];
                //++nz;
            // }
            // else {
            //     x[tmp] += A.values_[j];
            // }
        }
        for (int j = B.col_compressed_slice_[i]; j < B.col_compressed_slice_[i + 1]; ++j) {
            int tmp = B.row_indices_[j];
            if (w[tmp] < mark) {
                w[tmp] = mark;
                C.row_indices_[nz++] = tmp;
                x[tmp] = B.values_[j];
                //++nz;
            }
            else {
                x[tmp] += B.values_[j];
            }
        }
        for (int j = C.col_compressed_slice_[i]; j < nz; ++j) {
            C.values_[j] = x[C.row_indices_[j]];
        }
    }
    C.col_compressed_slice_[A.num_cols_] = nz;
        
    C.values_.erase(C.values_.begin() + nz, C.values_.end());
    C.row_indices_.erase(C.row_indices_.begin() + nz, C.row_indices_.end());

    return C;
};

template<typename T>
inline SparseMatrix<T> operator+ (const SparseMatrix<T>& A, SparseMatrix<T>&& B) {
    SparseMatrix<T> C(A.num_cols_);
    int tmp_size = A.col_compressed_slice_.back() + B.col_compressed_slice_.back();
    C.values_.resize(tmp_size);
    C.row_indices_.resize(tmp_size);

    std::vector<int> w(A.num_rows_);
    std::vector<T> x(A.num_rows_);

    int nz = 0;
    for (int i = 0; i < B.num_cols_; ++i) {
        C.col_compressed_slice_[i] = nz;
        int mark = i + 1;
        for (int j = A.col_compressed_slice_[i]; j < A.col_compressed_slice_[i + 1]; ++j) {
            int tmp = A.row_indices_[j];
            //if (w[tmp] < mark) {
                w[tmp] = mark;
                C.row_indices_[nz++] = tmp;
                x[tmp] = A.values_[j];
                //++nz;
            // }
            // else {
            //     x[tmp] += A.values_[j];
            // }
        }
        for (int j = B.col_compressed_slice_[i]; j < B.col_compressed_slice_[i + 1]; ++j) {
            int tmp = B.row_indices_[j];
            if (w[tmp] < mark) {
                w[tmp] = mark;
                C.row_indices_[nz++] = tmp;
                x[tmp] = B.values_[j];
                //++nz;
            }
            else {
                x[tmp] += B.values_[j];
            }
        }
        for (int j = C.col_compressed_slice_[i]; j < nz; ++j) {
            C.values_[j] = x[C.row_indices_[j]];
        }
    }
    C.col_compressed_slice_[A.num_cols_] = nz;
        
    C.values_.erase(C.values_.begin() + nz, C.values_.end());
    C.row_indices_.erase(C.row_indices_.begin() + nz, C.row_indices_.end());

    return C;
};

template<typename T>
inline SparseMatrix<T> operator+ (SparseMatrix<T>&& A, SparseMatrix<T>&& B) {
    SparseMatrix<T> C(A.num_cols_);
    int tmp_size = A.col_compressed_slice_.back() + B.col_compressed_slice_.back();
    C.values_.resize(tmp_size);
    C.row_indices_.resize(tmp_size);

    std::vector<int> w(A.num_rows_);
    std::vector<T> x(A.num_rows_);

    int nz = 0;
    for (int i = 0; i < B.num_cols_; ++i) {
        C.col_compressed_slice_[i] = nz;
        int mark = i + 1;
        for (int j = A.col_compressed_slice_[i]; j < A.col_compressed_slice_[i + 1]; ++j) {
            int tmp = A.row_indices_[j];
            //if (w[tmp] < mark) {
                w[tmp] = mark;
                C.row_indices_[nz++] = tmp;
                x[tmp] = A.values_[j];
                //++nz;
            // }
            // else {
            //     x[tmp] += A.values_[j];
            // }
        }
        for (int j = B.col_compressed_slice_[i]; j < B.col_compressed_slice_[i + 1]; ++j) {
            int tmp = B.row_indices_[j];
            if (w[tmp] < mark) {
                w[tmp] = mark;
                C.row_indices_[nz++] = tmp;
                x[tmp] = B.values_[j];
                //++nz;
            }
            else {
                x[tmp] += B.values_[j];
            }
        }
        for (int j = C.col_compressed_slice_[i]; j < nz; ++j) {
            C.values_[j] = x[C.row_indices_[j]];
        }
    }
    C.col_compressed_slice_[A.num_cols_] = nz;
        
    C.values_.erase(C.values_.begin() + nz, C.values_.end());
    C.row_indices_.erase(C.row_indices_.begin() + nz, C.row_indices_.end());

    return C;
};

template<typename T>
inline std::vector<T> operator* (const SparseMatrix<T>& A, const std::vector<T>& v) {
    std::vector<T> result;
    result.resize(v.size());

    for (int i = 0; i < A.num_cols_; ++i) {
        for (int j = A.col_compressed_slice_[i]; j < A.col_compressed_slice_[i + 1]; ++j) {
            result[A.row_indices_[j]] += A.values_[j] * v[i];
        }
    }

    return result;
};

template<typename T>
inline std::vector<T> operator* (const SparseMatrix<T>& A, T* v) {
    std::vector<T> result;
    result.resize(A.num_cols_);

    for (int i = 0; i < A.num_cols_; ++i) {
        for (int j = A.col_compressed_slice_[i]; j < A.col_compressed_slice_[i + 1]; ++j) {
            result[A.row_indices_[j]] += A.values_[j] * v[i];
        }
    }

    return result;
};

template<typename T>
inline SparseMatrix<T> operator* (const std::vector<T>& u, const std::vector<T>& v) {
    SparseMatrix<T> result(u.size());

    for (int i = 0; i < u.size(); ++i) {
        if (u[i] != 0.0) {
            for (int j = 0; j < v.size(); ++j) {
                if (v[j] != 0.0) result.AddElement(i, j, u[i] * v[j]);
            }
        }
    }

    result.CompressSpMat();
    result.DuplicateElements();

    return result;
};

template<typename T>
inline const std::vector<T> operator- (const std::vector<T>& u, const std::vector<T>& v) {
    std::vector<T> result(u.size());

    for (int i = 0; i < u.size(); ++i) {
        result[i] = u[i] - v[i];
    }

    return result;
};

template<typename T>
inline const std::vector<T> operator+ (const std::vector<T>& u, const std::vector<T>& v) {
    std::vector<T> result(u.size());

    for (int i = 0; i < u.size(); ++i) {
        result[i] = u[i] + v[i];
    }

    return result;
};



template<typename T>
inline SparseMatrix<std::complex<double> > operator+ (SparseMatrix<T>&& A, const std::complex<double>& s) {
    SparseMatrix<std::complex<double> > B(A.num_cols_);

    B.col_indices_ = A.col_indices_;
    B.row_indices_ = A.row_indices_;
    B.col_compressed_slice_ = A.col_compressed_slice_;
    B.diagonal_index_ = A.diagonal_index_;
    B.num_elements_ = A.num_elements_;
    B.values_.resize(B.col_compressed_slice_.back());

    for (int i = 0; i < B.values_.size(); ++i) {
        B.values_[i] = A.values_[i];
    }
    for (int i = 0; i < B.num_cols_; ++i) {
        int index = B.col_compressed_slice_[i] + B.diagonal_index_[i];
        B.values_[index] = s + A.values_[index];
    }

    return B;
};

template<typename T>
inline SparseMatrix<std::complex<double> > operator- (SparseMatrix<T>&& A, const std::complex<double>& s) {
    SparseMatrix<std::complex<double> > B(A.num_cols_);

    B.col_indices_ = A.col_indices_;
    B.row_indices_ = A.row_indices_;
    B.col_compressed_slice_ = A.col_compressed_slice_;
    B.diagonal_index_ = A.diagonal_index_;
    B.num_elements_ = A.num_elements_;
    B.values_.resize(B.col_compressed_slice_.back());

    for (int i = 0; i < B.values_.size(); ++i) {
        B.values_[i] = A.values_[i];
    }
    for (int i = 0; i < B.num_cols_; ++i) {
        int index = B.col_compressed_slice_[i] + B.diagonal_index_[i];
        B.values_[index] = A.values_[index] - s;
    }

    return B;
};

#endif