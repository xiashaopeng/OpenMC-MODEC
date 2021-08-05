#ifndef SYMBOLIC_SPARSE_H_
#define SYMBOLIC_SPARSE_H_

#include "sparse_matrix.h"

class SymbolicSpMat {
public:
    std::vector<int> l_col_compressed_slice_;
    std::vector<int> l_row_indice_;
    
    std::vector<int> u_col_compressed_slice_;
    std::vector<int> u_row_indice_;

    //int mat_size_;
    int tot_nozero_ele_ = 0;

    //bool is_symbolic_factorization_ = false;

    template<typename T>
    SymbolicSpMat(const SparseMatrix<T>& spmat) {
        SymbolicFactorization(spmat);
    }
    template<typename T>
    SymbolicSpMat(SparseMatrix<T>&& spmat) {
        SymbolicFactorization(std::forward<SparseMatrix<T> >(spmat));
    }

    SymbolicSpMat() : tot_nozero_ele_(0) {};

    template<typename T>
    inline void SymbolicFactorization(const SparseMatrix<T>& spmat) {
        if (tot_nozero_ele_ == spmat.col_compressed_slice_.back()) return; 

        std::vector<std::vector<int> > u_edag(spmat.num_cols_);
        std::vector<std::vector<int> > l_edag(spmat.num_cols_);
        std::vector<int> l_support(spmat.num_cols_);

        Initial(spmat.num_cols_);
        ConstructEdags(spmat.num_cols_, spmat.col_compressed_slice_, spmat.row_indices_, u_edag, l_edag, l_support);


        tot_nozero_ele_ = spmat.col_compressed_slice_.back();
    }
	template<typename T>
	inline void SymbolicFactorization(SparseMatrix<T>&& spmat) {
        if (tot_nozero_ele_ == spmat.col_compressed_slice_.back()) return;

        std::vector<std::vector<int> > u_edag(spmat.num_cols_);
        std::vector<std::vector<int> > l_edag(spmat.num_cols_);
        std::vector<int> l_support(spmat.num_cols_);

		Initial(std::forward<int>(spmat.num_cols_));
		ConstructEdags(spmat.num_cols_, std::forward<std::vector<int> >(spmat.col_compressed_slice_), std::forward<std::vector<int> >(spmat.row_indices_), u_edag, l_edag, l_support);
        tot_nozero_ele_ = spmat.col_compressed_slice_.back();
	}

    //inline void PrintSymL(const std::string& filename) {
    //    FILE* output_file = fopen(filename.c_str(), "w");
    //    fprintf(output_file, "%d-by-%d, tot-elements: %d\n", l_col_compressed_slice_.size() - 1, l_col_compressed_slice_.size() - 1, l_col_compressed_slice_.back());
    //    for (int i = 0; i < l_col_compressed_slice_.size() - 1; ++i) {
    //        fprintf(output_file, "    col %d : locations %d to %d\n", i, l_col_compressed_slice_[i], l_col_compressed_slice_[i + 1] - 1);
    //        for (int j = l_col_compressed_slice_[i]; j < l_col_compressed_slice_[i + 1]; ++j) {
    //            fprintf(output_file, "        %d : %d\n", l_row_indice_[j], 1);
    //        }
    //    }

    //    fclose(output_file);    
    //}

    //inline void PrintSymU(const std::string& filename) {
    //    FILE* output_file = fopen(filename.c_str(), "w");
    //    fprintf(output_file, "%d-by-%d, tot-elements: %d\n", u_col_compressed_slice_.size() - 1, u_col_compressed_slice_.size() - 1, u_col_compressed_slice_.back());
    //    for (int i = 0; i < u_col_compressed_slice_.size() - 1; ++i) {
    //        fprintf(output_file, "    col %d : locations %d to %d\n", i, u_col_compressed_slice_[i], u_col_compressed_slice_[i + 1] - 1);
    //        for (int j = u_col_compressed_slice_[i]; j < u_col_compressed_slice_[i + 1]; ++j) {
    //            fprintf(output_file, "        %d : %d\n", u_row_indice_[j], 1);
    //        }
    //    }

    //    fclose(output_file);    
    //}

private:
    inline void Initial(const int& mat_size) {
        l_col_compressed_slice_.clear();
        u_col_compressed_slice_.clear();
        l_row_indice_.clear();
        u_row_indice_.clear();

        l_col_compressed_slice_.resize(mat_size + 1);
        u_col_compressed_slice_.resize(mat_size + 1);
    }
	inline void Initial(int&& mat_size) {
        l_col_compressed_slice_.clear();
        u_col_compressed_slice_.clear();
        l_row_indice_.clear();
        u_row_indice_.clear();

		l_col_compressed_slice_.resize(mat_size + 1);
		u_col_compressed_slice_.resize(mat_size + 1);
	}

    inline void ConstructEdags(const int mat_size, const std::vector<int>& col_compressed_slice, const std::vector<int>& row_indices, std::vector<std::vector<int> >& u_edag, std::vector<std::vector<int> >& l_edag, std::vector<int>& l_support) {

        for (int i = 0; i < mat_size; ++i) {
            int j = col_compressed_slice[i];
            for (; j < col_compressed_slice[i + 1]; ++j) {
                if (row_indices[j] == i) break;
            }
            //auto dis = std::lower_bound(row_indices.begin() + col_compressed_slice[i], row_indices.begin() + col_compressed_slice[i + 1], i);
            //newUEdag(row_indices.begin() + col_compressed_slice[i], dis + 0, i);
            //newLEdag(dis + 1, row_indices.begin() + col_compressed_slice[i + 1], i);
            newUEdag(mat_size, row_indices.begin() + col_compressed_slice[i], row_indices.begin() + j, i, u_edag, l_edag, l_support);
            newLEdag(mat_size, row_indices.begin() + j + 1, row_indices.begin() + col_compressed_slice[i + 1], i, u_edag, l_edag, l_support);
        }
    }

	inline void ConstructEdags(const int mat_size, std::vector<int>&& col_compressed_slice, std::vector<int>&& row_indices, std::vector<std::vector<int> >& u_edag, std::vector<std::vector<int> >& l_edag, std::vector<int>& l_support) {
		//// std::cout << "here\n";
		for (int i = 0; i < mat_size; ++i) {
			int j = col_compressed_slice[i];
			for (; j < col_compressed_slice[i + 1]; ++j) {
				if (row_indices[j] == i) break;
			}
			//auto dis = std::lower_bound(row_indices.begin() + col_compressed_slice[i], row_indices.begin() + col_compressed_slice[i + 1], i);
			//newUEdag(row_indices.begin() + col_compressed_slice[i], dis + 0, i);
			//newLEdag(dis + 1, row_indices.begin() + col_compressed_slice[i + 1], i);
			newUEdag(mat_size, row_indices.begin() + col_compressed_slice[i], row_indices.begin() + j, i, u_edag, l_edag, l_support);
			newLEdag(mat_size, row_indices.begin() + j + 1, row_indices.begin() + col_compressed_slice[i + 1], i, u_edag, l_edag, l_support);
		}
	}

    //inline void newLEdag(const std::vector<int>& Aj, const int& num_j)
    inline void newLEdag(const int mat_size, std::vector<int>::const_iterator&& Aj_begin, std::vector<int>::const_iterator&& Aj_end, const int num_j, std::vector<std::vector<int> >& u_edag, std::vector<std::vector<int> >& l_edag, std::vector<int>& l_support) {
		std::vector<int> mark_vec;
		mark_vec.resize(mat_size);

        // ����������Ԫ��ӦΪ�Խ���
        l_row_indice_.push_back(num_j);
        std::for_each(Aj_begin, Aj_end, [&mark_vec, this](const int ele) { 
            mark_vec[ele] = 1; 
            l_row_indice_.push_back(ele);
        });

        int dis = Aj_end - Aj_begin;
        l_col_compressed_slice_[num_j + 1] = l_col_compressed_slice_[num_j] + dis + 1;

        for (auto ele : u_edag[num_j]) {
            int i = l_col_compressed_slice_[ele] + 1; // �Խ���Ԫ�ز�����
            for (; i < l_col_compressed_slice_[ele + 1]; ++i) {
                if (l_row_indice_[i] > num_j) break;
            }
            for (; i < l_col_compressed_slice_[ele + 1]; ++i) {
				if (mark_vec[l_row_indice_[i]] == 0) {
					mark_vec[l_row_indice_[i]] = 1;
                    l_row_indice_.push_back(l_row_indice_[i]);
                    ++l_col_compressed_slice_[num_j + 1];
				}
            }
        }
        std::sort(l_row_indice_.begin() + l_col_compressed_slice_[num_j] + 1, l_row_indice_.begin() + l_col_compressed_slice_[num_j + 1]); 
        l_support[num_j] = 1;

        std::vector<int> row;
        for (int i = 0; i < num_j; ++i) {
            if (l_support[i] < (l_col_compressed_slice_[i + 1] - l_col_compressed_slice_[i]) && l_row_indice_[l_col_compressed_slice_[i] + l_support[i]] == num_j) {
                row.push_back(i);
                ++l_support[i];
            }
        }

		/*int mark;

        for (int i = 0; i < row.size(); ++i) {
            std::vector<int>&& tmp = DfsEdags(l_edag, {row[i]});
			mark = 1;
            for (int j = 0; j < tmp.size(); ++j) {
                if (std::binary_search(row.begin() + i + 1, row.end(), tmp[j])) {
					mark = 0;
					break;
                }
            }
			if (mark == 1)
				l_edag[row[i]].emplace_back(num_j);
        }*/

		std::vector<int> mark_row;
		mark_row.resize(row.size(), 1);
		DfsEdags(l_edag, row, mark_row);
		for (int i = 0; i < row.size(); ++i) {
			if (mark_row[i]) l_edag[row[i]].emplace_back(num_j);
		}
    }

    inline void newUEdag(const int mat_size, std::vector<int>::const_iterator&& Aj_begin, std::vector<int>::const_iterator&& Aj_end, const int num_j, std::vector<std::vector<int> >& u_edag, std::vector<std::vector<int> >& l_edag, std::vector<int>& l_support) {
		std::vector<int> mark_vec(mat_size);

        /*for (auto ele = Aj_begin; ele != Aj_end; ++ele) {
            mark_vec[*ele] = 1;
        }*/
        std::for_each(Aj_begin, Aj_end, [&mark_vec](const int ele) { mark_vec[ele] = 1; });

        //std::vector<int>&& u_col = DfsEdags(l_edag, Aj);
        //u_struct_[num_j] = u_col;

		//u_struct_[num_j] = DfsEdags(l_edag, Aj_begin, Aj_end);
        std::vector<int>&& u_col = DfsEdags(l_edag, Aj_begin, Aj_end);
        u_col_compressed_slice_[num_j + 1] = u_col_compressed_slice_[num_j] + u_col.size() + 1;
        u_row_indice_.insert(u_row_indice_.end(), std::make_move_iterator(u_col.begin()), std::make_move_iterator(u_col.end()));
        u_row_indice_.push_back(num_j);


        //int mark_size = u_col_compressed_slice_[num_j + 1] - u_col_compressed_slice_[num_j] - 1;
        //std::vector<int> mark;
        //mark.resize(mark_size, 1);

        //for (int i = mark_size - 1; i >= 0; --i) {
        //    if (mark[i] == 0) continue;
        //    u_edag[num_j].push_back(u_row_indice_[u_col_compressed_slice_[num_j] + i]);
        //    
        //    //std::vector<int>&& tmp = DfsEdags(u_edag.begin(), u_edag.begin() + num_j, { u_row_indice_[u_col_compressed_slice_[num_j] + i] });
        //    std::vector<int>&& tmp = DfsEdags(u_edag.begin(), num_j, { u_row_indice_[u_col_compressed_slice_[num_j] + i] });
        //    for(int j = 0; j < tmp.size(); ++j) {
        //        //auto dis = std::find(u_row_indice_.begin() + u_col_compressed_slice_[num_j], u_row_indice_.begin() + u_col_compressed_slice_[num_j] + i, tmp[j]);
        //        auto dis = std::lower_bound(u_row_indice_.begin() + u_col_compressed_slice_[num_j], u_row_indice_.begin() + u_col_compressed_slice_[num_j] + i, tmp[j]);
        //        if (dis != u_row_indice_.begin() + u_col_compressed_slice_[num_j] + i) 
        //            mark[dis - (u_row_indice_.begin() + u_col_compressed_slice_[num_j])] = 0;
        //    }
        //}
        //std::sort(u_edag[num_j].begin(), u_edag[num_j].end());

        int col_begin = u_col_compressed_slice_[num_j];
        int col_end = u_col_compressed_slice_[num_j + 1] - 1;
        std::vector<int> mark(num_j);
        for (int i = col_begin; i < col_end; ++i) {
            mark[u_row_indice_[i]] = 1;
        }
        for (int i = col_end - 1; i >= col_begin; --i) {
            if (mark[u_row_indice_[i]] == 0) continue;
            u_edag[num_j].push_back(u_row_indice_[i]);
            DfsEdags(u_edag.begin(), num_j, u_row_indice_[i], mark);
        }    
    }

    inline std::vector<int> DfsEdags(const std::vector<std::vector<int> >& Edags, std::vector<int>&& entry) {
        std::vector<int> elements;
        std::vector<int> mark_vec;
        mark_vec.resize(Edags.size());

        std::vector<int> mark_ele;
        mark_ele.resize(Edags.size());

        std::vector<int> traverse_elements;
        for (auto& ele : entry) {
            traverse_elements.emplace_back(ele);
            while(!traverse_elements.empty()) {
                int new_entry = traverse_elements.back();
                if (mark_ele[new_entry] == 0) {
                    mark_ele[new_entry] = 1;
                    elements.emplace_back(new_entry);
                }
                if (mark_vec[new_entry] == Edags[new_entry].size()) {
                    traverse_elements.pop_back();
                }
                else traverse_elements.emplace_back(Edags[new_entry][mark_vec[new_entry]++]);
            }
        }

        std::sort(elements.begin(), elements.end());

        return elements;
    }                                                                                     
	
	// ����0
	inline void DfsEdags(const std::vector<std::vector<int> >& Edags, std::vector<int>& entry, std::vector<int>& mark_row) {
		std::vector<int> mark_vec;
		mark_vec.resize(Edags.size());

		for (int i = 0; i < entry.size(); ++i) {
			int ele = entry[i];

			std::vector<int> traverse_elements;
			traverse_elements.emplace_back(ele);

			while (!traverse_elements.empty()) {
				int new_entry = traverse_elements.back();
				if (std::binary_search(entry.begin() + i + 1, entry.end(), new_entry)) {
					mark_row[i] = 0;
					break;
				}
				if (mark_vec[new_entry] == Edags[new_entry].size()) {
					traverse_elements.pop_back();
				}
				else traverse_elements.emplace_back(Edags[new_entry][mark_vec[new_entry]++]);
			}
		}
	}

    inline std::vector<int> DfsEdags(const std::vector<std::vector<int> >& Edags, const std::vector<int>::const_iterator& entry_begin, const std::vector<int>::const_iterator& entry_end) {
        std::vector<int> elements;
        std::vector<int> mark_vec;
        mark_vec.resize(Edags.size());

        std::vector<int> mark_ele;
        mark_ele.resize(Edags.size());

        std::vector<int> traverse_elements;
        for (auto ele = entry_begin; ele != entry_end; ++ele) {
            traverse_elements.emplace_back(*ele);
            while (!traverse_elements.empty()) {
                int new_entry = traverse_elements.back();
                if (mark_ele[new_entry] == 0) {
                    mark_ele[new_entry] = 1;
                    elements.emplace_back(new_entry);
                }
                if (mark_vec[new_entry] == Edags[new_entry].size()) {
                    traverse_elements.pop_back();
                }
                else traverse_elements.emplace_back(Edags[new_entry][mark_vec[new_entry]++]);
            }
        }

        std::sort(elements.begin(), elements.end()); 

        return elements;
    }

    inline std::vector<int> DfsEdags(std::vector<std::vector<int> >::const_iterator&& Edags_begin, const int length, std::vector<int>&& entry) {
        std::vector<int> elements;
        std::vector<int> mark_vec;
        mark_vec.resize(length);

        std::vector<int> mark_ele;
        mark_ele.resize(length);

        std::vector<int> traverse_elements;
        for (auto& ele : entry) {
            traverse_elements.emplace_back(ele);
            while (!traverse_elements.empty()) {
                int new_entry = traverse_elements.back();
                if (mark_ele[new_entry] == 0) {
                    mark_ele[new_entry] = 1;
                    elements.emplace_back(new_entry);
                }
                if (mark_vec[new_entry] == Edags_begin[new_entry].size()) {
                    traverse_elements.pop_back();
                }
                else traverse_elements.emplace_back(Edags_begin[new_entry][mark_vec[new_entry]++]);
            }
        }

        //std::sort(elements.begin(), elements.end()); 

        return elements;
    }

    inline void DfsEdags(std::vector<std::vector<int> >::const_iterator&& Edags_begin, const int& length, const int& entry, std::vector<int>& mark_row) {
        std::vector<int> mark_vec;
        mark_vec.resize(length);

        //std::stack<int> traverse_elements;
        std::vector<int> traverse_elements; // using std::vector to simulate std::stack
        //for (auto& ele : entry) {
            traverse_elements.emplace_back(entry);
            while (!traverse_elements.empty()) {
                int new_entry = traverse_elements.back();
                if (mark_row[new_entry] == 1) {
                    mark_row[new_entry] = 0;
                }
                if (mark_vec[new_entry] == Edags_begin[new_entry].size()) {
                    traverse_elements.pop_back();
                }
                else traverse_elements.emplace_back(Edags_begin[new_entry][mark_vec[new_entry]++]);
            }
        //}

    }
};

#endif