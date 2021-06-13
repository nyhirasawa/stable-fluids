#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H
#include <iostream>
#include <vector>

namespace stable_fluids{
//CSR形式の疎行列
class sparse_matrix{
public:
    int num_element, num_row, num_column;
    int input_row_num;
    //非ゼロ成分を記録する
    std::vector<double> element_value;
    //非ゼロ成分が何列目にあるかを保持する
    std::vector<int> column_index;
    //各行に何個非ゼロ成分があるかを保持する
    std::vector<int> num_nonzero_element_in_row;

    sparse_matrix(int nr, int nc);
    //row行目、column列目の値elementを入力する関数
    void input_element(int row, int column, double element);
    //dense形式で行列を表示する関数
    void print_matrix();
    void print_information();
};

//CSR形式の疎行列クラスsparse_matrixに加えて対角成分の値を専用のベクターに保持するクラス
//(sparse_matrixクラスだけで疎行列の情報としては完全なのだが、対角成分の値を保持したベクターを
//別で持っておくとGauss-Seidelのアルゴリズムで便利なのでこのクラスを定義した)
class sparse_matrix_with_diagonal_element : public sparse_matrix{
public:
    std::vector<double> diagonal_element_value;

    sparse_matrix_with_diagonal_element(int nr, int nc);
    void input_element(int row, int column, double element);
};

void mat_vec_product(const sparse_matrix A,const std::vector<double> b, std::vector<double>& x, int n);
}//namespace stable_fluids
#endif
