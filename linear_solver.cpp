#include "linear_solver.h"
#include <iostream>
namespace stable_fluids{
//CG法Denseバージョン
void conjugate_gradient(const std::vector<std::vector<double>> &A, const std::vector<double> &b, std::vector<double> &x, int n){
    std::vector<double> error(n), direction(n), Ap(n);
    for(int i=0;i<n; i++){
        x[i]=0.0;
    }
    //初期error
    for(int ix=0;ix<n;ix++){
        double Ax_tmp=0.0;
        for(int iy=0;iy<n;iy++){
            Ax_tmp+=A[ix][iy]*x[iy];
        }
        error[ix]=b[ix]-Ax_tmp;
        direction[ix]=error[ix];
    }
    //メインの計算部分
    for(int i=0;i<100000;i++){
        for(int ix=0;ix<n;ix++){
            Ap[ix]=0.0;
            for(int iy=0;iy<n;iy++){
                Ap[ix]+=A[ix][iy]*direction[iy];
            }
        }
        //修正係数coeffの計算
        double errerr_prev=0.0, diry=0.0;
        for(int ix=0;ix<n;ix++){
            errerr_prev+=error[ix]*error[ix];
            diry+=Ap[ix]*direction[ix];
        }
        double coeff_0=errerr_prev/diry;
        for(int ix=0;ix<n;ix++){
            //解の近似値を計算
            x[ix]+=coeff_0*direction[ix];
            //errorの計算
            error[ix]-=coeff_0*Ap[ix];
        }
        //directionの計算
        double errerr_next=0.0;
        for(int ix=0;ix<n;ix++){
            errerr_next+=error[ix]*error[ix];
        }
        //終了条件
        if(errerr_next<0.01*0.01){
            return;
        }
        //directionの修正係数coeff_0の計算
        double coeff_dir=errerr_next/errerr_prev;
        for(int ix=0;ix<n;ix++){
            direction[ix]=error[ix]+coeff_dir*direction[ix];
        }
    }
}

//CG法sparseバージョン
void conjugate_gradient(const sparse_matrix &A, const std::vector<double> &b, std::vector<double> &x, int n, int max_itr, double eps){
    std::vector<double> error(n), direction(n), Adir(n);

    for(int i=0;i<n; i++){
        x[i]=0.0;
    }
    //初期error
    std::vector<double> Ax_tmp(n);
    mat_vec_product(A, x, Ax_tmp, n);
    for(int ix=0;ix<n;ix++){
        error[ix]=b[ix]-Ax_tmp[ix];
        direction[ix]=error[ix];
    }
    //初期エラーの大きさ
    double errerr_0=0.0;
    for(int ix=0;ix<n;ix++){
        errerr_0+=error[ix]*error[ix];
    }
    //メインの計算部分
    for(int i=0;i<max_itr;i++){
        mat_vec_product(A, direction, Adir, n);

        //解の修正係数coeff_0の計算
        double errerr_prev=0.0, diry=0.0;
        for(int ix=0;ix<n;ix++){
            errerr_prev+=error[ix]*error[ix];
            diry+=Adir[ix]*direction[ix];
        }
        double coeff_0=errerr_prev/diry;
        for(int ix=0;ix<n;ix++){
            //解の近似値を計算
            x[ix]+=coeff_0*direction[ix];
            //errorの計算
            error[ix]-=coeff_0*Adir[ix];
        }
        //directionの計算
        double errerr_next=0.0;
        for(int ix=0;ix<n;ix++){
            errerr_next+=error[ix]*error[ix];
        }
        //終了条件
        if((errerr_next/errerr_0)<eps*eps){
            return;
        }
/*
        if(errerr_next<eps*eps){
            return;
        }
*/
        //directionの修正係数coeff_0の計算
        double coeff_dir=errerr_next/errerr_prev;
        for(int ix=0;ix<n;ix++){
            direction[ix]=error[ix]+coeff_dir*direction[ix];
        }
    }
}

//gauss_seidel法sparseバージョン
void gauss_seidel(const sparse_matrix_with_diagonal_element &A, const std::vector<double> &b, std::vector<double> &x, int N, int max_itr){
    for(int ix=0;ix<N;ix++){
        x[ix]=0.0;
    }
    for(int k=0;k<max_itr;k++){
        int irow=0;
        int ielem=0;
        for(auto itr=A.num_nonzero_element_in_row.begin();itr!=A.num_nonzero_element_in_row.end();itr++){
            double Ax_tmp=0.0;
            for(int i=0;i<*itr;i++){
                if(irow!=A.column_index[ielem]){
                    Ax_tmp+=A.element_value[ielem]*x[A.column_index[ielem]];
                }
                ielem++;
            }
            x[irow]=(b[irow]-Ax_tmp)/A.diagonal_element_value[irow];
            irow++;
        }
    }
}

}//namespace stable_fluids
