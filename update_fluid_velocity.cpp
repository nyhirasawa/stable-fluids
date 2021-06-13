#include "update_fluid_velocity.h"

#include <math.h>
#include<iostream>

#include "linear_solver.h"
#include "physical_const.h"
#include "sparse_matrix.h"

namespace stable_fluids{

//壁における速度場をセットする関数
void set_boundary_velocity(Grid& all_grid){
    for(int i=1;i<stable_fluids::physical_const::kGrid_num_x-1;i++){
        all_grid.velocity[0][i][0]=-all_grid.velocity[1][i][0];
        all_grid.velocity[stable_fluids::physical_const::kGrid_num_x-1][i][0]
            =-all_grid.velocity[stable_fluids::physical_const::kGrid_num_x-2][i][0];
    }
    for(int i=1;i<stable_fluids::physical_const::kGrid_num_y-1;i++){
        all_grid.velocity[i][0][1]=-all_grid.velocity[i][1][1];
        all_grid.velocity[i][stable_fluids::physical_const::kGrid_num_y-1][1]
            =-all_grid.velocity[i][stable_fluids::physical_const::kGrid_num_y-2][1];
    }
    //四隅の速度場は周辺から線形補間
    for(int j=0;j<2;j++){
        all_grid.velocity[0][0][j]=(all_grid.velocity[1][0][j]+all_grid.velocity[0][1][j])/2.0;
        all_grid.velocity[stable_fluids::physical_const::kGrid_num_x-1][0][j]
        =(all_grid.velocity[stable_fluids::physical_const::kGrid_num_x-2][0][j]
         +all_grid.velocity[stable_fluids::physical_const::kGrid_num_x-1][1][j])/2.0;
        all_grid.velocity[0][stable_fluids::physical_const::kGrid_num_y-1][j]
        =(all_grid.velocity[0][stable_fluids::physical_const::kGrid_num_y-2][j]
         +all_grid.velocity[1][stable_fluids::physical_const::kGrid_num_y-1][j])/2.0;
        all_grid.velocity[stable_fluids::physical_const::kGrid_num_x-1][stable_fluids::physical_const::kGrid_num_y-1][j]
        =(all_grid.velocity[stable_fluids::physical_const::kGrid_num_x-2][stable_fluids::physical_const::kGrid_num_y-1][j]
         +all_grid.velocity[stable_fluids::physical_const::kGrid_num_x-1][stable_fluids::physical_const::kGrid_num_y-2][j])/2.0;
     }
}

//外力項の計算
void add_force_fluid(Grid& all_grid){
    for(int ix=0;ix<stable_fluids::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<stable_fluids::physical_const::kGrid_num_y;iy++){
            all_grid.velocity[ix][iy][1]-=9.8*stable_fluids::physical_const::kDt;
        }
    }
}

//advect項の計算
//速度場を時間 -dt だけバックトレースしてadvect項を計算する
void advect_fluid(Grid& all_grid){
    double velocity_after_advect[stable_fluids::physical_const::kGrid_num_x][stable_fluids::physical_const::kGrid_num_y][2];
    for(int ix=1;ix<stable_fluids::physical_const::kGrid_num_x-1;ix++){
        for(int iy=1;iy<stable_fluids::physical_const::kGrid_num_y-1;iy++){
            //バックトレース先の位置
            double advected_x=ix-((all_grid.velocity[ix][iy][0]*stable_fluids::physical_const::kDt)/stable_fluids::physical_const::kCell_length);
            double advected_y=iy-((all_grid.velocity[ix][iy][1]*stable_fluids::physical_const::kDt)/stable_fluids::physical_const::kCell_length);
            //バックトレース先の座標のindex
            int advected_index_x=(int)(advected_x);
            int advected_index_y=(int)(advected_y);
            //バックトレース先の座標のindexが系の外に出てしまった場合の処理
            if(advected_index_x<0){
                advected_index_x=0;
                advected_x=(advected_index_x+0.5);
            }
            if(advected_index_x>=stable_fluids::physical_const::kGrid_num_x-1){
                advected_index_x=stable_fluids::physical_const::kGrid_num_x-2;
                advected_x=(advected_index_x+0.5);
            }
            if(advected_index_y<0){
                advected_index_y=0;
                advected_y=(advected_index_y+0.5);
            }
            if(advected_index_y>=stable_fluids::physical_const::kGrid_num_y-1){
                advected_index_y=stable_fluids::physical_const::kGrid_num_y-2;
                advected_y=(advected_index_y+0.5);
            }

            //バックトレース先の速度を線形補間する
            double a0, a1, b0, b1;
            a0=advected_x-advected_index_x;
            a1=1-a0;
            b0=advected_y-advected_index_y;
            b1=1-b0;
            velocity_after_advect[ix][iy][0]=a1*(b1*all_grid.velocity[advected_index_x][advected_index_y][0]
                                            +b0*all_grid.velocity[advected_index_x][advected_index_y+1][0])
                                            +a0*(b1*all_grid.velocity[advected_index_x+1][advected_index_y][0]
                                            +b0*all_grid.velocity[advected_index_x+1][advected_index_y+1][0]);
            velocity_after_advect[ix][iy][1]=a1*(b1*all_grid.velocity[advected_index_x][advected_index_y][1]
                                            +b0*all_grid.velocity[advected_index_x][advected_index_y+1][1])
                                            +a0*(b1*all_grid.velocity[advected_index_x+1][advected_index_y][1]
                                            +b0*all_grid.velocity[advected_index_x+1][advected_index_y+1][1]);
        }
    }
    //計算結果をコピー
    for(int ix=1;ix<stable_fluids::physical_const::kGrid_num_x-1;ix++){
        for(int iy=1;iy<stable_fluids::physical_const::kGrid_num_y-1;iy++){
            all_grid.velocity[ix][iy][0]=velocity_after_advect[ix][iy][0];
            all_grid.velocity[ix][iy][1]=velocity_after_advect[ix][iy][1];
        }
    }
}


inline int get_one_d_index(int ix, int iy){
    return ix*stable_fluids::physical_const::kGrid_num_y+iy;
}
//diffuse項の計算
void diffuse_fluid(Grid& all_grid){
    //グリッドセルの総数
    const int N=stable_fluids::physical_const::kGrid_num_x
               *stable_fluids::physical_const::kGrid_num_y;

    //\nu * \Delta t /(h^2)
    double tmp_coefficient=stable_fluids::physical_const::kDynamic_viscosity_coefficient
                          *stable_fluids::physical_const::kDt
                          /(stable_fluids::physical_const::kCell_length
                           *stable_fluids::physical_const::kCell_length);
    //係数行列の計算
//    sparse_matrix A(N,N);
    sparse_matrix_with_diagonal_element A(N,N);
    for(int ix=0;ix<stable_fluids::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<stable_fluids::physical_const::kGrid_num_y;iy++){
            if(ix-1>=0){
                A.input_element(get_one_d_index(ix, iy),get_one_d_index(ix-1, iy),-tmp_coefficient);
            }
            if(iy-1>=0){
                A.input_element(get_one_d_index(ix, iy),get_one_d_index(ix, iy-1),-tmp_coefficient);
            }
            A.input_element(get_one_d_index(ix, iy),get_one_d_index(ix, iy),1+4*tmp_coefficient);
            if(iy+1<=stable_fluids::physical_const::kGrid_num_x-1){
                A.input_element(get_one_d_index(ix, iy),get_one_d_index(ix, iy+1),-tmp_coefficient);
            }
            if(ix+1<=stable_fluids::physical_const::kGrid_num_y-1){
                A.input_element(get_one_d_index(ix, iy),get_one_d_index(ix+1, iy),-tmp_coefficient);
            }
        }
    }

    //x成分とy成分について走査
    for(int icomp=0;icomp<2;icomp++){
        //連立方程式の右辺のベクトルを計算する
        std::vector<double> b(N);
        for(int ix=0;ix<stable_fluids::physical_const::kGrid_num_x;ix++){
            for(int iy=0;iy<stable_fluids::physical_const::kGrid_num_y;iy++){
                if(ix-1>=0&&ix+1<stable_fluids::physical_const::kGrid_num_x
                 &&iy-1>=0&&iy+1<stable_fluids::physical_const::kGrid_num_y){
                     b[get_one_d_index(ix, iy)]=all_grid.velocity[ix][iy][icomp];
                 }
                 else{
                     b[get_one_d_index(ix, iy)]=0.0;
                 }
            }
        }

        //係数行列がtmp_coefficientの値によって正定値にならない(?)からCG法ではなく
        //Gauss-Seidel法を使う
        ////gauss seidel法により更新後の速度場を得る
        std::vector<double> velocity_after_diffuse(N);
        gauss_seidel(A,b,velocity_after_diffuse,N,200);
        //結果をgridにコピー
        for(int ix=1;ix<stable_fluids::physical_const::kGrid_num_x-1;ix++){
            for(int iy=1;iy<stable_fluids::physical_const::kGrid_num_y-1;iy++){
                all_grid.velocity[ix][iy][icomp]=velocity_after_diffuse[get_one_d_index(ix, iy)];
            }
        }

    }
}



//divergence freeな空間へのproject
void project_fluid(Grid& all_grid){
    //グリッドの総数
    const int N=stable_fluids::physical_const::kGrid_num_x
               *stable_fluids::physical_const::kGrid_num_y;
    //係数行列の計算
    sparse_matrix A(N,N);
    for(int ix=0;ix<stable_fluids::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<stable_fluids::physical_const::kGrid_num_y;iy++){
            if(ix-1>=0){
                A.input_element(get_one_d_index(ix, iy),get_one_d_index(ix-1, iy),1);
            }
            if(iy-1>=0){
                A.input_element(get_one_d_index(ix, iy),get_one_d_index(ix, iy-1),1);
            }
            A.input_element(get_one_d_index(ix, iy),get_one_d_index(ix, iy),-4);
            if(iy+1<=stable_fluids::physical_const::kGrid_num_y-1){
                A.input_element(get_one_d_index(ix, iy),get_one_d_index(ix, iy+1),1);
            }
            if(ix+1<=stable_fluids::physical_const::kGrid_num_x-1){
                A.input_element(get_one_d_index(ix, iy),get_one_d_index(ix+1, iy),1);
            }
        }
    }
    //連立方程式の右辺のベクトルの計算
    std::vector<double> b(N);
    for(int ix=0;ix<stable_fluids::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<stable_fluids::physical_const::kGrid_num_y;iy++){
            if(ix-1>=0&&ix+1<stable_fluids::physical_const::kGrid_num_x
             &&iy-1>=0&&iy+1<stable_fluids::physical_const::kGrid_num_y){
                 b[get_one_d_index(ix, iy)]=stable_fluids::physical_const::kCell_length
                                       *(all_grid.velocity[ix+1][iy][0]-all_grid.velocity[ix-1][iy][0]
                                        +all_grid.velocity[ix][iy+1][1]-all_grid.velocity[ix][iy-1][1])
                                        /(2.0);
             }
             else{
                 b[get_one_d_index(ix, iy)]=0.0;
             }
        }
    }
    //CG法によりスカラー場 qを得る
    std::vector<double> scalar_q(N);
    conjugate_gradient(A,b,scalar_q,N,1000,0.001);

    //gauss seidel法を使う場合(Aはsparse_matrix_with_diagonal_elementにする)
//    gauss_seidel(A,b,scalar_q,N,200);

    //scalar field qからproject後の速度を計算
    for(int ix=1;ix<stable_fluids::physical_const::kGrid_num_x-1;ix++){
        for(int iy=1;iy<stable_fluids::physical_const::kGrid_num_y-1;iy++){
            all_grid.velocity[ix][iy][0]-=(scalar_q[get_one_d_index(ix+1, iy)]-scalar_q[get_one_d_index(ix-1, iy)])
            /(2.0*stable_fluids::physical_const::kCell_length);
            all_grid.velocity[ix][iy][1]-=(scalar_q[get_one_d_index(ix, iy+1)]-scalar_q[get_one_d_index(ix, iy-1)])
            /(2.0*stable_fluids::physical_const::kCell_length);
        }
    }

}

//流体の 1 time step
void update_fluid_velocity(Grid& all_grid){
//固定境界条件の場合
//    add_force_fluid(all_grid);
    set_boundary_velocity(all_grid);
    advect_fluid(all_grid);
    set_boundary_velocity(all_grid);
    diffuse_fluid(all_grid);
    set_boundary_velocity(all_grid);
    project_fluid(all_grid);
    set_boundary_velocity(all_grid);

/*
//    add_force_fluid(all_grid);
    set_boundary_velocity(all_grid);
    diffuse_fluid(all_grid);
    set_boundary_velocity(all_grid);
    project_fluid(all_grid);
    set_boundary_velocity(all_grid);
    advect_fluid(all_grid);
    set_boundary_velocity(all_grid);
    project_fluid(all_grid);
    set_boundary_velocity(all_grid);
*/
}
}
