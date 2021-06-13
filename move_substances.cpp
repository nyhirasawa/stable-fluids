#include "move_substances.h"

#include <iostream>

#include "linear_solver.h"
#include "physical_const.h"

namespace stable_fluids{

void set_periodic_boundary_condition_substance(Grid& all_grid){
    for(int i=0;i<stable_fluids::physical_const::kGrid_num_x;i++){
        double tmp=(all_grid.substance_density[0][i]+all_grid.substance_density[stable_fluids::physical_const::kGrid_num_x-1][i])/2.0;
        all_grid.substance_density[0][i]=tmp;
        all_grid.substance_density[stable_fluids::physical_const::kGrid_num_x-1][i]=tmp;
    }
    for(int i=0;i<stable_fluids::physical_const::kGrid_num_y;i++){
        double tmp=(all_grid.substance_density[i][0]+all_grid.substance_density[i][stable_fluids::physical_const::kGrid_num_y-1])/2.0;
        all_grid.substance_density[i][0]=tmp;
        all_grid.substance_density[i][stable_fluids::physical_const::kGrid_num_y-1]=tmp;
    }
    //四隅の速度場は周辺から線形補間
    all_grid.substance_density[0][0]=(all_grid.substance_density[1][0]+all_grid.substance_density[0][1])/2.0;
    all_grid.substance_density[stable_fluids::physical_const::kGrid_num_x-1][0]
    =(all_grid.substance_density[stable_fluids::physical_const::kGrid_num_x-2][0]
     +all_grid.substance_density[stable_fluids::physical_const::kGrid_num_x-1][1])/2.0;
    all_grid.substance_density[0][stable_fluids::physical_const::kGrid_num_y-1]
    =(all_grid.substance_density[0][stable_fluids::physical_const::kGrid_num_y-2]
     +all_grid.substance_density[1][stable_fluids::physical_const::kGrid_num_y-1])/2.0;
    all_grid.substance_density[stable_fluids::physical_const::kGrid_num_x-1][stable_fluids::physical_const::kGrid_num_y-1]
    =(all_grid.substance_density[stable_fluids::physical_const::kGrid_num_x-2][stable_fluids::physical_const::kGrid_num_y-1]
     +all_grid.substance_density[stable_fluids::physical_const::kGrid_num_x-1][stable_fluids::physical_const::kGrid_num_y-2])/2.0;
}


//外力項の計算
void add_force_substances(Grid& all_grid){

}

void transport_substances_DFT(Grid& all_grid){
    double substance_density_after_advect[stable_fluids::physical_const::kGrid_num_x][stable_fluids::physical_const::kGrid_num_y];
    for(int ix=0;ix<stable_fluids::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<stable_fluids::physical_const::kGrid_num_y;iy++){
            //バックトレース先の位置
            double advected_x=ix-((all_grid.velocity[ix][iy][0]*stable_fluids::physical_const::kDt)/stable_fluids::physical_const::kCell_length);
            double advected_y=iy-((all_grid.velocity[ix][iy][1]*stable_fluids::physical_const::kDt)/stable_fluids::physical_const::kCell_length);
            if (advected_x<0.5) advected_x=0.5;
            if (advected_x>stable_fluids::physical_const::kGrid_num_x-0.5) advected_x=stable_fluids::physical_const::kGrid_num_x-0.5;
            if (advected_y<0.5) advected_y=0.5;
            if (advected_y>stable_fluids::physical_const::kGrid_num_y-0.5) advected_y=stable_fluids::physical_const::kGrid_num_y-0.5;

            //バックトレース先の座標のindex
            int advected_index_x0=(int)(advected_x);
            int advected_index_y0=(int)(advected_y);
            double a0, a1, b0, b1;
            a0=advected_x-advected_index_x0;
            a1=1.0-a0;
            b0=advected_y-advected_index_y0;
            b1=1.0-b0;

            //バックトレース先の座標のindexが系の外に出てしまった場合の処理
            advected_index_x0=(advected_index_x0%stable_fluids::physical_const::kGrid_num_x+stable_fluids::physical_const::kGrid_num_x)%(stable_fluids::physical_const::kGrid_num_x);
            advected_index_y0=(advected_index_y0%stable_fluids::physical_const::kGrid_num_y+stable_fluids::physical_const::kGrid_num_y)%(stable_fluids::physical_const::kGrid_num_y);
            int advected_index_x1=(advected_index_x0+1)%stable_fluids::physical_const::kGrid_num_x;
            int advected_index_y1=(advected_index_y0+1)%stable_fluids::physical_const::kGrid_num_y;
            //バックトレース先の速度を線形補間する
            substance_density_after_advect[ix][iy]=a1*(b1*all_grid.substance_density[advected_index_x0][advected_index_y0]
                                                  +b0*all_grid.substance_density[advected_index_x0][advected_index_y1])
                                                  +a0*(b1*all_grid.substance_density[advected_index_x1][advected_index_y0]
                                                  +b0*all_grid.substance_density[advected_index_x1][advected_index_y1]);
//            線形補間しない場合
//            substance_density_after_advect[ix][iy]=all_grid.substance_density[advected_index_x][advected_index_y];

        }
    }
    //計算結果をコピー
    for(int ix=0;ix<stable_fluids::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<stable_fluids::physical_const::kGrid_num_y;iy++){
            all_grid.substance_density[ix][iy]=substance_density_after_advect[ix][iy];
        }
    }
}


//流体の速度場によってsubstanceが運ばれる項(流体のadvect項に相当)
//速度場を時間 -dt だけバックトレースして計算する
void transport_substances(Grid& all_grid){
    double substance_density_after_advect[stable_fluids::physical_const::kGrid_num_x][stable_fluids::physical_const::kGrid_num_y];
    for(int ix=1;ix<stable_fluids::physical_const::kGrid_num_x-1;ix++){
        for(int iy=1;iy<stable_fluids::physical_const::kGrid_num_y-1;iy++){
            //バックトレース先の座標
            double advected_x=ix-((all_grid.velocity[ix][iy][0]*stable_fluids::physical_const::kDt)/stable_fluids::physical_const::kCell_length);
            double advected_y=iy-((all_grid.velocity[ix][iy][1]*stable_fluids::physical_const::kDt)/stable_fluids::physical_const::kCell_length);
            //バックトレース先の座標のindex
            int advected_index_x=(int)(advected_x);
            int advected_index_y=(int)(advected_y);
            if(advected_index_x<0){
                advected_index_x=0;
                advected_x=(advected_index_x+0.5);
            }
            if(advected_index_x>=stable_fluids::physical_const::kGrid_num_x-1){
                advected_index_x=stable_fluids::physical_const::kGrid_num_x;
                advected_x=(advected_index_x+0.5);
            }
            if(advected_index_y<0){
                advected_index_y=0;
                advected_y=(advected_index_y+0.5);
            }
            if(advected_index_y>=stable_fluids::physical_const::kGrid_num_y-1){
                advected_index_y=stable_fluids::physical_const::kGrid_num_y;
                advected_y=(advected_index_y+0.5);
            }
            //バックトレース先の速度を線形補間する
            double a0, a1, b0, b1;
            a0=advected_x-advected_index_x;
            a1=1-a0;
            b0=advected_y-advected_index_y;
            b1=1-b0;
            substance_density_after_advect[ix][iy]=a1*(b1*all_grid.substance_density[advected_index_x][advected_index_y]
                                                  +b0*all_grid.substance_density[advected_index_x][advected_index_y+1])
                                                  +a0*(b1*all_grid.substance_density[advected_index_x+1][advected_index_y]
                                                  +b0*all_grid.substance_density[advected_index_x+1][advected_index_y+1]);
//            線形補間しない場合
//            substance_density_after_advect[ix][iy]=all_grid.substance_density[advected_index_x][advected_index_y];

        }
    }
    //計算結果をコピー
    for(int ix=1;ix<stable_fluids::physical_const::kGrid_num_x-1;ix++){
        for(int iy=1;iy<stable_fluids::physical_const::kGrid_num_y-1;iy++){
            all_grid.substance_density[ix][iy]=substance_density_after_advect[ix][iy];
        }
    }
}

inline int get_one_d_index(int ix, int iy){
    return ix*stable_fluids::physical_const::kGrid_num_y+iy;
}

//diffuse項の計算
void diffuse_substances(Grid& all_grid){
    //グリッドセルの総数
    const int N=stable_fluids::physical_const::kGrid_num_x
               *stable_fluids::physical_const::kGrid_num_y;

    //\kappa_a * \Delta t /(h^2)
    double tmp_coefficient=stable_fluids::physical_const::kDiffusion_const
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
            if(iy+1<=stable_fluids::physical_const::kGrid_num_y-1){
                A.input_element(get_one_d_index(ix, iy),get_one_d_index(ix, iy+1),-tmp_coefficient);
            }
            if(ix+1<=stable_fluids::physical_const::kGrid_num_x-1){
                A.input_element(get_one_d_index(ix, iy),get_one_d_index(ix+1, iy),-tmp_coefficient);
            }
        }
    }
    //連立方程式の右辺のベクトルを計算する
    std::vector<double> b(N);
    for(int ix=0;ix<stable_fluids::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<stable_fluids::physical_const::kGrid_num_y;iy++){
            if(ix-1>=0&&ix+1<stable_fluids::physical_const::kGrid_num_x
             &&iy-1>=0&&iy+1<stable_fluids::physical_const::kGrid_num_y){
                 b[get_one_d_index(ix, iy)]=all_grid.substance_density[ix][iy];
             }
             else{
                 b[get_one_d_index(ix, iy)]=0.0;
             }
        }
    }
    //係数行列がtmp_coefficientの値によって正定値にならない(?)からCG法ではなく
    //Gauss-Seidel法を使う
    ////gauss seidel法により更新後の速度場を得る
    std::vector<double> substance_density_after_diffuse(N);
//    conjugate_gradient(A,b,substance_density_after_diffuse,N,1000,0.0001);
    gauss_seidel(A,b,substance_density_after_diffuse,N,200);
    //結果をgridにコピー
    for(int ix=1;ix<stable_fluids::physical_const::kGrid_num_x-1;ix++){
        for(int iy=1;iy<stable_fluids::physical_const::kGrid_num_y-1;iy++){
            all_grid.substance_density[ix][iy]=substance_density_after_diffuse[get_one_d_index(ix, iy)];
        }
    }
}

//dissipate項の計算
void dissipate_substances(Grid& all_grid){
    for(int ix=0;ix<stable_fluids::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<stable_fluids::physical_const::kGrid_num_y;iy++){
            all_grid.substance_density[ix][iy]*=1.0/(1.0+stable_fluids::physical_const::kDissipation_rate*stable_fluids::physical_const::kDt);
        }
    }
}

//上の4ステップをまとめただけの関数(substance densityの1時間ステップ分の更新に相当)
void move_substances(Grid& all_grid){
//    add_force_substances(all_grid);
    diffuse_substances(all_grid);
    set_periodic_boundary_condition_substance(all_grid);
    transport_substances_DFT(all_grid);
    set_periodic_boundary_condition_substance(all_grid);

//    transport_substances(all_grid);
//    dissipate_substances(all_grid);
}
}
