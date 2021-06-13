#include "update_fluid_velocity_dft.h"
#include "update_fluid_velocity.h"

#include <math.h>
#include <iostream>

namespace stable_fluids{
//周期的境界条件を課すための関数
void set_periodic_boundary_condition_fluid(Grid& all_grid){
    for(int i=0;i<stable_fluids::physical_const::kGrid_num_y;i++){
        double tmp=(all_grid.velocity[0][i][0]+all_grid.velocity[stable_fluids::physical_const::kGrid_num_x-1][i][0])/2.0;
        all_grid.velocity[0][i][0]=tmp;
        all_grid.velocity[stable_fluids::physical_const::kGrid_num_x-1][i][0]=tmp;
    }
    for(int i=0;i<stable_fluids::physical_const::kGrid_num_x;i++){
        double tmp=(all_grid.velocity[i][0][1]+all_grid.velocity[i][stable_fluids::physical_const::kGrid_num_y-1][1])/2.0;
        all_grid.velocity[i][0][1]=tmp;
        all_grid.velocity[i][stable_fluids::physical_const::kGrid_num_y-1][1]=tmp;
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

//advect項の計算
//速度場を時間 -dt だけバックトレースしてadvect項を計算する
//周期的境界条件を適用したバージョン
void advect_fluid_PBC(Grid& all_grid){
    double velocity_after_advect[stable_fluids::physical_const::kGrid_num_x][stable_fluids::physical_const::kGrid_num_y][2];
    for(int ix=0;ix<stable_fluids::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<stable_fluids::physical_const::kGrid_num_y;iy++){
            //バックトレース先の位置
            double advected_x=ix-((all_grid.velocity[ix][iy][0]*stable_fluids::physical_const::kDt)/stable_fluids::physical_const::kCell_length);
            double advected_y=iy-((all_grid.velocity[ix][iy][1]*stable_fluids::physical_const::kDt)/stable_fluids::physical_const::kCell_length);

            //バックトレース先の座標のindex
            int advected_index_x0=(int)(advected_x);
            int advected_index_y0=(int)(advected_y);
            double a0, a1, b0, b1;
            a0=advected_x-advected_index_x0;
            a1=1-a0;
            b0=advected_y-advected_index_y0;
            b1=1-b0;

            //バックトレース先の座標のindexが系の外に出てしまった場合の処理
            advected_index_x0=(advected_index_x0%stable_fluids::physical_const::kGrid_num_x+stable_fluids::physical_const::kGrid_num_x)%stable_fluids::physical_const::kGrid_num_x;
            advected_index_y0=(advected_index_y0%stable_fluids::physical_const::kGrid_num_y+stable_fluids::physical_const::kGrid_num_y)%stable_fluids::physical_const::kGrid_num_y;
            int advected_index_x1=(advected_index_x0+1)%stable_fluids::physical_const::kGrid_num_x;
            int advected_index_y1=(advected_index_y0+1)%stable_fluids::physical_const::kGrid_num_y;

            //バックトレース先の速度を線形補間する
            velocity_after_advect[ix][iy][0]=a1*(b1*all_grid.velocity[advected_index_x0][advected_index_y0][0]
                                            +b0*all_grid.velocity[advected_index_x0][advected_index_y1][0])
                                            +a0*(b1*all_grid.velocity[advected_index_x1][advected_index_y0][0]
                                            +b0*all_grid.velocity[advected_index_x1][advected_index_y1][0]);
            velocity_after_advect[ix][iy][1]=a1*(b1*all_grid.velocity[advected_index_x0][advected_index_y0][1]
                                            +b0*all_grid.velocity[advected_index_x0][advected_index_y1][1])
                                            +a0*(b1*all_grid.velocity[advected_index_x1][advected_index_y0][1]
                                            +b0*all_grid.velocity[advected_index_x1][advected_index_y1][1]);
        }
    }
    //計算結果をコピー
    for(int ix=0;ix<stable_fluids::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<stable_fluids::physical_const::kGrid_num_y;iy++){
            all_grid.velocity[ix][iy][0]=velocity_after_advect[ix][iy][0];
            all_grid.velocity[ix][iy][1]=velocity_after_advect[ix][iy][1];
        }
    }
}

//sign == 1 ... DFT
//sign == -1 ... IDFT
//FFTを使ってないバージョンの離散フーリエ変換(テスト用)
void DFT(const std::vector<std::vector<std::vector<std::complex<double>>>>& input, std::vector<std::vector<std::vector<std::complex<double>>>>& output, int sign){
    for(int ic=0;ic<2;ic++){
        for(int ix=0;ix<(stable_fluids::physical_const::kGrid_num_x);ix++){
            for(int iy=0;iy<(stable_fluids::physical_const::kGrid_num_y);iy++){
                output[ix][iy][ic]=std::complex<double>(0.0, 0.0);
            }
        }
    }
    for(int ic=0;ic<2;ic++){
        for(int ix0=0;ix0<(stable_fluids::physical_const::kGrid_num_x);ix0++){
            for(int iy0=0;iy0<(stable_fluids::physical_const::kGrid_num_y);iy0++){
                for(int ix1=0;ix1<(stable_fluids::physical_const::kGrid_num_x);ix1++){
                    for(int iy1=0;iy1<(stable_fluids::physical_const::kGrid_num_y);iy1++){
                        output[ix0][iy0][ic]+=input[ix1][iy1][ic]
                                             *std::exp(2.0*M_PI*std::complex<double>(0, -sign)
                                             *( (double)ix0*(double)ix1/(double)(stable_fluids::physical_const::kGrid_num_x)
                                             + (double)iy0*(double)iy1/(double)(stable_fluids::physical_const::kGrid_num_y) ));
                    }
                }
            }
        }
    }
}

//波数空間におけるdifusse項の計算とdivergence freeな空間へのproject
//FFTを使ってないバージョン(テスト用)
void diffuse_and_project_fluid_DFT(std::vector<std::vector<std::vector<std::complex<double>>>>& velocity_in_k_space){
    for(int ix=0;ix<stable_fluids::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<stable_fluids::physical_const::kGrid_num_y;iy++){
            //波数ベクトル
            double k[2];
            if(ix<=stable_fluids::physical_const::kGrid_num_x/2){
                k[0]=ix;
            }
            else{
                k[0]=ix-stable_fluids::physical_const::kGrid_num_x;
            }
            if(iy<=stable_fluids::physical_const::kGrid_num_y/2){
                k[1]=iy;
            }
            else{
                k[1]=iy-stable_fluids::physical_const::kGrid_num_y;
            }
            double k_pow2 = k[0] * k[0] + k[1] * k[1];
            //ゼロ除算を回避するための処置
    		if (k_pow2>0.00001){
                //波数空間におけるdiffuse項
                double diffuse_factor=1.0/(1.0+stable_fluids::physical_const::kDynamic_viscosity_coefficient
                                      *stable_fluids::physical_const::kDt
                                      *k_pow2);
//                double diffuse_factor = exp(-k_pow2*stable_fluids::physical_const::kDt*stable_fluids::physical_const::kDynamic_viscosity_coefficient);
                std::complex<double> velocity_in_k_space_x_previous=velocity_in_k_space[ix][iy][0];
                std::complex<double> velocity_in_k_space_y_previous=velocity_in_k_space[ix][iy][1];
                //diffuse項とprojectを一緒に計算する
                velocity_in_k_space[ix][iy][0]
                    =std::complex<double>(diffuse_factor*((1-k[0]*k[0]/k_pow2)*velocity_in_k_space_x_previous.real()-(k[0]*k[1]/k_pow2)*velocity_in_k_space_y_previous.real())
                                         ,diffuse_factor*((1-k[0]*k[0]/k_pow2)*velocity_in_k_space_x_previous.imag()-(k[0]*k[1]/k_pow2)*velocity_in_k_space_y_previous.imag()));
                velocity_in_k_space[ix][iy][1]
                    =std::complex<double>(diffuse_factor*((-k[0]*k[1]/k_pow2)*velocity_in_k_space_x_previous.real()+(1-k[1]*k[1]/k_pow2)*velocity_in_k_space_y_previous.real())
                                         ,diffuse_factor*((-k[0]*k[1]/k_pow2)*velocity_in_k_space_x_previous.imag()+(1-k[1]*k[1]/k_pow2)*velocity_in_k_space_y_previous.imag()));
            }
        }
    }
}

//波数空間におけるdifusse項の計算とdivergence freeな空間へのproject
void diffuse_and_project_fluid_FFT(fftw_complex *velocity_in_k_space_x, fftw_complex *velocity_in_k_space_y){
    for(int ix=0;ix<stable_fluids::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<(stable_fluids::physical_const::kGrid_num_y/2)+1;iy++){
            //波数ベクトル
            double k[2];
            if(ix<=stable_fluids::physical_const::kGrid_num_x/2){
                k[0]=ix;
            }
            else{
                k[0]=ix-stable_fluids::physical_const::kGrid_num_x;
            }
            k[1]=iy;
            double k_pow2 = k[0] * k[0] + k[1] * k[1];
            //ゼロ除算を回避するための処置
    		if (k_pow2>0.00001){
                //波数空間におけるdiffuse項
                double diffuse_factor=1.0/(1.0+stable_fluids::physical_const::kDynamic_viscosity_coefficient
                                      *stable_fluids::physical_const::kDt
                                      *k_pow2);
//                double diffuse_factor = exp(-k_pow2*stable_fluids::physical_const::kDt*stable_fluids::physical_const::kDynamic_viscosity_coefficient);
                double velocity_in_k_space_x_previous[2], velocity_in_k_space_y_previous[2];
                velocity_in_k_space_x_previous[0]=velocity_in_k_space_x[iy+((stable_fluids::physical_const::kGrid_num_y/2)+1)*ix][0];
                velocity_in_k_space_x_previous[1]=velocity_in_k_space_x[iy+((stable_fluids::physical_const::kGrid_num_y/2)+1)*ix][1];
                velocity_in_k_space_y_previous[0]=velocity_in_k_space_y[iy+((stable_fluids::physical_const::kGrid_num_y/2)+1)*ix][0];
                velocity_in_k_space_y_previous[1]=velocity_in_k_space_y[iy+((stable_fluids::physical_const::kGrid_num_y/2)+1)*ix][1];
                //diffuse項とprojectを一緒に計算する
                velocity_in_k_space_x[iy+((stable_fluids::physical_const::kGrid_num_y/2)+1)*ix][0]
                    =diffuse_factor*((1-k[0]*k[0]/k_pow2)*velocity_in_k_space_x_previous[0]-(k[0]*k[1]/k_pow2)*velocity_in_k_space_y_previous[0]);
                velocity_in_k_space_x[iy+((stable_fluids::physical_const::kGrid_num_y/2)+1)*ix][1]
                    =diffuse_factor*((1-k[0]*k[0]/k_pow2)*velocity_in_k_space_x_previous[1]-(k[0]*k[1]/k_pow2)*velocity_in_k_space_y_previous[1]);
                velocity_in_k_space_y[iy+((stable_fluids::physical_const::kGrid_num_y/2)+1)*ix][0]
                    =diffuse_factor*((-k[0]*k[1]/k_pow2)*velocity_in_k_space_x_previous[0]+(1-k[1]*k[1]/k_pow2)*velocity_in_k_space_y_previous[0]);
                velocity_in_k_space_y[iy+((stable_fluids::physical_const::kGrid_num_y/2)+1)*ix][1]
                    =diffuse_factor*((-k[0]*k[1]/k_pow2)*velocity_in_k_space_x_previous[1]+(1-k[1]*k[1]/k_pow2)*velocity_in_k_space_y_previous[1]);
            }
        }
    }
}

//流体の 1 time step
//FFTを使ってないバージョン(テスト用)
void update_fluid_velocity_DFT(Grid& all_grid){
    //グリッドの総数
    const int N=(stable_fluids::physical_const::kGrid_num_x)
               *(stable_fluids::physical_const::kGrid_num_y);

//周期境界条件の場合
//    add_force_fluid(all_grid);
    advect_fluid_PBC(all_grid);
    set_periodic_boundary_condition_fluid(all_grid);

    std::vector<std::vector<std::vector<std::complex<double>>>> velocity_in_real_space(
        stable_fluids::physical_const::kGrid_num_x,std::vector<std::vector<std::complex<double>>>(
        stable_fluids::physical_const::kGrid_num_y,std::vector<std::complex<double>>(2)));
    //波数空間での速度
    std::vector<std::vector<std::vector<std::complex<double>>>> velocity_in_k_space(
        stable_fluids::physical_const::kGrid_num_x,std::vector<std::vector<std::complex<double>>>(
        stable_fluids::physical_const::kGrid_num_y,std::vector<std::complex<double>>(2)));
    for(int ic=0;ic<2;ic++){
        for(int ix=0;ix<stable_fluids::physical_const::kGrid_num_x;ix++){
            for(int iy=0;iy<stable_fluids::physical_const::kGrid_num_y;iy++){
                velocity_in_real_space[ix][iy][ic]=std::complex<double>(all_grid.velocity[ix][iy][ic], 0.0);
                velocity_in_k_space[ix][iy][ic]=std::complex<double>(0.0, 0.0);
            }
        }
    }
    //離散フーリエ変換
    DFT(velocity_in_real_space, velocity_in_k_space, 1);
    //波数空間におけるdifusse項の計算とdivergence freeな空間へのproject
    diffuse_and_project_fluid_DFT(velocity_in_k_space);
    //逆離散フーリエ変換
    DFT(velocity_in_k_space, velocity_in_real_space, -1);
    //正規化係数を掛けて結果を格納
    for(int ic=0;ic<2;ic++){
        for(int ix=0;ix<stable_fluids::physical_const::kGrid_num_x;ix++){
            for(int iy=0;iy<stable_fluids::physical_const::kGrid_num_y;iy++){
                all_grid.velocity[ix][iy][ic]=velocity_in_real_space[ix][iy][ic].real()/(N);
            }
        }
    }
    set_periodic_boundary_condition_fluid(all_grid);
}

inline int get_one_d_index(int ix, int iy){
    return ix*stable_fluids::physical_const::kGrid_num_y+iy;
}

//流体の 1 time step
//FFTを使ったバージョン
void update_fluid_velocity_FFT(Grid&all_grid,
                               double *velocity_in_real_space_x,
                               double *velocity_in_real_space_y,
                               fftw_complex *velocity_in_k_space_x,
                               fftw_complex *velocity_in_k_space_y,
                               fftw_plan fft_plan_x,
                               fftw_plan fft_plan_y,
                               fftw_plan ifft_plan_x,
                               fftw_plan ifft_plan_y){

//    add_force_fluid(all_grid);
    //advect項の計算
    advect_fluid_PBC(all_grid);
    set_periodic_boundary_condition_fluid(all_grid);
    for(int ix=0;ix<stable_fluids::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<stable_fluids::physical_const::kGrid_num_y;iy++){
            velocity_in_real_space_x[iy+stable_fluids::physical_const::kGrid_num_y*ix]=all_grid.velocity[ix][iy][0];
            velocity_in_real_space_y[iy+stable_fluids::physical_const::kGrid_num_y*ix]=all_grid.velocity[ix][iy][1];
        }
    }
    //速度場の離散フーリエ変換
    fftw_execute(fft_plan_x);
    fftw_execute(fft_plan_y);
    //波数空間におけるdifusse項の計算とdivergence freeな空間へのproject
    diffuse_and_project_fluid_FFT(velocity_in_k_space_x, velocity_in_k_space_y);
    //速度場の逆離散フーリエ変換
    fftw_execute(ifft_plan_x);
    fftw_execute(ifft_plan_y);
    //正規化係数を掛けて結果を格納
    double normalize_factor = 1.0/(stable_fluids::physical_const::kGrid_num_x*stable_fluids::physical_const::kGrid_num_y);
    for (int ix=0;ix<stable_fluids::physical_const::kGrid_num_x;ix++){
        for (int iy=0;iy<stable_fluids::physical_const::kGrid_num_y;iy++){
            all_grid.velocity[ix][iy][0]=velocity_in_real_space_x[iy+stable_fluids::physical_const::kGrid_num_y*ix]*normalize_factor;
			all_grid.velocity[ix][iy][1]=velocity_in_real_space_y[iy+stable_fluids::physical_const::kGrid_num_y*ix]*normalize_factor;
		}
	}
    set_periodic_boundary_condition_fluid(all_grid);

}

}
