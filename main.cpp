#include <opencv2/opencv.hpp>
#include <iostream>
#include <vector>
#include "fftw3.h"

#include "grid.h"
#include "initialize_grid.h"
#include "physical_const.h"
#include "draw_substance_density.h"
#include "update_fluid_velocity.h"
#include "update_fluid_velocity_dft.h"
#include "move_substances.h"

int main(int argc, char* argv[]){
    //系の情報が乗ったグリッド
    stable_fluids::Grid all_grid(stable_fluids::physical_const::kGrid_num_x, stable_fluids::physical_const::kGrid_num_y);
    //グリッドを初期化
    stable_fluids::initialize_grid(all_grid);
    const int scale=4;

    //出力する動画のフォーマットを設定
    int fourcc, width, height;
    //フレームレート
    double fps=1.0/stable_fluids::physical_const::kDt;
    //出力形式('M', 'P', '4', 'V' ならmp4形式で出力する)
    fourcc=cv::VideoWriter::fourcc('M', 'P', '4', 'V');
    width=stable_fluids::physical_const::kGrid_num_x*scale;
    height=stable_fluids::physical_const::kGrid_num_y*scale;
    cv::VideoWriter writer;
    writer.open("result.mp4", fourcc, fps, cv::Size(width, height), false);

    const int N=stable_fluids::physical_const::kGrid_num_x*stable_fluids::physical_const::kGrid_num_y;
    double velocity_in_real_space_x[N], velocity_in_real_space_y[N];
    fftw_complex *velocity_in_k_space_x, *velocity_in_k_space_y;
	fftw_plan fft_plan_x, fft_plan_y, ifft_plan_x, ifft_plan_y;
	velocity_in_k_space_x = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) *N*N);
	velocity_in_k_space_y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) *N*N);
	fft_plan_x = fftw_plan_dft_r2c_2d(stable_fluids::physical_const::kGrid_num_x, stable_fluids::physical_const::kGrid_num_y,velocity_in_real_space_x,velocity_in_k_space_x,FFTW_ESTIMATE);
	fft_plan_y = fftw_plan_dft_r2c_2d(stable_fluids::physical_const::kGrid_num_x, stable_fluids::physical_const::kGrid_num_y,velocity_in_real_space_y,velocity_in_k_space_y,FFTW_ESTIMATE);
	ifft_plan_x = fftw_plan_dft_c2r_2d(stable_fluids::physical_const::kGrid_num_x, stable_fluids::physical_const::kGrid_num_y,velocity_in_k_space_x,velocity_in_real_space_x, FFTW_ESTIMATE);
    ifft_plan_y = fftw_plan_dft_c2r_2d(stable_fluids::physical_const::kGrid_num_x, stable_fluids::physical_const::kGrid_num_y,velocity_in_k_space_y,velocity_in_real_space_y, FFTW_ESTIMATE);

    //メインの計算部分
    const int num_frame=1000;
    for(int i=0;i<num_frame;i++){
        std::cout<<i<<"/"<<num_frame<<std::endl;
        //substance densityの出力
        stable_fluids::draw_substance_density(all_grid, scale, writer);
        //流体の速度場の更新
//        stable_fluids::update_fluid_velocity(all_grid);
        //流体の速度場の更新(離散フーリエ変換を使うバージョン(FFTを使っていないテスト用))
//        stable_fluids::update_fluid_velocity_DFT(all_grid);
        //流体の速度場の更新(離散フーリエ変換を使うバージョン(FFTを使ったバージョン))
        stable_fluids::update_fluid_velocity_FFT(all_grid, velocity_in_real_space_x, velocity_in_real_space_y, velocity_in_k_space_x, velocity_in_k_space_y, fft_plan_x, fft_plan_y, ifft_plan_x, ifft_plan_y);
        //substance densityの更新
        stable_fluids::move_substances(all_grid);
    }

    return 0;
}
