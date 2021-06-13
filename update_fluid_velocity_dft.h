#ifndef UPDATE_FLUID_VELOCITY_DFT_H
#define UPDATE_FLUID_VELOCITY_DFT_H
#include <vector>
#include <complex>
#include "fftw3.h"

#include "grid.h"
#include "physical_const.h"
#include "update_fluid_velocity.h"

namespace stable_fluids{
//周期的境界条件を課すための関数
void set_periodic_boundary_condition_fluid(Grid& all_grid);

//advect項の計算
//速度場を時間 -dt だけバックトレースしてadvect項を計算する
//周期的境界条件を適用したバージョン
void advect_fluid_PBC(Grid& all_grid);

//FFTを使ってないバージョンの離散フーリエ変換(テスト用)
void DFT(const std::vector<std::vector<std::vector<std::complex<double>>>>& input, std::vector<std::vector<std::vector<std::complex<double>>>>& output, int sign);

//波数空間におけるdifusse項の計算とdivergence freeな空間へのproject
//FFTを使ってないバージョン(テスト用)
void diffuse_and_project_fluid_DFT(std::vector<std::vector<std::vector<std::complex<double>>>>& velocity_in_k_space);

//波数空間におけるdifusse項の計算とdivergence freeな空間へのproject
//FFTを使ったバージョン
void diffuse_and_project_fluid_FFT(fftw_complex *velocity_in_k_space_x, fftw_complex *velocity_in_k_space_y);

//流体の 1 time step
//FFTを使ってないバージョン(テスト用)
void update_fluid_velocity_DFT(Grid& all_grid);

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
                               fftw_plan ifft_plan_y);
}

#endif //UPDATE_FLUID_VELOCITY_H
