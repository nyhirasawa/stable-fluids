#ifndef UPDATE_FLUID_VELOCITY_H
#define UPDATE_FLUID_VELOCITY_H
#include <vector>
#include "fftw3.h"

#include "grid.h"
#include "physical_const.h"

namespace stable_fluids{
void set_boundary_velocity(Grid& all_grid);
//外力項の計算
void add_force_fluid(Grid& all_grid);
//advect項の計算
void advect_fluid(Grid& all_grid);
//diffuse項の計算
void diffuse_fluid(Grid& all_grid);
//divergence freeな空間へのproject
void project_fluid(Grid& all_grid);
//上記4ステップをまとめただけの関数
void update_fluid_velocity(Grid& all_grid);

void diffuse_fluid_FFT(fftw_complex *velocity_in_k_space_x, fftw_complex *velocity_in_k_space_y);
void project_fluid_FFT(fftw_complex *velocity_in_k_space_x, fftw_complex *velocity_in_k_space_y);


}

#endif //UPDATE_FLUID_VELOCITY_H
