#include "initialize_grid.h"
#include "physical_const.h"
#include <math.h>


namespace stable_fluids{
//グリッドを初期化する関数
void initialize_grid(Grid& all_grid){
    for(int ix=0;ix<stable_fluids::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<stable_fluids::physical_const::kGrid_num_y;iy++){
            //境界条件の設定
/*
            if(ix==0||ix==stable_fluids::physical_const::kGrid_num_x+1
             ||iy==0||iy==stable_fluids::physical_const::kGrid_num_y+1){
                all_grid.boundary_condition[ix][iy]=WALL;
            }
            else{
                all_grid.boundary_condition[ix][iy]=FLUID;
            }
*/

            //速度場と物質密度の初期化
            if(ix>(6*stable_fluids::physical_const::kGrid_num_x/15)
             &&ix<(9*stable_fluids::physical_const::kGrid_num_x/15)
             &&iy>(6*stable_fluids::physical_const::kGrid_num_y/15)
             &&iy<(9*stable_fluids::physical_const::kGrid_num_y/15)){
                all_grid.velocity[ix][iy][0]=10.0;
                all_grid.velocity[ix][iy][1]=0.0;
                all_grid.substance_density[ix][iy]=10.0;
            }
            else{
                all_grid.substance_density[ix][iy]=0.0;
                all_grid.velocity[ix][iy][0]=0.0;
                all_grid.velocity[ix][iy][1]=0.0;
            }


//            double x=(ix-3.5*stable_fluids::physical_const::kGrid_num_x/15.0)/*/(1.0*stable_fluids::physical_const::kGrid_num_x/15.0)*/;
//            double y=(iy-7.5*stable_fluids::physical_const::kGrid_num_y/15.0)/*/(1.0*stable_fluids::physical_const::kGrid_num_y/15.0)*/;
/*
//            all_grid.velocity[ix][iy][0]=5.0*(ix/(stable_fluids::physical_const::kGrid_num_x+2));
//            all_grid.velocity[ix][iy][1]=5.0*(iy/(stable_fluids::physical_const::kGrid_num_y+2));
            all_grid.velocity[ix][iy][0]=500.0*exp(-(x*x+y*y));
            all_grid.velocity[ix][iy][1]=0.0;

            //速度場と物質密度の初期化
            if(ix>(6*stable_fluids::physical_const::kGrid_num_x/15)
             &&ix<(9*stable_fluids::physical_const::kGrid_num_x/15)
             &&iy>(6*stable_fluids::physical_const::kGrid_num_y/15)
             &&iy<(9*stable_fluids::physical_const::kGrid_num_y/15)){
                all_grid.substance_density[ix][iy]=3.0;
            }
            else{
                all_grid.substance_density[ix][iy]=0.0;
            }
*/
        }
    }
}
}//namespace stable_fluids
