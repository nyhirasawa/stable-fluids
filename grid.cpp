#include "grid.h"

namespace stable_fluids{

Grid::Grid(int nx, int ny): Grid_num_x(nx), Grid_num_y(ny){
    //substance_densityのメモリ確保
    substance_density = new double*[nx];
    for(int i=0; i<nx;i++){
        substance_density[i]=new double[ny];
    }

    //BoundaryConditionのメモリ確保
    boundary_condition = new BoundaryCondition*[nx];
    for(int i=0; i<nx;i++){
        boundary_condition[i]=new BoundaryCondition[ny];
    }

    //velocityのメモリ確保
    velocity=new double**[nx];
    for(int i=0;i<nx; i++){
        velocity[i]=new double*[ny];
        for(int j=0;j<ny; j++){
            velocity[i][j]=new double[3];
        }
    }
}

Grid::~Grid(){
    //substance_densityのメモリ解放
    for(int i=0; i<Grid_num_x; i++){
        delete[] substance_density[i];
    }
    delete[] substance_density;

    //BoundaryConditionのメモリ解放
    for(int i=0; i<Grid_num_x;i++){
        delete[] boundary_condition[i];
    }
    delete[] boundary_condition;

    //velocityのメモリ解放
    for(int i=0; i<Grid_num_x; i++){
        for(int j=0; j<Grid_num_y; j++){
            delete[] velocity[i][j];
        }
        delete velocity[i];
    }
    delete[] velocity;
}

}
