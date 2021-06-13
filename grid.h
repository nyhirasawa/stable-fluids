#ifndef STABLE_FLUIDS_GRID_H
#define STABLE_FLUIDS_GRID_H

namespace stable_fluids{
//セルの境界条件
enum BoundaryCondition{
    FLUID=0,
    WALL=1,
};

//系の物理量が乗るグリッドの定義
class Grid{
public:
    const int Grid_num_x, Grid_num_y;
    double ***velocity;
    double **substance_density;
    BoundaryCondition **boundary_condition;
    Grid(int nx, int ny);
    ~Grid();
};
}
#endif//GRID_H
