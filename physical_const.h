#ifndef PHYSICAL_CONST_H
#define PHYSICAL_CONST_H

namespace stable_fluids{
struct physical_const{
//x, y方向のグリッドセルの個数(壁は含まない)
static const int kGrid_num_x=128;
static const int kGrid_num_y=128;
//static const int kGrid_num=(kGrid_num_x+2)*(kGrid_num_y+2);
//セルの一辺の長さ
static constexpr double kCell_length=1.0/kGrid_num_x;
//1時間ステップの長さ
static constexpr double kDt=0.01;
//流体の動的粘性係数
static constexpr double kDynamic_viscosity_coefficient=0.001;
//物質の拡散係数
static constexpr double kDiffusion_const=0.001;
//物質の散逸係数
static constexpr double kDissipation_rate=0.0;
};
} // namespace stable_fluids

#endif
