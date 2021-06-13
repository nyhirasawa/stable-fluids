#ifndef DRAW_SUBSTANCE_DENSITY
#define DRAW_SUBSTANCE_DENSITY
#include <vector>
#include <opencv2/opencv.hpp>
#include "grid.h"
#include "physical_const.h"


namespace stable_fluids{
//物質密度を描画する関数
void draw_substance_density(Grid& all_grid, int scale, cv::VideoWriter& writer);
}

#endif //DRAW_SUBSTANCE_DENSITY
