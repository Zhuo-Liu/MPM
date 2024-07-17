#ifndef GRID_H
#define GRID_H

#include "border.h"

struct Grid{
	float mi_;											
	glm::vec2 xi_;										
	glm::vec2 vi_;
	glm::vec2 vi_star_;									
	glm::vec2 vi_col_;									
	glm::vec2 vi_fri_;									
	glm::vec2 fi_;		

	void apply_collision(std::vector<Plane> planes, float dt) {
		vi_col_ = vi_star_;
		for (int b = 0; b < planes.size(); b++) {

			float phi_xn = glm::dot(planes[b].normal_, xi_) - planes[b].offset_;
			glm::vec2 x_np1 = xi_ + dt * vi_star_;
			float phi_xnp1 = glm::dot(planes[b].normal_, x_np1) - planes[b].offset_;

			// two cases:
			// if phi_xn < 0, and if phi_xnp1 < phi_xn, we need update, the velocity correction is (phi_xn - phi_xnp1)/dt
			// if phi_xn >= 0, but phi_xnp1 < 0 , we nned update, the velocity correction is (0 - phi_xnp1)/dt
			// combine the two, \delta v = (min(phi_xn,0) - phi_xnp1) / dt
			float phi_min = std::min(phi_xn, 0.0f);

			if (phi_xnp1 < phi_min){
				vi_col_ = vi_star_ + (phi_min - phi_xnp1) / dt * planes[b].normal_;
			}
		}
	}

	void apply_friction(std::vector<Plane> planes) {
		vi_fri_ = vi_col_;
		for (int b = 0; b < planes.size(); b++) {
			if (vi_col_ == vi_star_) {
				break; //early termination if no collision happens!
			}

			glm::vec2 vi_n = glm::dot(planes[b].normal_, vi_col_) * planes[b].normal_;
			glm::vec2 vi_t = vi_col_ - vi_n;
			glm::vec2 t = glm::normalize(vi_t);

			float delta_v_fric = 0.5f * glm::length(vi_col_ - vi_star_);
			if (glm::length(vi_t) > 1e-7)
			{
				if (delta_v_fric > glm::length(vi_t)) {
					vi_fri_ = vi_n;
				}
				else {
					vi_fri_ = vi_col_ - delta_v_fric * t;
				}
			}
		}
	}
};

//class of grids, intialized
class GridNode {
public:
	GridNode() {
		for (int j = 0; j < grid_num_y; j++) {
			for (int i = 0; i < grid_num_x; i++) {
				Grid grid;
				grid.mi_ = 0.0f;
				grid.vi_ = glm::vec2(0.0f);
				grid.vi_star_ = glm::vec2(0.0f);
				grid.vi_col_ = glm::vec2(0.0f);
				grid.vi_fri_ = glm::vec2(0.0f);
				grid.fi_ = glm::vec2(0.0f);
				grid.xi_ = glm::vec2((float)i, (float)j);
				grids_.push_back(grid);
			}
		}
	}
	~GridNode() {};

	std::vector<Grid> grids_;

	const int grid_num_x = 101;
	const int grid_num_y = 101;
private:
};
#endif
