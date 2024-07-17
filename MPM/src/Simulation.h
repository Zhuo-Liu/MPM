#ifndef SIMULATION_H
#define SIMULATION_H

#include<iostream>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "grid.h"
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "grid.h"
#include "particle.h"

#include "include/Shader.h"
#include "include/Mesh.h"

#include <set>

#include <set>

class Simulation {
public:
	Simulation(std::vector<Water> particles, float dt);
	~Simulation() {};

	void P2G();
	void update_grid();
	void G2P();
	void update_particle();
	void reset_grid();

	void forward_one_step();

	void render_particles(Shader& shader, std::shared_ptr<Mesh> sphere, glm::mat4 projection, glm::mat4 view);
	void render_borders(Shader& shader, glm::mat4 projection, glm::mat4 view);

	void add_particle(Water water) {
		particles_.push_back(water);
	}

	int get_grid_index_3d(int x, int y, int z);

	float get_wip(const glm::vec3& deltax_ip)
	{
		return kernel_(deltax_ip.x) * kernel_(deltax_ip.y) * kernel_(deltax_ip.z);
	}

	//equation (2) (3)
	glm::vec3 get_gradient_wip(const glm::vec3& deltax_ip)
	{
		return glm::vec3(kernel_gradient_(deltax_ip.x) * kernel_(deltax_ip.y) * kernel_(deltax_ip.z), kernel_(deltax_ip.x) * kernel_gradient_(deltax_ip.y) * kernel_(deltax_ip.z), kernel_(deltax_ip.x) * kernel_(deltax_ip.y) * kernel_gradient_(deltax_ip.z));
	}

private:
	std::vector<Water> particles_;
	std::shared_ptr<GridNode> gridnode_;
	std::shared_ptr<Border> borders_;

	std::set<int> active_grids_;

	float dt_;

	//equation (1)
	float kernel_(float deltax_ip)	//\hat{N(x)}
	{
		float wip;
		float deltax_ip_abs = fabs(deltax_ip);

		if (deltax_ip_abs < 1.0f) {
			wip = (0.5f * deltax_ip_abs * deltax_ip_abs * deltax_ip_abs - deltax_ip_abs * deltax_ip_abs + 2.0f / 3.0f);
		}
		else if (deltax_ip_abs < 2.0f) {
			wip = 1.0f / 6.0f * (2 - deltax_ip_abs) * (2 - deltax_ip_abs) * (2 - deltax_ip_abs);
		}
		else {
			wip = 0.0f;
		}

		return wip;
	}

	//equation (1)
	float kernel_gradient_(float deltax_ip)
	{
		float dwip;
		float deltax_ip_abs = fabs(deltax_ip);

		if (deltax_ip_abs < 1.0f) {
			dwip = 1.5f * deltax_ip * deltax_ip_abs - 2.0f * deltax_ip;
		}
		else if (deltax_ip_abs < 2.0f) {
			dwip = -0.5f * deltax_ip * deltax_ip_abs + 2.0f * deltax_ip - 2.0f * deltax_ip / deltax_ip_abs;
		}
		else {
			dwip = 0.0f;
		}

		return dwip;
	}
};

#endif