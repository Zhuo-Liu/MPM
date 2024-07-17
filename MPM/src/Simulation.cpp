#include "Simulation.h"

#include <algorithm>

//=======================Utilities=======================================//
int Simulation::get_grid_index_3d(int x, int y, int z)
{
	return x + gridnode_->grid_num_x * y + gridnode_->grid_num_x * gridnode_->grid_num_y * z;
}

//=======================End of Utilities=================================//

//Simulation::Simulation(std::vector<Water> particles, float dt){
//	borders_ = std::make_shared<Border>();
//	gridnode_ = std::make_shared<GridNode>();
//	particles_ = particles;
//
//	dt_ = dt;
//}

Simulation::Simulation(std::vector<Water> particles, float dt) {
	borders_ = std::make_shared<Border>();
	gridnode_ = std::make_shared<GridNode>();
	particles_ = particles;

	dt_ = dt;
}

void Simulation::P2G() {
#pragma omp parallel for 
	for (size_t i = 0; i < particles_.size(); i++) {
		int x0 = int(particles_[i].xp_.x);
		int y0 = int(particles_[i].xp_.y);
		int z0 = int(particles_[i].xp_.z);
		int index0 = get_grid_index_3d(x0, y0, z0);

		//water
		float force = particles_[i].compute_force_per_particle();


		//loop over active grids
		for (int z = z0 - 1; z < z0 + 3; z++) {
			for (int y = y0 - 1; y < y0 + 3; y++) {
				for (int x = x0 - 1; x < x0 + 3; x++)
				{
					int index = get_grid_index_3d(x, y, z); //grid index

					//add weight to particle for numerical purpose
					glm::vec3 xpmxi = particles_[i].xp_ - gridnode_->grids_[index].xi_;
					float wip = get_wip(xpmxi);
					glm::vec3 dwip = get_gradient_wip(xpmxi);

					//store for later use
					particles_[i].wip_list_.push_back(wip);
					particles_[i].dwip_list_.push_back(dwip);
					particles_[i].ximxp_list_.push_back(xpmxi);

					// P2G
					// 1. update grid mass
					// 2. update grid velocity
					// 3. update grid force
					float m_ip = wip * particles_[i].mp_;
					glm::mat3 Dp_inv = 3.0f * glm::mat3(1.0f);
					glm::vec3 v_ip = wip * particles_[i].mp_ * (particles_[i].vp_ + particles_[i].Bp_ * Dp_inv * (-xpmxi)); //NOTE: not divided by mass here!!!
					glm::vec3 f_ip = -force * dwip; //NO gravity !!!!

					gridnode_->grids_[index].mi_ += m_ip;
					gridnode_->grids_[index].vi_ += v_ip;
					gridnode_->grids_[index].fi_ += f_ip;

					if (m_ip > 0) {
						active_grids_.insert(index);
					}
				}
			} // end of nearby node iteration for each particle
		}
	} //end of particle iteration
}

void Simulation::update_grid()
{
	const glm::vec3 G = glm::vec3(0.0f, -10.0f, 0.0f);

	std::set<int>::iterator it = active_grids_.begin();
	// iterate till the end of set
#pragma omp parallel for schedule (dynamic)
	while (it != active_grids_.end())
	{
		int i = *it;
		Grid the_grid = gridnode_->grids_[i];
		//first complete the things left in P2G,  this can't be done in P2G
		gridnode_->grids_[i].vi_ = the_grid.vi_ / the_grid.mi_;
		gridnode_->grids_[i].fi_ = the_grid.fi_ + the_grid.mi_ * G;

		//Update grids
		//1. update v_star_
		//2. update v_col_
		//3. update v_fri_
		gridnode_->grids_[i].vi_star_ = gridnode_->grids_[i].vi_ + (gridnode_->grids_[i].fi_ / the_grid.mi_) * dt_;
		// Apply collisions and frictions
		/*TODO: change the location of collision and friction , for example, move the collision function to grid, return function!*/
		const std::vector<Plane>& planes = borders_->GetPlanes();
		gridnode_->grids_[i].apply_collision(planes, dt_);
		gridnode_->grids_[i].apply_friction(planes);

		//Increment the iterator
		it++;
	}
}

void Simulation::G2P() {
#pragma omp parallel for 
	for (int i = 0; i < particles_.size(); i++) {
		int x0 = int(particles_[i].xp_.x);
		int y0 = int(particles_[i].xp_.y);
		int z0 = int(particles_[i].xp_.z);
		int index0 = get_grid_index_3d(x0, y0, z0);

		glm::vec3 vpnp1(0.0f);
		glm::mat3 Bpnp1(0.0f);

		for (int z = z0 - 1; z < z0 + 3; z++) {
			for (int y = y0 - 1; y < y0 + 3; y++) {
				for (int x = x0 - 1; x < x0 + 3; x++)
				{
					int index = get_grid_index_3d(x, y, z); //grid index

					//Get wip and ximxp from particle
					int wip_index = (z - z0 + 1) * 16 + (y - y0 + 1) * 4 + (x - x0 + 1);
					float wip = particles_[i].wip_list_[wip_index];
					glm::vec3 xpmxi = particles_[i].ximxp_list_[wip_index];

					//G2P
					//1. update vp^n+1
					//2. update Bp^n+1
					vpnp1 += wip * gridnode_->grids_[index].vi_fri_;
					Bpnp1 += wip * glm::outerProduct(gridnode_->grids_[index].vi_fri_, -xpmxi);
				}
			}
		}

		particles_[i].vp_ = vpnp1;
		particles_[i].Bp_ = Bpnp1;
	}//end of loop of particle
}

void Simulation::update_particle()
{
#pragma omp parallel for 
	for (int i = 0; i < particles_.size(); i++) {
		int x0 = int(particles_[i].xp_.x);
		int y0 = int(particles_[i].xp_.y);
		int z0 = int(particles_[i].xp_.z);
		int index0 = get_grid_index_3d(x0, y0, z0);

		glm::vec3 Xp_buff = particles_[i].xp_;
		glm::vec3 xpnp1(0.0f);
		glm::mat3 gradient_vp(0.0f);

		for (int z = z0 - 1; z < z0 + 3; z++) {
			for (int y = y0 - 1; y < y0 + 3; y++) {
				for (int x = x0 - 1; x < x0 + 3; x++)
				{
					int index = get_grid_index_3d(x, y, z); //grid index

					//Get wip and ximxp from particle
					int wip_index = (z - z0 + 1) * 16 + (y - y0 + 1) * 4 + (x - x0 + 1);
					float wip = particles_[i].wip_list_[wip_index];
					glm::vec3 dwip = particles_[i].dwip_list_[wip_index];

					xpnp1 = xpnp1 + wip * (gridnode_->grids_[index].xi_ + dt_ * gridnode_->grids_[index].vi_col_);
					gradient_vp = gradient_vp + glm::outerProduct(gridnode_->grids_[index].vi_col_, dwip);
				}
			}
		}

		particles_[i].xp_ = xpnp1;

		// Update deformation, probably should be wrapped up!!!
		float trace_gradient_vp = gradient_vp[0][0] + gradient_vp[1][1] + +gradient_vp[2][2];
		particles_[i].Jp_ = (1.0f + dt_ * trace_gradient_vp) * particles_[i].Jp_;

	}

}

void Simulation::reset_grid() {
	int p_num = particles_.size();

	std::set<int>::iterator it = active_grids_.begin();
#pragma omp parallel for schedule (dynamic)
	while (it != active_grids_.end())
	{
		int i = *it;
		gridnode_->grids_[i].mi_ = 0;
		gridnode_->grids_[i].vi_ = glm::vec3(0.0f);
		gridnode_->grids_[i].vi_star_ = glm::vec3(0.0f);
		gridnode_->grids_[i].vi_col_ = glm::vec3(0.0f);
		gridnode_->grids_[i].vi_fri_ = glm::vec3(0.0f);
		gridnode_->grids_[i].fi_ = glm::vec3(0.0f);
		it++;
	}

#pragma omp parallel for schedule (dynamic)
	for (int j = 0; j < p_num; j++) {
		particles_[j].Resetwip();
	}

	active_grids_.clear();
}

void Simulation::forward_one_step() {
	P2G();
	update_grid();
	G2P();
	update_particle();
}


void Simulation::render_particles(Shader& shader, std::shared_ptr<Mesh> sphere, glm::mat4 projection, glm::mat4 view) {
	size_t p_num = particles_.size();

	glm::mat4* modelMatrices;
	modelMatrices = new glm::mat4[p_num];
	int index = 0;

	for (int i = 0; i < p_num; i++) {
		glm::mat4 model = glm::mat4(1.0f);
		model = glm::translate(model, 0.1f * particles_[i].xp_);
		modelMatrices[index++] = model;
	}

	unsigned int instanceVBO;
	glGenBuffers(1, &instanceVBO);
	glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(glm::mat4) * p_num, &modelMatrices[0], GL_STATIC_DRAW);

	glBindVertexArray(sphere->VAO);

	glEnableVertexAttribArray(3);
	glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), (void*)0);
	glEnableVertexAttribArray(4);
	glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), (void*)(sizeof(glm::vec4)));
	glEnableVertexAttribArray(5);
	glVertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), (void*)(2 * sizeof(glm::vec4)));
	glEnableVertexAttribArray(6);
	glVertexAttribPointer(6, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), (void*)(3 * sizeof(glm::vec4)));

	glVertexAttribDivisor(3, 1);
	glVertexAttribDivisor(4, 1);
	glVertexAttribDivisor(5, 1);
	glVertexAttribDivisor(6, 1);

	glBindVertexArray(0);

	shader.use();
	shader.setMat4("projection", projection);
	shader.setMat4("view", view);
	//glBindVertexArray(quadVAO);
	glBindVertexArray(sphere->VAO);
	//glDrawArraysInstanced(GL_TRIANGLES, 0, 6, 100); // 100 triangles of 6 vertices each
	glDrawElementsInstanced(GL_TRIANGLES, sphere->indices.size(), GL_UNSIGNED_INT, 0, p_num);
	glBindVertexArray(0);
}

void Simulation::render_borders(Shader& shader, glm::mat4 projection, glm::mat4 view) {
	//float vertices[] = {
	//	 -1.0f, 0.0f, 2.0f, // left  
	//	 2.0f, 0.0f, 2.0f, // right 
	//	 2.0f, 0.0f, -1.0f, // left  
	//	 -1.0f, 0.0f, -1.0f, // right 
	//};
	float vertices[] = {
	 3.2f, 0.2f, 6.8f, // left  
	 6.8f, 0.2f, 6.8f, // right 
	 6.8f, 0.2f, 3.2f, // left  
	 3.2f, 0.2f, 3.2f, // right 
	};

	unsigned int indices[] = {  // note that we start from 0!
	0, 1, 3,   // first triangle
	1, 2, 3    // second triangle
	};

	unsigned int VBO, VAO, EBO;
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &EBO);
	// bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
	glBindVertexArray(VAO);

	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);


	shader.use();
	shader.setMat4("projection", projection);
	shader.setMat4("view", view);
	glBindVertexArray(VAO);
	glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

	glBindVertexArray(0);
}