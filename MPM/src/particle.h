#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

class Water {
public:
	Water(glm::vec3 x_p, glm::vec3 v_p, glm::mat3 B_p, float mass, float v0) {
		xp_ = x_p;
		vp_ = v_p;
		Bp_ = B_p;
		mp_ = mass;
		v0_ = v0;

		Jp_ = 1.0f;
	}

	float compute_force_per_particle();

	float Jp_; //Deformation

	glm::vec3 xp_; //particle location
	glm::vec3 vp_; //particle velocity
	glm::mat3 Bp_; //particle velocity field

	float mp_; //particle mass
	float v0_; //particle initial volumn

	//calculation purpose
	std::vector<float> wip_list_;
	std::vector<glm::vec3> dwip_list_;
	std::vector<glm::vec3> ximxp_list_;

	void Resetwip() {
		wip_list_.clear();
		dwip_list_.clear();
		ximxp_list_.clear();
	}
private:

};

#endif