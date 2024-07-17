#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>

//class Particle
//{
//public:
//	glm::vec2 xp_;
//	glm::vec2 vp_;
//	glm::mat2 Bp_;
//
//	float mp_;
//	float v0_;
//
//	//void UpdateDeformation(); //made it virtual, realize it later
//	//VPOD compute_force_per_particle(); // make it virtual
//};

const float K_water = 50.0f;
const float Gamma_water = 3.0f;

/* Dry Sand */
const float RHO_dry_sand = 1600.0;				// Density
const float E_dry_sand = 3.537e5;				// Young's modulus
const float V_dry_sand = 0.3;					// Poisson's ratio
const float LAMBDA_dry_sand =					// Lame parameters
E_dry_sand * V_dry_sand / (1.0 + V_dry_sand) / (1.0 - 2.0 * V_dry_sand);
const float MU_dry_sand = E_dry_sand / (1.0 + V_dry_sand) / 2.0;

const float PI = 3.14159265358979323846;
const float H0 = 35 * PI / 180.0;				// Hardening parameters
const float H1 = 9 * PI / 180.0;
const float H2 = 0.2;
const float H3 = 10 * PI / 180.0;

class Water {
public:
	Water(glm::vec2 x_p, glm::vec2 v_p, glm::mat2 B_p, float mass, float v0) {
		xp_ = x_p;
		vp_ = v_p;
		Bp_ = B_p;
		mp_ = mass;
		v0_ = v0;

		Jp_ = 1.0f;
	}

	float compute_force_per_particle();

	float Jp_; //Deformation

	glm::vec2 xp_; //particle location
	glm::vec2 vp_; //particle velocity
	glm::mat2 Bp_; //particle velocity field

	float mp_; //particle mass
	float v0_; //particle initial volumn

	//calculation purpose
	std::vector<float> wip_list_;
	std::vector<glm::vec2> dwip_list_;
	std::vector<glm::vec2> ximxp_list_;

	void Resetwip() {
		wip_list_.clear();
		dwip_list_.clear();
		ximxp_list_.clear();
	}
private:

};


class Sand {
public:
	Sand(glm::vec2 x_p, glm::vec2 v_p, glm::mat2 B_p, float mass, float v0) {
		xp_ = x_p;
		vp_ = v_p;
		Bp_ = B_p;
		mp_ = mass;
		v0_ = v0;

		Ap_ = glm::mat2(0.0f);
		Fe_ = glm::mat2(1.0f);
		FeTr_ = glm::mat2(1.0f);
		Fp_ = glm::mat2(1.0f);
		FpTr_ = glm::mat2(1.0f);
		q_ = 0.0;

		float phi = H0 + (H1 * q_ - H3) * exp(-H2 * q_);
		alpha_ = sqrt(2.0 / 3.0) * 2 * sin(phi) / (3.0 - sin(phi));
	}

	glm::mat2 compute_force_per_particle();
	void UpdateDeformation(glm::mat2& T, float dt);

	glm::mat2 Fe_;
	glm::mat2 FeTr_;
	glm::mat2 Fp_;
	glm::mat2 FpTr_;
	glm::mat2 Ap_;
	float q_;
	float alpha_;

	glm::vec2 xp_; //particle location
	glm::vec2 vp_; //particle velocity
	glm::mat2 Bp_; //particle velocity field

	float mp_; //particle mass
	float v0_; //particle initial volumn

	//calculation purpose
	std::vector<float> wip_list_;
	std::vector<glm::vec2> dwip_list_;
	std::vector<glm::vec2> ximxp_list_;

	void Resetwip() {
		wip_list_.clear();
		dwip_list_.clear();
		ximxp_list_.clear();
	}
private:
};
#endif