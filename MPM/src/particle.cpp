# include "particle.h"

const float K_water = 50.0f;
const float Gamma_water = 3;

float Water::compute_force_per_particle() {
	float force = -K_water * (1.0 / pow(Jp_, Gamma_water) - 1.0) * v0_ * Jp_;
	return force;
}
