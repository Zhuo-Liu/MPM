# include "particle.h"
#include <math.h>  

float Water::compute_force_per_particle() {
	float force =  - K_water * (1.0 / pow(Jp_, Gamma_water) - 1.0) * v0_ * Jp_;
	return force;
}

glm::mat2 Sand::compute_force_per_particle() {

	//Performing SVD for Fe_
    Eigen::Matrix2f U, V, Eps, cov;
    Eigen::Vector2f w;
    Eigen::JacobiSVD<Eigen::Matrix2f> svd;
    cov << Fe_[0][0], Fe_[0][1], Fe_[1][0], Fe_[1][1];

    svd.compute(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
    U = svd.matrixU();
    V = svd.matrixV();
    w = svd.singularValues();
    Eps << w(0), 0.0f, 0.0f, w(1);
	//end of SVD

	Eigen::Matrix2f dFe = 2.0f * MU_dry_sand * Eps.inverse() * Eps.log() + LAMBDA_dry_sand * (Eps.log().sum()) * Eps.inverse();

    glm::mat2 A;
    glm::mat2 glmV(V(0, 0), V(1, 0), V(0, 1), V(1, 1));
    glm::mat2 glmU(U(0, 0), U(1, 0), U(0, 1), U(1, 1));
    glm::mat2 glmdFe(dFe(0, 0), dFe(1, 0), dFe(0, 1), dFe(1, 1));

    //in glm
    A = v0_ * ((glmU * glmdFe) * glm::transpose(glmV)) * glm::transpose(Fe_);
    return A;
}

void Sand::UpdateDeformation(glm::mat2& T, float dt)
{
	FeTr_ = (glm::mat2(1.0f) + dt * T) * Fe_;
	FpTr_ = Fp_;

	Eigen::Matrix2f U, V, Eps, cov;
	Eigen::Vector2f w;
	Eigen::JacobiSVD<Eigen::Matrix2f> svd;
	cov << FeTr_[0][0], FeTr_[0][1], FeTr_[1][0], FeTr_[1][1];

	svd.compute(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
	U = svd.matrixU();
	V = svd.matrixV();
	w = svd.singularValues();
	Eps << w(0), 0.0f, 0.0f, w(1);

	glm::vec2 T2; 
	double dq;

	Eigen::Matrix2f e = Eps.log();
	glm::vec2 glme(e(0, 0), e(1, 1));
	float e_sum = glme[0] + glme[1];

	glm::vec2 e_c = glme - e_sum / 2.0f * glm::vec2(1.0f, 1.0f);

	if ( glm::length(e_c) < 1e-8  || e_sum > 0) {
		T2 = glm::vec2(1.0f, 1.0f);
		dq = glm::length(glme);
	}
	else {
		float dg = glm::length(e_c) + (LAMBDA_dry_sand + MU_dry_sand) / MU_dry_sand * e_sum * alpha_;

		if (dg <= 0) {
			T2 = glm::vec2(Eps(0, 0), Eps(1, 1));
			dq = 0;
		}
		else {
			glm::vec2 Hm = glme - dg * e_c / glm::length(e_c);
			T2 = glm::vec2(exp(Hm[0]), exp(Hm[1]));
			dq = dg;
		}
	}


	glm::mat2 glmV(V(0, 0), V(1, 0), V(0, 1), V(1, 1));
	glm::mat2 glmU(U(0, 0), U(1, 0), U(0, 1), U(1, 1));
	glm::mat2 glmdT2(T2[0], 0.0f, 0.0f, T2[1]);
	glm::mat2 glmEps(Eps(0, 0), 0.0f, 0.0f, Eps(1, 1));

	// Elastic and plastic state
	Fe_ = glmU * glmdT2 * glm::transpose(glmV);
	Fp_ = glmV * glm::inverse(glmdT2) * glmEps * glm::transpose(glmV) * FpTr_;

	// hardening
	q_ += dq;
	float phi = H0 + (H1 * q_ - H3) * exp(-H2 * q_);
	alpha_ = sqrt(2.0 / 3.0) * (2.0 * sin(phi)) / (3.0 - sin(phi));
}