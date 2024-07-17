#ifndef BORDER_H
#define BORDER_H

#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

struct Plane {
	glm::vec3 normal_;
	float offset_;
};

/*TODO: correct plane location according to the grids and particles*/
class Border {
public:
	Border() {
		Plane left_plane;
		left_plane.normal_ = glm::vec3(1.0f, 0.0f, 0.0f);
		left_plane.offset_ = 2.0f;

		Plane right_plane;
		right_plane.normal_ = glm::vec3(-1.0f, 0.0f, 0.0f);
		right_plane.offset_ = -98.0f;

		Plane bottom_plane;
		bottom_plane.normal_ = glm::vec3(0.0f, 1.0f, 0.0f);
		bottom_plane.offset_ = 2.0f;

		Plane up_plane;
		up_plane.normal_ = glm::vec3(0.0f, -1.0f, 0.0f);
		up_plane.offset_ = -8.0f;

		Plane front_plane;
		front_plane.normal_ = glm::vec3(0.0f, 0.0f, 1.0f);
		front_plane.offset_ = 2.0f;

		Plane rear_plane;
		rear_plane.normal_ = glm::vec3(0.0f, 0.0f, -1.0f);
		rear_plane.offset_ = -98.0f;

		borders_.push_back(left_plane);
		borders_.push_back(right_plane);
		borders_.push_back(bottom_plane);
		borders_.push_back(up_plane);
		borders_.push_back(front_plane);
		borders_.push_back(rear_plane);
	}

	/*TODO: is & need here?*/
	const std::vector<Plane>& GetPlanes() {
		return borders_;
	}

private:
	std::vector<Plane> borders_;

};

#endif //  Border_H

