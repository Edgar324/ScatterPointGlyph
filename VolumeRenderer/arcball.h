#ifndef ARCBALL_H
#define ARCBALL_H

#include "glm/glm.hpp"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/quaternion.hpp>

#include <glm/gtx/transform.hpp>
#include <glm/gtx/transform2.hpp>
#include <glm/gtx/euler_angles.hpp>

class ArcBall {
public:
	ArcBall() : mWidth(0), mHeight(0) {}

	void Begin(int x, int y);
	void Drag(int x, int y);
	void SetWindowSize(int width, int height);
	const glm::mat3 Rotation() const { return glm::mat3_cast(mQuatNow); }

private:
	glm::vec3 ScreenToVector(float screenX, float screenY);

	int mWidth, mHeight;
	glm::quat mQuatNow, mQuatDown;
	glm::vec3 mVecNow, mVecDown;
};

#endif