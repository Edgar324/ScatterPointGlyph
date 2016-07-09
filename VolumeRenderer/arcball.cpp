#include "arcball.h"

void ArcBall::Begin(int x, int y) {
	mVecDown = ScreenToVector(static_cast<float>(x), static_cast<float>(y));
	mQuatDown = mQuatNow;
}

void ArcBall::Drag(int x, int y) {
	mVecNow = ScreenToVector(static_cast<float>(x), static_cast<float>(y));

	glm::vec3 p = glm::cross(mVecDown, mVecNow);

	float l = glm::length(p);
	if (glm::length(p) < 1e-5) return;
	if (glm::length(p) > 1e-5) {
		glm::quat q = glm::quat(glm::dot(mVecDown, mVecNow), p);
		mQuatNow = glm::cross(glm::normalize(q), mQuatDown);
	}
	else {
		mQuatNow = glm::cross(glm::quat(), mQuatDown);
	}
}
//-----------------------------------------------------------------------------
void ArcBall::SetWindowSize(int width, int height) {
	mWidth = width;
	mHeight = height;
}
//-----------------------------------------------------------------------------
glm::vec3 ArcBall::ScreenToVector(float screenX, float screenY) {
	glm::vec2 v;
	v.x = ((screenX / ((mWidth - 1) / 2)) - 1);
	v.y = -((screenY / ((mHeight - 1) / 2)) - 1);

	float len = glm::length(v);
	if (len > 1.0f)
		return glm::vec3(v / sqrt(len), 0);

	return glm::vec3(v, sqrt(1.0f - len));
}