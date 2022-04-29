#pragma once
#include <cmath>

struct Vector3f
{
	float x;
	float y;
	float z;

public:
	Vector3f(float _x, float _y, float _z) :x(_x), y(_y), z(_z) {}
	Vector3f() :Vector3f(0.0f, 0.0f, 0.0f) {}



	void set(float _x, float _y, float _z)
	{
		x = _x, y = _y, z = _z;
	}

	Vector3f normalize() const
	{
		float l = std::sqrt(x * x + y * y + z * z);
		return Vector3f(x / l, y / l, z / l);
	}

	Vector3f operator+(const Vector3f& v) const
	{
		return Vector3f(x + v.x, y + v.y, z + v.z);
	}
	Vector3f operator-(const Vector3f& v) const
	{
		return Vector3f(x - v.x, y - v.y, z - v.z);
	}
	Vector3f crossed(const Vector3f& v) const
	{
		return Vector3f(
			y * v.z - z * v.y,
			z * v.x - x * v.z,
			x * v.y - y * v.x
		);
	}
	Vector3f operator*(const float& t)const
	{
		return Vector3f(x * t, y * t, z * t);
	}
	Vector3f operator/(const float& t)const
	{
		return Vector3f(x / t, y / t, z / t);
	}

};

