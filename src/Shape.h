#pragma once

#include <algorithm>
#include "util/VecMath.h"
#include "SpatialDatastructures.h"

class Material
{
public:
    vec3 diffuse;
    vec3 specular;
    float shininess;
    vec3 ambient;
    vec3 emission;
    
    Material() = default;
    Material(const vec3& diffuse, const vec3& specular, float shininess, const vec3& ambient, const vec3& emission) : diffuse(diffuse), specular(specular), shininess(shininess), ambient(ambient), emission(emission) {}
};

class Shape : public Hittable
{
public:
    Material mat;
};

class Sphere : public Shape
{
private:
    mat4 transform;
    mat4 inverse_transform;
    
public:
    vec3 center;
    float radius;
    
    Sphere() : center(vec3(0.0f, 0.0f, 0.0f)), radius(1.0f)
    {
        set_bbox();
    }
    Sphere(const vec3& cen, float r, const Material& mat, const mat4& transform, const mat4& inverse_transform) : center(cen), radius(r), transform(transform), inverse_transform(inverse_transform)
    {
        this->mat = mat;
        set_bbox();
    }
    
    bool hit(const Ray& r, float t_min, float t_max, HitRecord& rec) override;
    
    void set_bbox();
};

class Triangle : public Shape
{
public:
    vec3 v1;
    vec3 v2;
    vec3 v3;
    vec3 n;
    float area;
    
    Triangle() : v1(vec3{}), v2(vec3{}), v3(vec3{}), n(vec3{}), area(0.0f)
    {
        set_bbox();
    }
    Triangle(const vec3& v1, const vec3& v2, const vec3& v3, const Material& mat) : v1(v1), v2(v2), v3(v3)
    {
        this->mat = mat;
        n = normalize(cross(v2 - v1, v3 - v1));
        area = triangle_area(v1, v2, v3);
        set_bbox();
    }
    
    bool hit(const Ray& r, float t_min, float t_max, HitRecord& rec) override;
    
    void set_bbox();
};

