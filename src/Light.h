#pragma once

#include "util/VecMath.h"
#include <numbers>

class Light
{
public:
    vec3 color;
    
    // TODO: remove this.
    virtual float calc_attenuation(vec3 pos) const = 0;
};

class DirnLight : public Light
{
public:
    vec3 dirn_to_light;
    
    DirnLight() : dirn_to_light(vec3{0.0f, 0.0f, 0.0f})
    {}
    
    DirnLight(vec3 color, vec3 direction_to_light) : dirn_to_light(direction_to_light)
    {
        this->color = color;
    }
    
    float calc_attenuation(vec3 pos) const override
    {
        return 1.0f;
    }
};

class PointLight : public Light
{
public:
    vec3 position;
    float atten[3]{};
    
    PointLight() : position(vec3{})
    {
        atten[0] = 1.0f;
        atten[1] = 0.0f;
        atten[2] = 0.0f;
    }
    
    PointLight(vec3 color, vec3 position, float constant, float linear, float quadratic) : position(position)
    {
        this->color = color;
        atten[0] = constant;
        atten[1] = linear;
        atten[2] = quadratic;
    }
    
    float calc_attenuation(vec3 obj_posn) const override
    {
        float distance = norm(obj_posn - position);
        return 1.0f / (atten[0] + atten[1] * distance + atten[2] * distance * distance);
    }
};

class AreaLight : public Light
{
public:
    // vertices in order: abdc
    vec3 a;
    vec3 b;
    vec3 c;
    vec3 d;
    float area;
    vec3 normal;
    
    AreaLight() {}
    AreaLight(vec3 color, vec3 a, vec3 b, vec3 c, vec3 d) : a(a), b(b), c(c), d(d)
    {
        this->color = color;
        area = triangle_area(a, b, d) + triangle_area(a, c, d);
        normal = cross(b - a, c - a);
        normal = normalize(normal);
    }
    
    float calc_attenuation(vec3 obj_posn) const override
    {
        return 1.0f;
    }
    
    vec3 calculate_analytic_radiance(const HitRecord &rec, const vec3 &hitpoint_diffuse) const
    {
        float theta1 = acos(dot(normalize(a - rec.hit_pos), normalize(b - rec.hit_pos)));
        float theta2 = acos(dot(normalize(b - rec.hit_pos), normalize(d - rec.hit_pos)));
        float theta3 = acos(dot(normalize(d - rec.hit_pos), normalize(c - rec.hit_pos)));
        float theta4 = acos(dot(normalize(c - rec.hit_pos), normalize(a - rec.hit_pos)));
        vec3 gamma1 = cross(a - rec.hit_pos, b - rec.hit_pos);
        gamma1 = normalize(gamma1);
        vec3 gamma2 = cross(b - rec.hit_pos, d - rec.hit_pos);
        gamma2 = normalize(gamma2);
        vec3 gamma3 = cross(d - rec.hit_pos, c - rec.hit_pos);
        gamma3 = normalize(gamma3);
        vec3 gamma4 = cross(c - rec.hit_pos, a - rec.hit_pos);
        gamma4 = normalize(gamma4);
        
        vec3 phi = (gamma1 * theta1 + gamma2 * theta2 + gamma3 * theta3 + gamma4 * theta4) / 2.0f;
        return hitpoint_diffuse / std::numbers::pi * color * dot(phi, rec.normal);
    }
    
    bool hit(const Ray &r, HitRecord &rec, vec3 &color)
    {
        color = this->color;
        // TODO: maybe negate the normal?
        float t = dot(a - r.origin, normal) / dot(r.direction, normal);
        if (t < 0)
            return false;
        rec.t = t;
        rec.hit_pos = r(t);
        rec.normal = normal;
        if (point_in_triangle(rec.hit_pos, a, b, d))
            return true;
        if (point_in_triangle(rec.hit_pos, a, c, d))
            return true;
        return false;
    }
};
