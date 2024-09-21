#include <corecrt_math_defines.h>
#include "VecMath.h"

vec3::vec3(const vec4 &v)
{
    x = v.x;
    y = v.y;
    z = v.z;
}

float norm2(const vec2 &v)
{
    return dot(v, v);
}
float norm2(const vec3 &v)
{
    return dot(v, v);
}

float norm(const vec2 &v)
{
    return std::sqrt(norm2(v));
}
float norm(const vec3 &v)
{
    return std::sqrt(norm2(v));
}

vec3 normalize(const vec3 &v)
{
    return v / norm(v);
}
vec2 normalize(const vec2 &v)
{
    return v / norm(v);
}

vec3 cross(const vec3 &v1, const vec3 &v2)
{
    return vec3{v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x};
}

float dot(const vec2 &v1, const vec2 &v2)
{
    return v1.x * v2.x + v1.y * v2.y;
}
float dot(const vec3 &v1, const vec3 &v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

float triangle_area(const vec2 &v1, const vec2 &v2, const vec2 &v3)
{
    vec2 u1 = v2 - v1;
    vec2 u2 = v3 - v1;
    
    return abs(0.5f * (u1.x * u2.y - u1.y * u2.x));
}
float triangle_area(const vec3 &v1, const vec3 &v2, const vec3 &v3)
{
    return 0.5f * norm(cross(v2 - v1, v3 - v1));
}

mat4 transpose(const mat4 &m)
{
    return mat4{vec4{m.c1.x, m.c2.x, m.c3.x, m.c4.x}, vec4{m.c1.y, m.c2.y, m.c3.y, m.c4.y},
                vec4{m.c1.z, m.c2.z, m.c3.z, m.c4.z}, vec4{m.c1.w, m.c2.w, m.c3.w, m.c4.w}};
}

float radians(float degree)
{
    return degree * static_cast<float>(M_PI) / 180.0f;
}

mat4 get_translation_matrix(const vec3 &direction)
{
    return mat4{vec4{1.0f, 0.0f, 0.0f, 0.0f}, vec4{0.0f, 1.0f, 0.0f, 0.0f}, vec4{0.0f, 0.0f, 1.0f, 0.0f},
                vec4{direction, 1.0f}};
}

// TODO: change this to a more numerically stable method
mat4 get_rotation_matrix(const vec3 &axis, float angle)
{
    vec3 axis_normalized = normalize(axis);
    float cos_angle = cos(angle);
    float one_minus_cos_angle = 1.0f - cos_angle;
    float sin_angle = sin(angle);
    return mat4{vec4{cos_angle + one_minus_cos_angle * axis_normalized.x * axis_normalized.x,
                     one_minus_cos_angle * axis_normalized.x * axis_normalized.y + sin_angle * axis_normalized.z,
                     one_minus_cos_angle * axis_normalized.x * axis_normalized.z - sin_angle * axis_normalized.y, 0.0f},
                vec4{one_minus_cos_angle * axis_normalized.x * axis_normalized.y - sin_angle * axis_normalized.z,
                     cos_angle + one_minus_cos_angle * axis_normalized.y * axis_normalized.y,
                     one_minus_cos_angle * axis_normalized.y * axis_normalized.z + sin_angle * axis_normalized.x, 0.0f},
                vec4{one_minus_cos_angle * axis_normalized.x * axis_normalized.z + sin_angle * axis_normalized.y,
                     one_minus_cos_angle * axis_normalized.y * axis_normalized.z - sin_angle * axis_normalized.x,
                     cos_angle + one_minus_cos_angle * axis_normalized.z * axis_normalized.z, 0.0f},
                vec4{0.0f, 0.0f, 0.0f, 1.0f}};
}

vec3 reflect(const vec3 &dir, const vec3 &normal)
{
    return normal * 2 * dot(dir, normal) - dir;
}

bool point_in_triangle(vec3 p, vec3 v1, vec3 v2, vec3 v3)
{
    vec3 n = normalize(cross(v2 - v1, v3 - v1));
    return dot(n, cross(v2 - v1, p - v1)) >= 0.0f &&
           dot(n, cross(v3 - v2, p - v2)) >= 0.0f &&
           dot(n, cross(v1 - v3, p - v3)) >= 0.0f;
}

