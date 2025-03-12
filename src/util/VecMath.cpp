#include <corecrt_math_defines.h>
#include "VecMath.h"

vec3::vec3(const vec4 &v) {
    x = v.x;
    y = v.y;
    z = v.z;
}

float norm2(const vec2 &v) {
    return dot(v, v);
}

float norm2(const vec3 &v) {
    return dot(v, v);
}

float norm(const vec2 &v) {
    return std::sqrt(norm2(v));
}

float norm(const vec3 &v) {
    return std::sqrt(norm2(v));
}

vec3 normalize(const vec3 &v) {
    return v / norm(v);
}

vec2 normalize(const vec2 &v) {
    return v / norm(v);
}

vec3 cross(const vec3 &v1, const vec3 &v2) {
    return vec3{v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x};
}

float dot(const vec2 &v1, const vec2 &v2) {
    return v1.x * v2.x + v1.y * v2.y;
}

float dot(const vec3 &v1, const vec3 &v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

float triangle_area(const vec2 &v1, const vec2 &v2, const vec2 &v3) {
    vec2 u1 = v2 - v1;
    vec2 u2 = v3 - v1;

    return abs(0.5f * (u1.x * u2.y - u1.y * u2.x));
}

float triangle_area(const vec3 &v1, const vec3 &v2, const vec3 &v3) {
    return 0.5f * norm(cross(v2 - v1, v3 - v1));
}

mat4 transpose(const mat4 &m) {
    return mat4{
        vec4{m.c1.x, m.c2.x, m.c3.x, m.c4.x}, vec4{m.c1.y, m.c2.y, m.c3.y, m.c4.y},
        vec4{m.c1.z, m.c2.z, m.c3.z, m.c4.z}, vec4{m.c1.w, m.c2.w, m.c3.w, m.c4.w}
    };
}

float radians(float degree) {
    return degree * static_cast<float>(M_PI) / 180.0f;
}

mat4 get_translation_matrix(const vec3 &direction) {
    return mat4{
        vec4{1.0f, 0.0f, 0.0f, 0.0f}, vec4{0.0f, 1.0f, 0.0f, 0.0f}, vec4{0.0f, 0.0f, 1.0f, 0.0f},
        vec4{direction, 1.0f}
    };
}

// TODO: change this to a more numerically stable method
mat4 get_rotation_matrix(const vec3 &axis, float angle) {
    vec3 v = normalize(axis);

    vec3 s{};
    if (abs(v.x) <= abs(v.y) && abs(v.x) <= abs(v.z))
        s = vec3{0.0f, -v.z, v.y};
    else if (abs(v.y) <= abs(v.x) && abs(v.y) <= abs(v.z))
        s = vec3{-v.z, 0.0f, v.x};
    else if (abs(v.z) <= abs(v.x) && abs(v.z) <= abs(v.y))
        s = vec3{-v.y, v.x, 0.0f};
    vec3 t = cross(v, s);
    mat4 M = mat4{
        vec4{v, 0.0f},
        vec4{s, 0.0f},
        vec4{t, 0.0f},
        vec4{0.0f, 0.0f, 0.0f, 1.0f}
    };
    mat4 R = mat4{
        vec4{1.0f, 0.0f, 0.0f, 0.0f},
        vec4{0.0f, cos(angle), sin(angle), 0.0f},
        vec4{0.0f, -sin(angle), cos(angle), 0.0f},
        vec4{0.0f, 0.0f, 0.0f, 1.0f}
    };
    return M * R * transpose(M);

    //    float cos_angle = cos(angle);
    //    float one_minus_cos_angle = 1.0f - cos_angle;
    //    float sin_angle = sin(angle);
    //    return mat4{vec4{cos_angle + one_minus_cos_angle * v.x * v.x,
    //                     one_minus_cos_angle * v.x * v.y + sin_angle * v.z,
    //                     one_minus_cos_angle * v.x * v.z - sin_angle * v.y, 0.0f},
    //                vec4{one_minus_cos_angle * v.x * v.y - sin_angle * v.z,
    //                     cos_angle + one_minus_cos_angle * v.y * v.y,
    //                     one_minus_cos_angle * v.y * v.z + sin_angle * v.x, 0.0f},
    //                vec4{one_minus_cos_angle * v.x * v.z + sin_angle * v.y,
    //                     one_minus_cos_angle * v.y * v.z - sin_angle * v.x,
    //                     cos_angle + one_minus_cos_angle * v.z * v.z, 0.0f},
    //                vec4{0.0f, 0.0f, 0.0f, 1.0f}};
}

vec3 reflect(const vec3 &dir, const vec3 &normal) {
    return normal * 2 * dot(dir, normal) - dir;
}

bool point_in_triangle(vec3 p, vec3 v1, vec3 v2, vec3 v3) {
    vec3 n = normalize(cross(v2 - v1, v3 - v1));
    return dot(n, cross(v2 - v1, p - v1)) >= 0.0f &&
           dot(n, cross(v3 - v2, p - v2)) >= 0.0f &&
           dot(n, cross(v1 - v3, p - v3)) >= 0.0f;
}

vec3 get_random_unit_vector(std::mt19937 &generator, std::uniform_real_distribution<float> &uniform_dis) {
    // in the range [0, pi]
    //    float phi = std::numbers::pi * uniform_dis(generator);
    float phi = std::acos(uniform_dis(generator));
    // in the range [0, 2*pi]
    float theta = 2.0f * std::numbers::pi * uniform_dis(generator);

    return vec3{std::sin(phi) * std::cos(theta), std::sin(phi) * std::sin(theta), std::cos(phi)};
}

vec3 get_cosine_unit_vector(std::mt19937 &generator, std::uniform_real_distribution<float> &uniform_dis) {
    float phi = 2 * std::numbers::pi * uniform_dis(generator);
    float theta = std::acos(std::sqrt(uniform_dis(generator)));
    return normalize(vec3{std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta), std::cos(theta)});
}

vec3 get_brdf_unit_vector(std::mt19937 &generator, std::uniform_real_distribution<float> &uniform_dis, bool specular_brdf, float shininess) {
    float phi = 2.0f * std::numbers::pi * uniform_dis(generator);
    float theta = specular_brdf ? std::acos(std::pow(uniform_dis(generator), 1.0f / (shininess + 1.0f))) : std::acos(std::sqrt(uniform_dis(generator)));
    return normalize(vec3{std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta), std::cos(theta)});
}

vec3 get_vector_around_normal(std::mt19937 &generator, std::uniform_real_distribution<float> &uniform_dis, const vec3 &normal, vector_sampling_type type, bool specular_brdf, float shininess) {
    vec3 rnd_vec = type == vector_sampling_type::UNIFORM_HEMISPHERE
                       ? get_random_unit_vector(generator, uniform_dis)
                       : type == vector_sampling_type::COSINE ? get_cosine_unit_vector(generator, uniform_dis) : get_brdf_unit_vector(generator, uniform_dis, specular_brdf, shininess);
    if (type == vector_sampling_type::BRDF && rnd_vec.z < 0.0f)
        return vec3{};
    rnd_vec.z = std::abs(rnd_vec.z);

    return align_vector(rnd_vec, normal);
}

vec3 align_vector(const vec3 &vector, const vec3 &center)
{
    vec3 r = normalize(center);
    vec3 s{};
    if (std::abs(r.x) <= std::abs(r.y) && std::abs(r.x) <= std::abs(r.z))
        s = vec3{0.0f, r.z * -1, r.y};
    else if (std::abs(r.y) <= std::abs(r.x) && std::abs(r.y) <= std::abs(r.z))
        s = vec3{r.z * -1, 0.0f, r.x};
    else if (std::abs(r.z) <= std::abs(r.x) && std::abs(r.z) <= std::abs(r.y))
        s = vec3{r.y * -1, r.x, 0.0f};
    s = normalize(s);
    vec3 t = normalize(cross(r, s));

    return normalize(s * vector.x + t * vector.y + r * vector.z);
}

vec3 get_brdf_unit_ggx(std::mt19937 &generator, std::uniform_real_distribution<float> &uniform_dis, const vec3 &normal, bool specular_brdf, float roughness, const vec3 &wo)
{
    if (specular_brdf)
    {
        float phi = 2.0f * std::numbers::pi_v<float> * uniform_dis(generator);
        float rnd = uniform_dis(generator);
        float theta = std::atan(roughness * std::sqrt(rnd) / std::sqrt(1.0f - rnd));
        vec3 h{std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta), std::cos(theta)};
        vec3 h_rotated = align_vector(h, normal);
        vec3 wi = reflect(wo, h_rotated);
        if (dot(wi, normal) < 0.0f)
            return vec3{};
        return wi;
    }
    else
    {
        float phi = 2 * std::numbers::pi_v<float> * uniform_dis(generator);
        float theta = std::acos(std::sqrt(uniform_dis(generator)));
        vec3 sample{std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta), std::cos(theta)};
        return align_vector(sample, normal);
    }
}

