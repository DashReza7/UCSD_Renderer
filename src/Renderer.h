#pragma once

#include <random>
#include "Scene.h"

class Renderer
{
private:
    Scene *scene;
    Film *film;
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<float> uniform_dis{0.0f, 1.0f};
    
    bool hit_brute_force(const Ray &r, float t_min, float t_max, HitRecord &rec) const;
    
    // ray tracer Integrator
    vec3 get_pixel_color_raytrace(const Ray &r, uint32_t depth) const;
    
    // analytic direct light Integrator
    vec3 get_pixel_color_analyticDirect(const Ray &r) const;
    
    // Monte-Carlo direct area light
    vec3 get_pixel_color_direct(const Ray &r);
    
    // Path-tracer
    vec3 get_pixel_color_pathtrace(const Ray &r, uint32_t depth, vec3 incoming_throughput = vec3{1.0f, 1.0f, 1.0f});

    // handle which Integrator to use
    vec3 get_pixel_color(const Ray &r, uint32_t depth);

    void render_parallel();

    void render_sequential();

    vec3 get_brdf(const Shape *shape, const vec3 &w_o, const vec3 &w_i, const vec3 &normal);

    // vec3 get_bounced_ray_dirn() const;

    float get_pdf(ImportanceSamplingType importance_sampling_type, const vec3 &normal, const vec3 &w_i, const vec3 &w_o, const Shape *hit_shape) const;

    vec3 sample_dirn(bool sample_specular, const Shape *hit_shape, const vec3 &w_o, const vec3 &normal);

public:
    Renderer(Scene *s, Film *f)
    {
        scene = s;
        film = f;
        gen = std::mt19937(rd());
    }
    
    void render()
    {
        if (scene->parallel_run)
            render_parallel();
        else
            render_sequential();
    }
};
