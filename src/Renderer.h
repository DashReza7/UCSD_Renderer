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

    bool hit(const Ray &r, float t_min, float t_max, HitRecord &rec) const;

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

    vec3 get_brdf(const Shape *shape, const vec3 &w_o, const vec3 &w_i, const vec3 &normal) const;

    float get_pdf(ImportanceSamplingType sampling_type, const vec3 &normal, const vec3 &w_i, const vec3 &w_o, const Shape *hit_shape) const;

    vec3 get_sample_dirn(ImportanceSamplingType sampling_type, bool sample_specular, const Shape *hit_shape, const vec3 &w_o, const vec3 &normal);

    // pdf_nee for direct lighting in MIS setting
    float get_pdf_nee(const vec3 &w_i, const HitRecord &world_rec, bool is_world_hit) const;

    // iterate over all area lights, and return if the ray hits any (also taking shadows into account)
    AreaLight *hit_area_light(const Ray &r, vec3 &light_hit_pos) const;

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
