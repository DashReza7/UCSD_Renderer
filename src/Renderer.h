#pragma once

#include <random>
#include "Scene.h"
#include "Light.h"

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
    
    // Just a simple path-tracer
    vec3 get_pixel_color_pathtrace(const Ray &r, uint32_t depth);
    
    // Path tracer with Next Event Estimation
    vec3 get_pixel_color_pathtrace_NEE(const Ray &r, uint32_t depth);
    
    // handle which Integrator to use
    vec3 get_pixel_color(const Ray &r, uint32_t depth);

public:
    Renderer(Scene *s, Film *f)
    {
        scene = s;
        film = f;
        gen = std::mt19937(rd());
    }
    
    void render();
};
