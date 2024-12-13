#pragma once

#include <random>
#include "Scene.h"
#include "Film.h"
#include "Light.h"

class Renderer
{
private:
    Scene *scene;
    Film *film;
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<float> uniform_dis{0.0f, 1.0f};
    bool use_aliasing;
    uint32_t SPP;
    
    // analytic direct light Integrator
    vec3 get_pixel_color_analyticDirect(const Ray &r) const
    {
        HitRecord rec;
        bool is_hit = scene->bvh->hit(r, T_MIN, T_MAX, rec);
        
        // if hit area light before anything else, return the (clamped) color of the area light
        bool hit_area_light = false;
        float t_area_light = std::numeric_limits<float>::max();
        vec3 arealight_color;
        for (AreaLight al : scene->areaLights)
        {
            HitRecord foo_rec;
            vec3 foo_color;
            if (al.hit(r, foo_rec, foo_color) == true && foo_rec.t < t_area_light)
            {
                hit_area_light = true;
                t_area_light = foo_rec.t;
                arealight_color = foo_color;
            }
        }
        if ((hit_area_light && !is_hit) || (hit_area_light && t_area_light < rec.t))
        {
            arealight_color.clamp(0.0f, 1.0f);
            return arealight_color;
        }
        
        if (!is_hit)
        {
            return vec3{0.0f, 0.0f, 0.0f};
        }
        else
        {
            vec3 color = vec3{0.0f, 0.0f, 0.0f};
            
            auto hit_shape = (Shape *) rec.hit_obj;

//            color += hit_shape->mat.ambient;
//            color += hit_shape->mat.emission;
            
            for (AreaLight al: scene->areaLights)
                color += al.calculate_analytic_radiance(rec, hit_shape->mat.diffuse);
            
            color.clamp(0.0f, 1.0f);
            return color;
        }
    }
    
    // ray tracer Integrator
    vec3 get_pixel_color_raytrace(const Ray &r, uint32_t depth) const
    {
        if (depth == 0)
            return vec3{0.0f, 0.0f, 0.0f};
        
        HitRecord rec;
        bool is_hit = scene->bvh->hit(r, T_MIN, T_MAX, rec);
        if (!is_hit)
        {
            return vec3{0.0f, 0.0f, 0.0f};
        }
        else
        {
            vec3 color = vec3{0.0f, 0.0f, 0.0f};
            
            auto hit_shape = (Shape *) rec.hit_obj;
            
            color += hit_shape->mat.ambient;
            color += hit_shape->mat.emission;
            
            vec3 eye_dirn = normalize(scene->main_camera->position - rec.hit_pos);
            
            HitRecord foo_rec;
            if (scene->dirn_light.dirn_to_light != vec3{0.0f, 0.0f, 0.0f})
            {
                is_hit = scene->bvh->hit(Ray{rec.hit_pos + rec.normal * EPS, scene->dirn_light.dirn_to_light}, T_MIN,
                                         T_MAX, foo_rec);
                if (!is_hit)
                {
                    vec3 tmp_color = vec3{0.0f, 0.0f, 0.0f};
                    
                    tmp_color += hit_shape->mat.diffuse * fmax(0.0f, dot(rec.normal, scene->dirn_light.dirn_to_light));
                    vec3 half_dirn = normalize(scene->dirn_light.dirn_to_light - r.direction);
                    tmp_color += hit_shape->mat.specular *
                                 pow(fmax(0.0f, dot(half_dirn, rec.normal)), hit_shape->mat.shininess);
                    
                    color += tmp_color * scene->dirn_light.color;
                }
            }
            for (auto &pt_light: scene->pointlights)
            {
                float light_distance = norm(pt_light.position - rec.hit_pos);
                vec3 light_dirn = normalize(pt_light.position - rec.hit_pos);
                
                foo_rec.t = std::numeric_limits<float>::max();
                is_hit = scene->bvh->hit(
                        Ray{rec.hit_pos + rec.normal * EPS, normalize(pt_light.position - rec.hit_pos)}, T_MIN, T_MAX,
                        foo_rec);
                if (is_hit && foo_rec.t < light_distance)
                    continue;
                
                vec3 diffuse_color = hit_shape->mat.diffuse * fmax(0.0f, dot(light_dirn, rec.normal));
                
                vec3 half_dirn = normalize(light_dirn - r.direction);
                vec3 specular_color =
                        hit_shape->mat.specular * pow(fmax(0.0f, dot(half_dirn, rec.normal)), hit_shape->mat.shininess);
                
                color += (diffuse_color + specular_color) * pt_light.color * pt_light.calc_attenuation(rec.hit_pos);
            }
            vec3 reflected_dirn = reflect(r.direction * -1, rec.normal);
            Ray reflected_ray = Ray{rec.hit_pos + rec.normal * EPS, reflected_dirn};
            if (hit_shape->mat.specular != vec3{0.0f, 0.0f, 0.0f})
                color += hit_shape->mat.specular * get_pixel_color_raytrace(reflected_ray, depth - 1);
            
            color.clamp(0.0f, 1.0f);
            return color;
        }
    }
    
    // Monte-Carlo direct area light
    vec3 get_pixel_color_direct(const Ray &r)
    {
        HitRecord rec;
        bool is_hit = scene->bvh->hit(r, T_MIN, T_MAX, rec);
        
        // if hit area light before anything else, return the (clamped) color of the area light
        bool hit_area_light = false;
        float t_area_light = std::numeric_limits<float>::max();
        vec3 arealight_color;
        for (AreaLight al : scene->areaLights)
        {
            HitRecord foo_rec;
            vec3 foo_color;
            if (al.hit(r, foo_rec, foo_color) == true && foo_rec.t < t_area_light)
            {
                hit_area_light = true;
                t_area_light = foo_rec.t;
                arealight_color = foo_color;
            }
        }
        if ((hit_area_light && !is_hit) || (hit_area_light && t_area_light < rec.t))
        {
            arealight_color.clamp(0.0f, 1.0f);
            return arealight_color;
        }
        
        if (!is_hit)
        {
            return vec3{0.0f, 0.0f, 0.0f};
        }
        else
        {
            vec3 color{};
            vec3 reflected_ray_dir = reflect(r.direction * -1, rec.normal);
            auto hit_shape = (Shape *) rec.hit_obj;
            for (AreaLight al: scene->areaLights)
            {
                std::vector<vec3> rnd_light_points;
                if (!scene->light_stratify)
                    for (int i = 0; i < scene->light_samples; ++i)
                    {
                        float rnd1 = uniform_dis(gen);
                        float rnd2 = uniform_dis(gen);
                        rnd_light_points.push_back(al.a + (al.b - al.a) * rnd1 + (al.c - al.a) * rnd2);
                    }
                else
                {
                    int root_N = std::round(std::sqrt(scene->light_samples));
                    for (int i = 0; i < root_N; ++i)
                        for (int j = 0; j < root_N; ++j)
                        {
                           float rnd1 = uniform_dis(gen);
                           float rnd2 = uniform_dis(gen);
                            rnd_light_points.push_back(al.a + (al.b - al.a) * i / root_N + (al.c - al.a) * j / root_N + (al.b - al.a) * rnd1 / root_N + (al.c - al.a) * rnd2 / root_N);
                        }
                }
                
                vec3 sum_value;
                for (int i = 0; i < scene->light_samples; ++i)
                {
                    vec3 rnd_light_pos = rnd_light_points[i];
                    
                    vec3 surface_to_light_vec = rnd_light_pos - rec.hit_pos;
                    float light_dist = norm(surface_to_light_vec);
                    surface_to_light_vec = normalize(surface_to_light_vec);
                    
                    float cos_theta_i = dot(rec.normal, surface_to_light_vec);
                    cos_theta_i = std::clamp(cos_theta_i, 0.0f, 1.0f);
                    float cos_theta_o = dot(surface_to_light_vec, al.normal);
                    cos_theta_o = std::clamp(cos_theta_o, 0.0f, 1.0f);
                    
                    float geom_term =
                            cos_theta_i * cos_theta_o / dot(rec.hit_pos - rnd_light_pos, rec.hit_pos - rnd_light_pos);

                    HitRecord foo_rec;
                    is_hit = scene->bvh->hit(Ray{rec.hit_pos + rec.normal * EPS, surface_to_light_vec}, T_MIN, T_MAX,
                                             foo_rec);
                    if (is_hit && foo_rec.t < light_dist)
                        continue;
                    sum_value += (hit_shape->mat.diffuse / std::numbers::pi + hit_shape->mat.specular / (2 * std::numbers::pi) * (hit_shape->mat.shininess + 2) * std::pow(dot(reflected_ray_dir, surface_to_light_vec), hit_shape->mat.shininess)) * geom_term;
                }
                color += sum_value * al.area / scene->light_samples * al.color;
            }
            color.clamp(0.0f, 1.0f);
            return color;
        }
        
        return vec3{};
    }
    
    // handle which Integrator to use
    vec3 get_pixel_color(const Ray &r, uint32_t depth)
    {
        if (scene->integrator == "raytracer")
            return get_pixel_color_raytrace(r, depth);
        else if (scene->integrator == "analyticdirect")
            return get_pixel_color_analyticDirect(r);
        else if (scene->integrator == "direct")
            return get_pixel_color_direct(r);
        return vec3{};
    }

public:
    Renderer(Scene *s, Film *f, bool use_aliasing = false, uint32_t SPP = 10)
    {
        scene = s;
        film = f;
        gen = std::mt19937(rd());
        this->use_aliasing = use_aliasing;
        this->SPP = SPP;
    }
    
    void render()
    {
        for (uint32_t i = 0; i < scene->height; ++i)
        {
            for (uint32_t j = 0; j < scene->width; ++j)
            {
                if (use_aliasing)
                {
                    // TODO: implement blue noise instead of uniform noise
                    vec3 pixel_color = vec3{0.0f, 0.0f, 0.0f};
                    for (int k = 0; k < SPP; ++k)
                    {
                        float u = (static_cast<float>(j) + uniform_dis(gen)) / static_cast<float>(scene->width);
                        float v = (static_cast<float>(i) + uniform_dis(gen)) / static_cast<float>(scene->height);
                        pixel_color += get_pixel_color(scene->main_camera->get_ray(u, v), scene->maxdepth);
                    }
                    pixel_color = pixel_color / static_cast<float>(SPP);
                    film->commit(j, i, pixel_color);
                }
                else
                {
                    auto u = (static_cast<float>(j) + 0.5f) / static_cast<float>(scene->width);
                    auto v = (static_cast<float>(i) + 0.5f) / static_cast<float>(scene->height);
                    vec3 pixel_color = get_pixel_color(scene->main_camera->get_ray(u, v), scene->maxdepth);
                    film->commit(j, i, pixel_color);
                }
            }
            std::cout << std::format("Progress: {:.2f}%\r", (float) i / (float) scene->height * 100.0f);
        }
        std::cout << "Progress: 100.0%" << std::endl;
    }
};
