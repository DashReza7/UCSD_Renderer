#include <tbb/tbb.h>
#include <random>
#include "Scene.h"
#include "Light.h"
#include "Renderer.h"

void Renderer::render_parallel()
{
    std::atomic<int> render_progress = 0;
    tbb::parallel_for(0, static_cast<int>(scene->height), [&](int i)
    {
        for (uint32_t j = 0; j < scene->width; ++j)
        {
            vec3 pixel_color = vec3{0.0f, 0.0f, 0.0f};
            for (int k = 0; k < scene->spp; ++k)
            {
                float offset1 = k == 0 ? 0.5f : uniform_dis(gen);
                float offset2 = k == 0 ? 0.5f : uniform_dis(gen);
                float u = (static_cast<float>(j) + offset1) / static_cast<float>(scene->width);
                float v = (static_cast<float>(i) + offset2) / static_cast<float>(scene->height);
                pixel_color += get_pixel_color(scene->main_camera->get_ray(u, v), scene->maxdepth);
            }
            pixel_color = pixel_color / static_cast<float>(scene->spp);
            pixel_color.clamp(0.0f, 1.0f);
            film->commit(j, i, pixel_color);
        }
        render_progress.fetch_add(1, std::memory_order_relaxed);
        std::cout << std::format("Progress: {:.2f}%\r", (float) render_progress / (float) scene->height * 100.0f);
    });
    std::cout << "Progress: 100.0%" << std::endl;
}

void Renderer::render_sequential()
{
    for (uint32_t i = 0; i < scene->height; ++i)
    {
        for (uint32_t j = 0; j < scene->width; ++j)
        {
            vec3 pixel_color = vec3{0.0f, 0.0f, 0.0f};
            for (int k = 0; k < scene->spp; ++k)
            {
                float offset1 = k == 0 ? 0.5f : uniform_dis(gen);
                float offset2 = k == 0 ? 0.5f : uniform_dis(gen);
                float u = (static_cast<float>(j) + offset1) / static_cast<float>(scene->width);
                float v = (static_cast<float>(i) + offset2) / static_cast<float>(scene->height);
                pixel_color += get_pixel_color(scene->main_camera->get_ray(u, v), scene->maxdepth);
            }
            pixel_color = pixel_color / static_cast<float>(scene->spp);
            pixel_color.clamp(0.0f, 1.0f);
            film->commit(j, i, pixel_color);
        }
        std::cout << std::format("Progress: {:.2f}%\r", (float) i / (float) scene->height * 100.0f);
    }
    std::cout << "Progress: 100.0%" << std::endl;
}

bool Renderer::hit_brute_force(const Ray &r, float t_min, float t_max, HitRecord &rec) const
{
    rec.t = t_max;
    HitRecord tmp_rec;
    bool hit_flag = false;
    for (Shape *s: scene->shapes)
    {
        if (s->hit(r, t_min, t_max, tmp_rec) && tmp_rec.t < rec.t)
        {
            hit_flag = true;
            rec = tmp_rec;
            rec.hit_obj = s;
        }
    }

    if (hit_flag)
        return true;
    return false;
}

// ray tracer Integrator
vec3 Renderer::get_pixel_color_raytrace(const Ray &r, uint32_t depth) const
{
    if (depth == 0)
        return vec3{};

    HitRecord rec;
    bool is_hit = scene->use_bvh ? scene->bvh->hit(r, T_MIN, T_MAX, rec) : hit_brute_force(r, T_MIN, T_MAX, rec);
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
            Ray shadow_ray{rec.hit_pos + rec.normal * EPS, scene->dirn_light.dirn_to_light};
            is_hit = scene->use_bvh
                         ? scene->bvh->hit(shadow_ray, T_MIN, T_MAX, foo_rec)
                         : hit_brute_force(r, T_MIN,
                                           T_MAX,
                                           foo_rec);
            if (!is_hit)
            {
                vec3 tmp_color = vec3{0.0f, 0.0f, 0.0f};

                tmp_color += hit_shape->mat.diffuse * fmax(0.0f, dot(rec.normal, scene->dirn_light.dirn_to_light));
                vec3 half_dirn = normalize(scene->dirn_light.dirn_to_light - r.direction);
                tmp_color +=
                        hit_shape->mat.specular * pow(fmax(0.0f, dot(half_dirn, rec.normal)), hit_shape->mat.shininess);

                color += tmp_color * scene->dirn_light.color;
            }
        }
        for (auto &pt_light: scene->pointlights)
        {
            float light_distance = norm(pt_light.position - rec.hit_pos);
            vec3 light_dirn = normalize(pt_light.position - rec.hit_pos);

            foo_rec.t = std::numeric_limits<float>::max();
            Ray shadow_ray{rec.hit_pos + rec.normal * EPS, normalize(pt_light.position - rec.hit_pos)};
            is_hit = scene->use_bvh
                         ? scene->bvh->hit(shadow_ray, T_MIN, T_MAX, foo_rec)
                         : hit_brute_force(shadow_ray,
                                           T_MIN, T_MAX,
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

// analytic direct light Integrator
vec3 Renderer::get_pixel_color_analyticDirect(const Ray &r) const
{
    HitRecord rec;
    bool is_hit = scene->use_bvh ? scene->bvh->hit(r, T_MIN, T_MAX, rec) : hit_brute_force(r, T_MIN, T_MAX, rec);

    // if hit area light before anything else, return the (clamped) color of the area light
    bool hit_area_light = false;
    float t_area_light = std::numeric_limits<float>::max();
    vec3 arealight_color;
    for (AreaLight al: scene->areaLights)
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

// Monte-Carlo direct area light
vec3 Renderer::get_pixel_color_direct(const Ray &r)
{
    HitRecord rec;
    bool is_hit = scene->use_bvh ? scene->bvh->hit(r, T_MIN, T_MAX, rec) : hit_brute_force(r, T_MIN, T_MAX, rec);

    // if hit area light before anything else, return the (clamped) color of the area light
    bool hit_area_light = false;
    float t_area_light = std::numeric_limits<float>::max();
    vec3 arealight_color;
    for (AreaLight al: scene->areaLights)
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

    if (is_hit)
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
                        rnd_light_points.push_back(al.a + (al.b - al.a) * i / root_N + (al.c - al.a) * j / root_N +
                                                   (al.b - al.a) * rnd1 / root_N + (al.c - al.a) * rnd2 / root_N);
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
                Ray shadow_ray{rec.hit_pos + rec.normal * EPS, surface_to_light_vec};
                is_hit = scene->use_bvh
                             ? scene->bvh->hit(shadow_ray, T_MIN, T_MAX, foo_rec)
                             : hit_brute_force(
                                 shadow_ray, T_MIN, T_MAX, foo_rec);
                if (is_hit && foo_rec.t < light_dist)
                    continue;
                sum_value += (hit_shape->mat.diffuse / std::numbers::pi +
                              hit_shape->mat.specular / (2 * std::numbers::pi) * (hit_shape->mat.shininess + 2) *
                              std::pow(dot(reflected_ray_dir, surface_to_light_vec), hit_shape->mat.shininess)) *
                        geom_term;
            }
            color += sum_value * al.area / scene->light_samples * al.color;
        }
        color.clamp(0.0f, 1.0f);
        return color;
    }

    return vec3{};
}

// Path-tracer
vec3 Renderer::get_pixel_color_pathtrace(const Ray &r, uint32_t depth, vec3 incoming_throughput)
{
    if (depth == 0)
        return vec3{};

    HitRecord rec;
    bool is_hit = scene->use_bvh ? scene->bvh->hit(r, T_MIN, T_MAX, rec) : hit_brute_force(r, T_MIN, T_MAX, rec);

    bool hit_area_light = false;
    float t_area_light = std::numeric_limits<float>::max();
    vec3 arealight_color;
    for (AreaLight al: scene->areaLights)
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
        if (!scene->next_event_estimation || depth == scene->maxdepth)
            return arealight_color;
        return vec3{};
    }
    if (!is_hit)
        return vec3{};

    auto hit_shape = (Shape *) rec.hit_obj;

    vec3 color{};
    vec3 reflected_ray_dir = reflect(r.direction * -1, rec.normal);

    // direct light
    if (scene->next_event_estimation)
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
                        rnd_light_points.push_back(al.a + (al.b - al.a) * i / root_N + (al.c - al.a) * j / root_N +
                                                   (al.b - al.a) * rnd1 / root_N + (al.c - al.a) * rnd2 / root_N);
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
                Ray shadow_ray{rec.hit_pos + rec.normal * EPS, surface_to_light_vec};
                is_hit = scene->use_bvh
                             ? scene->bvh->hit(shadow_ray, T_MIN, T_MAX, foo_rec)
                             : hit_brute_force(
                                 shadow_ray, T_MIN, T_MAX, foo_rec);
                if (is_hit && foo_rec.t < light_dist)
                    continue;
                sum_value += (hit_shape->mat.diffuse / std::numbers::pi +
                              hit_shape->mat.specular / (2 * std::numbers::pi) * (hit_shape->mat.shininess + 2) *
                              std::pow(dot(reflected_ray_dir, surface_to_light_vec), hit_shape->mat.shininess)) *
                        geom_term;
            }
            color += sum_value * al.area / scene->light_samples * al.color;
        }

    // indirect light
    float reflectivity = hit_shape->mat.specular.x + hit_shape->mat.specular.y + hit_shape->mat.specular.z;
    reflectivity = reflectivity / (reflectivity + hit_shape->mat.diffuse.x + hit_shape->mat.diffuse.y + hit_shape->mat.diffuse.z);
    bool sample_specular = uniform_dis(gen) < reflectivity;
    vec3 center_vec = (scene->importance_sampling_type == vector_sampling_type::BRDF && sample_specular)
                          ? reflected_ray_dir
                          : rec.normal;
    Ray reflected_ray{
        rec.hit_pos + rec.normal * EPS,
        get_vector_around_normal(gen, uniform_dis, center_vec, scene->importance_sampling_type,
                                 sample_specular, hit_shape->mat.shininess)
    };

    if (reflected_ray.direction == vec3{})
        return color;

    vec3 brdf = (hit_shape->mat.diffuse / std::numbers::pi + hit_shape->mat.specular * (
                     hit_shape->mat.shininess + 2) / (2 * std::numbers::pi) * std::pow(
                     dot(reflected_ray_dir, reflected_ray.direction), hit_shape->mat.shininess));
    vec3 throughput{};
    if (scene->importance_sampling_type == vector_sampling_type::UNIFORM_HEMISPHERE)
        throughput = brdf * (2 * std::numbers::pi) * dot(rec.normal, reflected_ray.direction);
    else if (scene->importance_sampling_type == vector_sampling_type::COSINE)
        throughput = brdf * std::numbers::pi;
    else if (scene->importance_sampling_type == vector_sampling_type::BRDF)
    {
        float pdf = 0.0f;
        // if (uniform_dis(gen) <= reflectivity)
        //     pdf = (hit_shape->mat.shininess + 1) / (2 * std::numbers::pi) * std::pow(
        //               dot(reflected_ray_dir, reflected_ray.direction), hit_shape->mat.shininess);
        // else
        //     pdf = dot(rec.normal, reflected_ray.direction) / std::numbers::pi;
        pdf = (1.0f - reflectivity) * dot(rec.normal, reflected_ray.direction) / std::numbers::pi + reflectivity * (
                  hit_shape->mat.shininess + 1) / (2 * std::numbers::pi) * std::pow(
                  dot(reflected_ray_dir, reflected_ray.direction), hit_shape->mat.shininess);
        throughput = brdf * dot(rec.normal, reflected_ray.direction) / pdf;
    }
    incoming_throughput *= throughput;
    float termination_probability = scene->russian_roulette
                                        ? 1.0f - std::min(1.0f, std::max(incoming_throughput.x,
                                                                         std::max(incoming_throughput.y,
                                                                             incoming_throughput.z)))
                                        : 0.0f;
    if (!scene->russian_roulette || uniform_dis(gen) > termination_probability)
    {
        vec3 traced_color = get_pixel_color_pathtrace(reflected_ray, depth - 1, incoming_throughput / (1 - termination_probability));
        color += traced_color * throughput / (1 - termination_probability);
    }
    return color;
}

vec3 Renderer::get_pixel_color_ggx(const Ray &r, uint32_t depth, vec3 incoming_throughput)
{
    if (depth == 0)
        return vec3{};

    HitRecord rec;
    bool is_hit = scene->use_bvh ? scene->bvh->hit(r, T_MIN, T_MAX, rec) : hit_brute_force(r, T_MIN, T_MAX, rec);

    bool hit_area_light = false;
    float t_area_light = std::numeric_limits<float>::max();
    vec3 arealight_color;
    for (AreaLight al: scene->areaLights)
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
        if (!scene->next_event_estimation || depth == scene->maxdepth)
            return arealight_color;
        return vec3{};
    }
    if (!is_hit)
        return vec3{};

    auto hit_shape = (Shape *) rec.hit_obj;

    vec3 color{};
    vec3 reflected_ray_dir = reflect(r.direction * -1, rec.normal);

    // direct light
    if (scene->next_event_estimation)
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
                        rnd_light_points.push_back(al.a + (al.b - al.a) * i / root_N + (al.c - al.a) * j / root_N +
                                                   (al.b - al.a) * rnd1 / root_N + (al.c - al.a) * rnd2 / root_N);
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
                Ray shadow_ray{rec.hit_pos + rec.normal * EPS, surface_to_light_vec};
                is_hit = scene->use_bvh
                             ? scene->bvh->hit(shadow_ray, T_MIN, T_MAX, foo_rec)
                             : hit_brute_force(
                                 shadow_ray, T_MIN, T_MAX, foo_rec);
                if (is_hit && foo_rec.t < light_dist)
                    continue;
                vec3 brdf{};
                if (hit_shape->mat.brdf_type == BRDF_TYPE::Phong)
                {
                    brdf = hit_shape->mat.diffuse / std::numbers::pi + hit_shape->mat.specular / (2 * std::numbers::pi) * (hit_shape->mat.shininess + 2) * std::pow(dot(reflected_ray_dir, surface_to_light_vec), hit_shape->mat.shininess);
                }
                else if (hit_shape->mat.brdf_type == BRDF_TYPE::GGX)
                {
                    vec3 half_vector = normalize(surface_to_light_vec - r.direction);
                    float tan_sqrd_h = dot(rec.normal, half_vector) != 0.0f ? 1.0f / (std::pow(dot(rec.normal, half_vector), 2)) - 1.0f : 1000.0f;
                    float alpha_sqrd = std::pow(hit_shape->mat.roughness, 2);
                    float mic_dis_func = 0.0f;
                    if (dot(rec.normal, half_vector) <= EPS || alpha_sqrd + tan_sqrd_h <= EPS)
                        mic_dis_func = 1.0f;
                    else
                        mic_dis_func = alpha_sqrd / (std::numbers::pi_v<float> * std::pow(dot(rec.normal, half_vector), 4) * (std::pow(alpha_sqrd + tan_sqrd_h, 2)));
                    float tan_sqrd_wi = dot(rec.normal, surface_to_light_vec) != 0 ? 1.0f / std::pow(dot(rec.normal, surface_to_light_vec), 2) - 1.0f : 1000.0f;
                    float tan_sqrd_wo = dot(rec.normal, r.direction * -1.0f) != 0 ? 1.0f / std::pow(dot(rec.normal, r.direction * -1.0f), 2) - 1.0f : 1000.0f;
                    float smith_g_func = (2.0f / (1.0f + std::sqrt(1.0f + alpha_sqrd * tan_sqrd_wi))) * (2.0f / (1.0f + std::sqrt(1.0f + alpha_sqrd * tan_sqrd_wo)));
                    vec3 fresnel_term = hit_shape->mat.specular + (vec3{1.0f, 1.0f, 1.0f} - hit_shape->mat.specular) * std::pow(1 - dot(surface_to_light_vec, half_vector), 5);
                    vec3 brdf_specular{};
                    if (std::abs(dot(r.direction * -1.0f, rec.normal)) <= EPS || std::abs(dot(surface_to_light_vec, rec.normal)) <= EPS)
                        brdf_specular = vec3{1.0f, 1.0f, 1.0f};
                    else
                        brdf_specular = (fresnel_term * smith_g_func * mic_dis_func) / (4.0f * dot(surface_to_light_vec, rec.normal) * dot(r.direction * -1.0f, rec.normal));
                    brdf = brdf_specular + hit_shape->mat.diffuse / std::numbers::pi_v<float>;
                }
                sum_value += brdf * geom_term;
            }
            color += sum_value * al.area / scene->light_samples * al.color;
        }

    // indirect light
    float reflectivity = 0.0f;
    if (hit_shape->mat.brdf_type == BRDF_TYPE::Phong)
    {
        reflectivity = hit_shape->mat.specular.x + hit_shape->mat.specular.y + hit_shape->mat.specular.z;
        if (reflectivity + hit_shape->mat.diffuse.x + hit_shape->mat.diffuse.y + hit_shape->mat.diffuse.z <= EPS)
            reflectivity = 1.0f;
        else
            reflectivity = reflectivity / (reflectivity + hit_shape->mat.diffuse.x + hit_shape->mat.diffuse.y + hit_shape->mat.diffuse.z);
    }
    else if (hit_shape->mat.brdf_type == BRDF_TYPE::GGX)
    {
        float numerator = hit_shape->mat.specular.x + hit_shape->mat.specular.y + hit_shape->mat.specular.z;
        float denom = numerator + hit_shape->mat.diffuse.x + hit_shape->mat.diffuse.y + hit_shape->mat.diffuse.z;
        if (denom < EPS)
            reflectivity = 1.0f;
        else
            reflectivity = std::max(0.25f, numerator / denom);
    }
    bool sample_specular = uniform_dis(gen) <= reflectivity;
    vec3 bounced_ray_dirn{};
    if (hit_shape->mat.brdf_type == BRDF_TYPE::Phong)
    {
        vec3 center_vec = (scene->importance_sampling_type == vector_sampling_type::BRDF && sample_specular) ? reflected_ray_dir : rec.normal;
        bounced_ray_dirn = get_vector_around_normal(gen, uniform_dis, center_vec, scene->importance_sampling_type, sample_specular, hit_shape->mat.shininess);
    }
    else if (hit_shape->mat.brdf_type == BRDF_TYPE::GGX)
    {
        bounced_ray_dirn = get_brdf_unit_ggx(gen, uniform_dis, rec.normal, sample_specular, hit_shape->mat.roughness, r.direction * -1);
    }
    Ray bounced_ray{rec.hit_pos + rec.normal * EPS, bounced_ray_dirn};

    if (bounced_ray.direction == vec3{})
        return color;

    vec3 brdf{};
    if (hit_shape->mat.brdf_type == BRDF_TYPE::Phong)
    {
        brdf = hit_shape->mat.diffuse / std::numbers::pi + hit_shape->mat.specular * (
                         hit_shape->mat.shininess + 2) / (2 * std::numbers::pi) * std::pow(
                         dot(reflected_ray_dir, bounced_ray.direction), hit_shape->mat.shininess);
    }
    else if (hit_shape->mat.brdf_type == BRDF_TYPE::GGX)
    {
        vec3 half_vector = normalize(bounced_ray.direction - r.direction);
        float tan_sqrd_h = dot(rec.normal, half_vector) != 0.0f ? 1.0f / (std::pow(dot(rec.normal, half_vector), 2)) - 1.0f : 1000.0f;
        float alpha_sqrd = std::pow(hit_shape->mat.roughness, 2);
        float mic_dis_func = 0.0f;
        if (dot(rec.normal, half_vector) <= EPS || alpha_sqrd + tan_sqrd_h <= EPS)
            mic_dis_func = 1.0f;
        else
            mic_dis_func = alpha_sqrd / (std::numbers::pi_v<float> * std::pow(dot(rec.normal, half_vector), 4) * (std::pow(alpha_sqrd + tan_sqrd_h, 2)));
        float tan_sqrd_wi = dot(rec.normal, bounced_ray.direction) != 0 ? 1.0f / std::pow(dot(rec.normal, bounced_ray.direction), 2) - 1.0f : 1000.0f;
        float tan_sqrd_wo = dot(rec.normal, r.direction * -1.0f) != 0 ? 1.0f / std::pow(dot(rec.normal, r.direction * -1.0f), 2) - 1.0f : 1000.0f;
        float smith_g_func = (2.0f / (1.0f + std::sqrt(1.0f + alpha_sqrd * tan_sqrd_wi))) * (2.0f / (1.0f + std::sqrt(1.0f + alpha_sqrd * tan_sqrd_wo)));
        vec3 fresnel_term = hit_shape->mat.specular + (vec3{1.0f, 1.0f, 1.0f} - hit_shape->mat.specular) * std::pow(1 - dot(bounced_ray.direction, half_vector), 5);
        vec3 brdf_specular{};
        if (std::abs(dot(r.direction * -1.0f, rec.normal)) <= EPS || std::abs(dot(bounced_ray.direction, rec.normal)) <= EPS)
            brdf_specular = vec3{1.0f, 1.0f, 1.0f};
        else
            brdf_specular = (fresnel_term * smith_g_func * mic_dis_func) / (4.0f * dot(bounced_ray.direction, rec.normal) * dot(r.direction * -1.0f, rec.normal));
        brdf = brdf_specular + hit_shape->mat.diffuse / std::numbers::pi_v<float>;
    }
    vec3 throughput{};
    if (scene->importance_sampling_type == vector_sampling_type::UNIFORM_HEMISPHERE)
    {
        throughput = brdf * (2 * std::numbers::pi) * dot(rec.normal, bounced_ray.direction);
    }
    else if (scene->importance_sampling_type == vector_sampling_type::COSINE)
    {
        throughput = brdf * std::numbers::pi;
    }
    else if (scene->importance_sampling_type == vector_sampling_type::BRDF)
    {
        float pdf = 0.0f;
        if (hit_shape->mat.brdf_type == BRDF_TYPE::Phong)
        {
            pdf = (1.0f - reflectivity) * dot(rec.normal, bounced_ray.direction) / std::numbers::pi + reflectivity * (
                      hit_shape->mat.shininess + 1) / (2 * std::numbers::pi) * std::pow(
                      dot(reflected_ray_dir, bounced_ray.direction), hit_shape->mat.shininess);
        }
        else if (hit_shape->mat.brdf_type == BRDF_TYPE::GGX)
        {
            vec3 half_vector = normalize(r.direction * -1 + bounced_ray.direction);
            float tan_sqrd_h = dot(rec.normal, half_vector) != 0 ? 1.0f / (std::pow(dot(rec.normal, half_vector), 2)) - 1.0f : 1000.0f;
            float alpha_sqrd = std::pow(hit_shape->mat.roughness, 2);
            float mic_dis_func = alpha_sqrd / (std::numbers::pi_v<float> * std::pow(dot(rec.normal, half_vector), 4) * (std::pow(alpha_sqrd + tan_sqrd_h, 2)));

            pdf = (1.0f - reflectivity) * dot(rec.normal, bounced_ray.direction) / std::numbers::pi + reflectivity * mic_dis_func * dot(rec.normal, half_vector) / (4.0f * dot(half_vector, bounced_ray.direction));
        }
        throughput = brdf * dot(rec.normal, bounced_ray.direction) / pdf;
    }
    incoming_throughput *= throughput;
    float termination_probability = scene->russian_roulette
                                        ? 1.0f - std::min(1.0f, std::max(incoming_throughput.x, std::max(incoming_throughput.y, incoming_throughput.z)))
                                        : 0.0f;
    if (!scene->russian_roulette || uniform_dis(gen) > termination_probability)
    {
        vec3 traced_color = get_pixel_color_ggx(bounced_ray, depth - 1, incoming_throughput / (1 - termination_probability));
        color += traced_color * throughput / (1 - termination_probability);
    }
    return color;
}

// handle which Integrator to use
vec3 Renderer::get_pixel_color(const Ray &r, uint32_t depth)
{
    return get_pixel_color_ggx(r, depth);

    if (scene->integrator == "raytracer")
        return get_pixel_color_raytrace(r, depth);
    else if (scene->integrator == "analyticdirect")
        return get_pixel_color_analyticDirect(r);
    else if (scene->integrator == "direct")
        return get_pixel_color_direct(r);
    else if (scene->integrator == "pathtracer")
        return get_pixel_color_pathtrace(r, depth);
    return vec3{};
}
