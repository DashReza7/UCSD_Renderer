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
        std::cout << std::format("Progress: {:.2f}%\r", static_cast<float>(render_progress) / static_cast<float>(scene->height) * 100.0f);
    });
    std::cout << "Progress: 100.0%" << std::endl;
}

void Renderer::render_sequential()
{
    for (uint32_t i = 0; i < scene->height; ++i)
    {
        for (uint32_t j = 0; j < scene->width; ++j)
        {
            // if (i != 240 || j != 240)
            //     continue;

            vec3 pixel_color = vec3{};
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
        std::cout << std::format("Progress: {:.2f}%\r", static_cast<float>(i) / static_cast<float>(scene->height) * 100.0f);
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

bool Renderer::hit(const Ray &r, float t_min, float t_max, HitRecord &rec) const
{
    return scene->use_bvh ? scene->bvh->hit(r, T_MIN, T_MAX, rec) : hit_brute_force(r, T_MIN, T_MAX, rec);
}

// ray tracer Integrator
vec3 Renderer::get_pixel_color_raytrace(const Ray &r, uint32_t depth) const
{
    if (depth <= 0)
        return vec3{};

    HitRecord rec;
    bool is_hit = hit(r, T_MIN, T_MAX, rec);
    if (!is_hit)
        return vec3{};

    vec3 color{};

    const auto hit_shape = dynamic_cast<Shape *>(rec.hit_obj);

    color += hit_shape->mat.ambient;
    color += hit_shape->mat.emission;

    // vec3 eye_dirn = normalize(scene->main_camera->position - rec.hit_pos);

    HitRecord foo_rec;
    if (!scene->dirn_light.dirn_to_light.is_exactly_zero())
    {
        Ray shadow_ray{rec.hit_pos + rec.normal * EPS, scene->dirn_light.dirn_to_light};
        is_hit = hit(shadow_ray, T_MIN, T_MAX, foo_rec);
        if (!is_hit)
        {
            vec3 tmp_color = vec3{0.0f, 0.0f, 0.0f};

            tmp_color += hit_shape->mat.diffuse * fmax(0.0f, dot(rec.normal, scene->dirn_light.dirn_to_light));
            vec3 half_dirn = normalize(scene->dirn_light.dirn_to_light - r.direction);
            tmp_color +=
                    hit_shape->mat.specular * std::powf(fmax(0.0f, dot(half_dirn, rec.normal)), hit_shape->mat.shininess);

            color += tmp_color * scene->dirn_light.color;
        }
    }
    for (auto &pt_light: scene->pointlights)
    {
        float light_distance = norm(pt_light.position - rec.hit_pos);
        vec3 light_dirn = normalize(pt_light.position - rec.hit_pos);

        foo_rec.t = std::numeric_limits<float>::max();
        Ray shadow_ray{rec.hit_pos + rec.normal * EPS, normalize(pt_light.position - rec.hit_pos)};
        is_hit = hit(shadow_ray, T_MIN, T_MAX, foo_rec);
        if (is_hit && foo_rec.t < light_distance)
            continue;

        vec3 diffuse_color = hit_shape->mat.diffuse * fmax(0.0f, dot(light_dirn, rec.normal));

        vec3 half_dirn = normalize(light_dirn - r.direction);
        vec3 specular_color =
                hit_shape->mat.specular * std::powf(fmax(0.0f, dot(half_dirn, rec.normal)), hit_shape->mat.shininess);

        color += (diffuse_color + specular_color) * pt_light.color * pt_light.calc_attenuation(rec.hit_pos);
    }
    vec3 reflected_dirn = reflect(r.direction * -1, rec.normal);
    Ray reflected_ray = Ray{rec.hit_pos + rec.normal * EPS, reflected_dirn};
    if (hit_shape->mat.specular != vec3{0.0f, 0.0f, 0.0f})
        color += hit_shape->mat.specular * get_pixel_color_raytrace(reflected_ray, depth - 1);

    color.clamp(0.0f, 1.0f);
    return color;
}

// analytic direct light Integrator
vec3 Renderer::get_pixel_color_analyticDirect(const Ray &r) const
{
    HitRecord rec;
    bool is_hit = hit(r, T_MIN, T_MAX, rec);

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

        const auto hit_shape = dynamic_cast<Shape *>(rec.hit_obj);

        //            color += hit_shape->mat.ambient;
        //            color += hit_shape->mat.emission;

        for (const AreaLight &al: scene->areaLights)
            color += al.calculate_analytic_radiance(rec, hit_shape->mat.diffuse);

        color.clamp(0.0f, 1.0f);
        return color;
    }
}

// Monte-Carlo direct area light
vec3 Renderer::get_pixel_color_direct(const Ray &r)
{
    HitRecord rec;
    bool is_hit = hit(r, T_MIN, T_MAX, rec);

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
        auto hit_shape = dynamic_cast<Shape *>(rec.hit_obj);
        for (const AreaLight &al: scene->areaLights)
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
                is_hit = hit(shadow_ray, T_MIN, T_MAX, foo_rec);
                if (is_hit && foo_rec.t < light_dist)
                    continue;
                sum_value += (hit_shape->mat.diffuse / std::numbers::pi +
                              hit_shape->mat.specular / (2 * std::numbers::pi) * (hit_shape->mat.shininess + 2) *
                              std::powf(dot(reflected_ray_dir, surface_to_light_vec), hit_shape->mat.shininess)) *
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
    bool is_hit = hit(r, T_MIN, T_MAX, rec);
    vec3 tmp_vec;
    AreaLight *hitted_area_light = hit_area_light(r, tmp_vec);
    if (hitted_area_light != nullptr)
    {
        if (dot(hitted_area_light->normal, r.direction * -1.0f) <= 0.0f)
            return vec3{};

        if (depth == scene->maxdepth)
            return vec3{};

        // return the light color, if either no NEE, or the ray is in its first iteration (in other iterations, direct light has already been taken into account)
        if (scene->nee_type == NextEventEstimationType::OFF || depth == scene->maxdepth)
            return hitted_area_light->color;
        return vec3{};
    }
    if (!is_hit)
        return vec3{};

    auto hit_shape = dynamic_cast<Shape *>(rec.hit_obj);
    vec3 color{};

    // direct light
    if (scene->nee_type == NextEventEstimationType::ON)
    {
        for (const AreaLight& al: scene->areaLights)
        {
            // sample random points on the area_light source
            std::vector<vec3> rnd_light_points;
            if (!scene->light_stratify)
            {
                for (int i = 0; i < scene->light_samples; ++i)
                {
                    float rnd1 = uniform_dis(gen);
                    float rnd2 = uniform_dis(gen);
                    rnd_light_points.push_back(al.a + (al.b - al.a) * rnd1 + (al.c - al.a) * rnd2);
                }
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
                // surface_to_light_vec
                vec3 bounced_ray_dirn = rnd_light_points[i] - rec.hit_pos;
                float light_dist = norm(bounced_ray_dirn);
                bounced_ray_dirn = normalize(bounced_ray_dirn);
                // cast shadow ray
                HitRecord foo_rec;
                Ray shadow_ray{rec.hit_pos + rec.normal * EPS, bounced_ray_dirn};
                bool is_foo_hit = hit(shadow_ray, T_MIN, T_MAX, foo_rec);
                vec3 tmp_light_hit_pos{};
                AreaLight *hitted_area_light = hit_area_light(Ray{rec.hit_pos + rec.normal * EPS, bounced_ray_dirn}, tmp_light_hit_pos);
                if (&al != hitted_area_light || dot(hitted_area_light->normal, bounced_ray_dirn * -1.0f) >= 0.0f)
                    continue;
                if (is_foo_hit && foo_rec.t < light_dist) // if the area_light is in shadow, ignore this area_light source
                    continue;

                float cos_theta_i = std::clamp(dot(rec.normal, bounced_ray_dirn), 0.0f, 1.0f);
                float cos_theta_o = std::clamp(dot(al.normal, bounced_ray_dirn), 0.0f, 1.0f);
                float geom_term = cos_theta_i * cos_theta_o / dot(rec.hit_pos - rnd_light_points[i], rec.hit_pos - rnd_light_points[i]);
                vec3 brdf = get_brdf(hit_shape, r.direction * -1.0f, bounced_ray_dirn, rec.normal);
                sum_value += brdf * geom_term;
            }
            color += sum_value * al.area / scene->light_samples * al.color;
        }
    }
    else if (scene->nee_type == NextEventEstimationType::MIS)
    {
        vec3 sum_value{};

        // BRDF sample
        if (true)
        {
            vec3 sample_dirn_brdf = get_sample_dirn(ImportanceSamplingType::BRDF, uniform_dis(gen) < hit_shape->mat.get_reflectivity(), hit_shape, r.direction * -1.0f, rec.normal);
            vec3 light_hit_pos{};
            AreaLight *hitted_area_light = hit_area_light(Ray{rec.hit_pos + rec.normal * EPS, sample_dirn_brdf}, light_hit_pos);
            if (hitted_area_light != nullptr && dot(sample_dirn_brdf, hitted_area_light->normal) >= 0.0f)
            {
                vec3 brdf = sample_dirn_brdf.is_exactly_zero() ? vec3{} : get_brdf(hit_shape, r.direction * -1.0f, sample_dirn_brdf, rec.normal);
                float pdf_brdf = sample_dirn_brdf.is_exactly_zero() ? 1.0f : get_pdf(ImportanceSamplingType::BRDF, rec.normal, sample_dirn_brdf, r.direction * -1.0f, hit_shape);
                float weight = std::powf(pdf_brdf, 2) / (std::powf(pdf_brdf, 2) + std::powf(get_pdf_nee(sample_dirn_brdf, rec, is_hit), 2));
                // float cos_theta_i = std::clamp(dot(rec.normal, sample_dirn_brdf), 0.0f, 1.0f);
                // float cos_theta_o= std::clamp(dot(hitted_area_light->normal, sample_dirn_brdf), 0.0f, 1.0f);
                // float geom_term = cos_theta_i * cos_theta_o / dot(rec.hit_pos - light_hit_pos, rec.hit_pos - light_hit_pos);
                // sum_value += brdf / pdf_brdf * weight * geom_term * hitted_area_light->color;
                sum_value += brdf / pdf_brdf * weight * hitted_area_light->color * dot(rec.normal, sample_dirn_brdf);
            }
        }
        // NEE sample for each area light
        if (true)
        for (const AreaLight &al: scene->areaLights)
        {
            vec3 rnd_light_point{al.a + (al.b - al.a) * uniform_dis(gen) + (al.c - al.a) * uniform_dis(gen)};
            // surface_to_light_vec
            vec3 bounced_ray_dirn = rnd_light_point - rec.hit_pos;
            float light_dist = norm(bounced_ray_dirn);
            bounced_ray_dirn = normalize(bounced_ray_dirn);
            // cast shadow ray
            HitRecord foo_rec;
            Ray shadow_ray{rec.hit_pos + rec.normal * EPS, bounced_ray_dirn};
            bool is_foo_hit = hit(shadow_ray, T_MIN, T_MAX, foo_rec);
            vec3 tmp_hit_light_pos{};
            AreaLight *hitted_area_light = hit_area_light(Ray{rec.hit_pos + rec.normal * EPS, bounced_ray_dirn}, tmp_hit_light_pos);
            if (hitted_area_light != &al || dot(hitted_area_light->normal, bounced_ray_dirn * -1.0f) >= 0.0f)
                continue;
            if (!is_foo_hit || foo_rec.t >= light_dist) // if the area_light is in shadow, ignore this area_light source
            {
                vec3 brdf = get_brdf(hit_shape, r.direction * -1.0f, bounced_ray_dirn, rec.normal);
                float pdf_nee = (light_dist * light_dist / (al.area * std::abs(dot(al.normal, bounced_ray_dirn)))) / scene->areaLights.size();
                float pdf_brdf = get_pdf(ImportanceSamplingType::BRDF, rec.normal, bounced_ray_dirn, r.direction * -1.0f, hit_shape);
                float weight = std::powf(pdf_nee, 2) / (std::powf(pdf_nee, 2) + std::powf(pdf_brdf, 2));
                sum_value += brdf / pdf_nee * dot(rec.normal, bounced_ray_dirn) * weight * al.color / scene->areaLights.size();

                // float cos_theta_i = std::clamp(dot(rec.normal, bounced_ray_dirn), 0.0f, 1.0f);
                // float cos_theta_o= std::clamp(dot(al.normal, bounced_ray_dirn), 0.0f, 1.0f);
                // float geom_term = cos_theta_i * cos_theta_o / dot(rec.hit_pos - light_hit_pos, rec.hit_pos - light_hit_pos);
                // sum_value += brdf / pdf_nee * weight * al.color / scene->areaLights.size() * geom_term;
            }
        }

        color += sum_value;
    }

    // indirect light
    Ray bounced_ray{rec.hit_pos + rec.normal * EPS, get_sample_dirn(scene->importance_sampling_type, uniform_dis(gen) <= hit_shape->mat.get_reflectivity(), hit_shape, r.direction * -1.0f, rec.normal)};
    if (bounced_ray.direction == vec3{}) // if bounced_ray.direction is all zero, it means this new ray should have no contribution
        return color;
    vec3 brdf = get_brdf(hit_shape, r.direction * -1.0f, bounced_ray.direction, rec.normal);
    float pdf = get_pdf(scene->importance_sampling_type, rec.normal, bounced_ray.direction, r.direction * -1.0f, hit_shape);
    vec3 throughput = brdf / pdf * dot(rec.normal, bounced_ray.direction);
    float termination_prob = !scene->russian_roulette ? 0.0f : 1.0f - std::min(1.0f, max_vec(incoming_throughput * throughput));
    if (!scene->russian_roulette || uniform_dis(gen) > termination_prob)
    {
        vec3 traced_color = get_pixel_color_pathtrace(bounced_ray, depth - 1, incoming_throughput * throughput / (1.0f - termination_prob));
        color += traced_color * throughput / (1 - termination_prob);
    }

    return color;
}

// handle which Integrator to use
vec3 Renderer::get_pixel_color(const Ray &r, uint32_t depth)
{
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

vec3 Renderer::get_brdf(const Shape *shape, const vec3 &w_o, const vec3 &w_i, const vec3 &normal) const
{
    vec3 brdf{};
    if (shape->mat.brdf_type == BRDF_TYPE::Phong)
    {
        vec3 mirror_reflected_wo = reflect(w_o, normal);
        brdf = shape->mat.diffuse / std::numbers::pi + shape->mat.specular / (2 * std::numbers::pi) * (shape->mat.shininess + 2) * std::powf(dot(mirror_reflected_wo, w_i), shape->mat.shininess);
    }
    else if (shape->mat.brdf_type == BRDF_TYPE::GGX)
    {
        if (dot(w_i, normal) <= 0.0f || dot(w_o, normal) <= 0.0f)
            return vec3{};

        vec3 half_vector = normalize(w_i + w_o);
        // float tan_sqrd_h = is_almost_zero(dot(normal, half_vector)) ? std::numeric_limits<float>::infinity() : 1.0f / (std::powf(dot(normal, half_vector), 2)) - 1.0f;
        float cos_sqrd_h = dot(normal, half_vector) * dot(normal, half_vector);
        float alpha_sqrd = std::powf(shape->mat.roughness, 2);
        // float mic_dis_func = alpha_sqrd / (std::numbers::pi_v<float> * std::powf(dot(normal, half_vector), 4) * (std::powf(alpha_sqrd + tan_sqrd_h, 2)));
        float mic_dis_func = alpha_sqrd / (std::numbers::pi_v<float> * (alpha_sqrd * alpha_sqrd * cos_sqrd_h * cos_sqrd_h + (1.0f - cos_sqrd_h) * (1.0f - cos_sqrd_h) + 2.0f * alpha_sqrd * cos_sqrd_h * (1.0f - cos_sqrd_h)));
        float tan_sqrd_wi = is_almost_zero(dot(normal, w_i)) ? std::numeric_limits<float>::infinity() : 1.0f / std::powf(dot(normal, w_i), 2) - 1.0f;
        float tan_sqrd_wo = is_almost_zero(dot(normal, w_o)) ? std::numeric_limits<float>::infinity() : 1.0f / std::powf(dot(normal, w_o), 2) - 1.0f;
        float smith_g_func = (2.0f / (1.0f + std::sqrt(1.0f + alpha_sqrd * tan_sqrd_wi))) * (2.0f / (1.0f + std::sqrt(1.0f + alpha_sqrd * tan_sqrd_wo)));
        vec3 fresnel_term = shape->mat.specular + (vec3{1.0f, 1.0f, 1.0f} - shape->mat.specular) * std::powf(1 - dot(w_i, half_vector), 5);
        vec3 brdf_specular{};
        // if (std::abs(dot(w_o, normal)) <= EPS || std::abs(dot(w_i, normal)) <= EPS)
        //     brdf_specular = vec3{1.0f, 1.0f, 1.0f};
        // else
            brdf_specular = fresnel_term * smith_g_func * mic_dis_func / (4.0f * dot(w_i, normal) * dot(w_o, normal));
        brdf = brdf_specular + shape->mat.diffuse / std::numbers::pi_v<float>;
    }

    return brdf;
}

float Renderer::get_pdf(ImportanceSamplingType sampling_type, const vec3 &normal, const vec3 &w_i, const vec3 &w_o, const Shape *hit_shape) const
{
    float pdf = 0.0f;
    if (sampling_type == ImportanceSamplingType::UNIFORM_HEMISPHERE)
    {
        pdf = 1.0f / (2.0f * std::numbers::pi);
    }
    else if (sampling_type == ImportanceSamplingType::COSINE)
    {
        // TODO: division by this pdf might raise DivisionByZero Error!
        pdf = dot(normal, w_i) / std::numbers::pi;
    }
    else if (sampling_type == ImportanceSamplingType::BRDF)
    {
        float reflectivity = hit_shape->mat.get_reflectivity();

        if (hit_shape->mat.brdf_type == BRDF_TYPE::Phong)
        {
            pdf = (1.0f - reflectivity) * dot(normal, w_i) / std::numbers::pi_v<float> + reflectivity * (hit_shape->mat.shininess + 1) / (2.0f * std::numbers::pi_v<float>) * std::powf(dot(w_i, reflect(w_o, normal)), hit_shape->mat.shininess);
        }
        else if (hit_shape->mat.brdf_type == BRDF_TYPE::GGX)
        {
            vec3 half_vector = normalize(w_o + w_i);
            float cos_sqrd_h = dot(normal, half_vector) * dot(normal, half_vector);
            float alpha_sqrd = std::powf(hit_shape->mat.roughness, 2);
            float mic_dis_func = alpha_sqrd / (std::numbers::pi_v<float> * (alpha_sqrd * alpha_sqrd * cos_sqrd_h * cos_sqrd_h + (1.0f - cos_sqrd_h) * (1.0f - cos_sqrd_h) + 2.0f * alpha_sqrd * cos_sqrd_h * (1.0f - cos_sqrd_h)));
    
            pdf = (1.0f - reflectivity) * dot(normal, w_i) / std::numbers::pi + reflectivity * mic_dis_func * dot(normal, half_vector) / (4.0f * dot(half_vector, w_i));
        }
    }

    return pdf;
}

vec3 Renderer::get_sample_dirn(ImportanceSamplingType sampling_type, bool sample_specular, const Shape *hit_shape, const vec3 &w_o, const vec3 &normal)
{
    vec3 bounced_ray_dirn{};
    if (sampling_type == ImportanceSamplingType::UNIFORM_HEMISPHERE)
    {
        // in the range [0, pi]
        //float phi = std::numbers::pi * uniform_dis(generator);
        float phi = std::acos(uniform_dis(gen));
        // in the range [0, 2*pi]
        float theta = 2.0f * std::numbers::pi * uniform_dis(gen);
        bounced_ray_dirn = align_vector(vec3{std::sin(phi) * std::cos(theta), std::sin(phi) * std::sin(theta), std::abs(std::cos(phi))}, normal);
    }
    else if (sampling_type == ImportanceSamplingType::COSINE)
    {
        float phi = 2 * std::numbers::pi * uniform_dis(gen);
        float theta = std::acos(std::sqrt(uniform_dis(gen)));
        bounced_ray_dirn = align_vector(vec3{std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta), std::abs(std::cos(theta))}, normal);
    }
    else if (sampling_type == ImportanceSamplingType::BRDF)
    {
        if (hit_shape->mat.brdf_type == BRDF_TYPE::Phong)
        {
            float phi = 2.0f * std::numbers::pi_v<float> * uniform_dis(gen);
            float theta = sample_specular ? std::acos(std::powf(uniform_dis(gen), 1.0f / (hit_shape->mat.shininess + 1.0f))) : std::acos(std::sqrtf(uniform_dis(gen)));
            vec3 tmp_dirn = vec3{std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta), std::cos(theta)};
            if (!sample_specular || tmp_dirn.z >= 0.0f)
                bounced_ray_dirn = sample_specular ? align_vector(tmp_dirn, reflect(w_o, normal)) : align_vector(tmp_dirn, normal);
        }
        else if (hit_shape->mat.brdf_type == BRDF_TYPE::GGX)
        {
            if (sample_specular)
            {
                float phi = 2.0f * std::numbers::pi_v<float> * uniform_dis(gen);
                float rnd = uniform_dis(gen);
                float theta = std::atan(hit_shape->mat.roughness * std::sqrtf(rnd) / std::sqrtf(1.0f - rnd));
                vec3 h{std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta), std::cos(theta)};
                vec3 h_rotated = align_vector(h, normal);
                vec3 wi = reflect(w_o, h_rotated);
                if (dot(wi, normal) >= 0.0f)
                    bounced_ray_dirn = wi;
            }
            else
            {
                float phi = 2 * std::numbers::pi_v<float> * uniform_dis(gen);
                float theta = std::acos(std::sqrtf(uniform_dis(gen)));
                vec3 sample{std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta), std::cos(theta)};
                bounced_ray_dirn = align_vector(sample, normal);
            }
        }
    }

    return bounced_ray_dirn;
}

float Renderer::get_pdf_nee(const vec3 &w_i, const HitRecord &world_rec, bool is_world_hit) const
{
    vec3 light_hit_pos{};
    AreaLight *hitted_area_light = hit_area_light(Ray{world_rec.hit_pos + world_rec.normal * EPS, w_i}, light_hit_pos);
    if (hitted_area_light == nullptr || dot(hitted_area_light->normal, w_i * -1.0f) >= 0.0f)
        return 0.0f;
    return (norm2(world_rec.hit_pos - light_hit_pos) / (hitted_area_light->area * std::abs(dot(hitted_area_light->normal, w_i)))) / scene->areaLights.size();
}

AreaLight *Renderer::hit_area_light(const Ray &r, vec3 &light_hit_pos) const
{
    HitRecord world_rec;
    bool is_world_hit = hit(r, T_MIN, T_MAX, world_rec);
    bool is_hit_area_light = false;
    float t_area_light = std::numeric_limits<float>::max();
    AreaLight *hitted_area_light = nullptr;
    for (auto &al: scene->areaLights)
    {
        HitRecord foo_rec;
        vec3 foo_color;
        if (al.hit(r, foo_rec, foo_color) == true && foo_rec.t < t_area_light)
        {
            is_hit_area_light = true;
            t_area_light = foo_rec.t;
            hitted_area_light = &al;
        }
    }

    if (!is_hit_area_light || (is_world_hit && t_area_light > world_rec.t))
        return nullptr;
    light_hit_pos = r.origin + r.direction * t_area_light;
    return hitted_area_light;
}
