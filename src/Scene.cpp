#include <iostream>
#include <fstream>
#include <sstream>
#include <format>
#include <stack>
#include "Scene.h"

void Scene::parse_scene_file(const char *input_filename, std::string &output_filename) {
    std::ifstream file{input_filename};
    if (!file)
        throw std::runtime_error(std::format("Error opening scene file: {}", input_filename));

    // set defaults
    width = 640;
    height = 480;
    maxdepth = 5;
    output_filename = "raytrace.png";
    float atten_const = 1.0f;
    float atten_linear = 0.0f;
    float atten_quadratic = 0.0f;
    int vertex_idx = 0;
    std::vector<vec3> vertexes;
    int vertnormal_idx = 0;
    std::vector<std::pair<vec3, vec3> > vertnormals;
    // TODO: Should I make this zero or not?
    vec3 ambient = vec3{0.2f, 0.2f, 0.2f};
    vec3 diffuse = vec3{0.0f, 0.0f, 0.0f};
    vec3 specular = vec3{0.0f, 0.0f, 0.0f};
    float shininess = 1.0f;
    vec3 emission = vec3{0.0f, 0.0f, 0.0f};
    // allow only one directional light
    bool is_dirn_light_set = false;
    this->integrator = "raytracer";

    // keep track of transformations
    mat4 cur_transform = mat4{1.0f};
    mat4 cur_inv_transform = mat4{1.0f};
    std::stack<mat4> transform_stack;
    std::stack<mat4> inv_transform_stack;

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss{line};
        std::vector<std::string> tokens;

        std::string token;
        while (iss >> token)
            tokens.push_back(token);

        if (tokens.empty() || tokens[0][0] == '#')
            continue;

        std::string command = tokens[0];
        if (command == "size") {
            width = std::stoi(tokens[1]);
            height = std::stoi(tokens[2]);
        }
        else if (command == "maxdepth") {
            maxdepth = std::stoi(tokens[1]);
        }
        else if (command == "output") {
            output_filename = tokens[1];
        }
        else if (command == "camera") {
            vec3 look_from = vec3{std::stof(tokens[1]), std::stof(tokens[2]), std::stof(tokens[3])};
            vec3 look_at = vec3{std::stof(tokens[4]), std::stof(tokens[5]), std::stof(tokens[6])};
            vec3 up = normalize(vec3{std::stof(tokens[7]), std::stof(tokens[8]), std::stof(tokens[9])});
            float fov = std::stof(tokens[10]);
            main_camera = new Camera{
                fov, static_cast<float>(width) / static_cast<float>(height), up, look_from,
                look_at
            };
        }
        else if (command == "sphere") {
            vec3 center{std::stof(tokens[1]), std::stof(tokens[2]), std::stof(tokens[3])};
            float radius = std::stof(tokens[4]);
            spheres.emplace_back(center, radius, Material{diffuse, specular, shininess, ambient, emission},
                                 cur_transform, cur_inv_transform);
        }
        else if (command == "maxverts") {
            vertexes = std::vector<vec3>{std::stoul(tokens[1])};
        }
        else if (command == "maxvertnormals") {
            vertnormals = std::vector<std::pair<vec3, vec3> >{std::stoul(tokens[1])};
        }
        else if (command == "vertex") {
            vertexes[vertex_idx++] = vec3{std::stof(tokens[1]), std::stof(tokens[2]), std::stof(tokens[3])};
        }
        else if (command == "vertexnormal") {
            vertnormals[vertnormal_idx++] = std::pair<vec3, vec3>{
                vec3{std::stof(tokens[1]), std::stof(tokens[2]), std::stof(tokens[3])},
                vec3{std::stof(tokens[4]), std::stof(tokens[5]), std::stof(tokens[6])}
            };
        }
        else if (command == "tri") {
            vec4 tr_v1 = cur_transform * vec4{vertexes[std::stoi(tokens[1])], 1.0f};
            vec4 tr_v2 = cur_transform * vec4{vertexes[std::stoi(tokens[2])], 1.0f};
            vec4 tr_v3 = cur_transform * vec4{vertexes[std::stoi(tokens[3])], 1.0f};
            triangles.emplace_back(vec3{tr_v1.x, tr_v1.y, tr_v1.z}, vec3{tr_v2.x, tr_v2.y, tr_v2.z},
                                   vec3{tr_v3.x, tr_v3.y, tr_v3.z},
                                   Material{diffuse, specular, shininess, ambient, emission});
        }
        else if (command == "trinormal") {
            throw std::runtime_error("trinormal not implemented yet!");
        }
        else if (command == "translate") {
            vec3 direction = vec3{std::stof(tokens[1]), std::stof(tokens[2]), std::stof(tokens[3])};
            cur_transform = cur_transform * get_translation_matrix(direction);
            cur_inv_transform = get_translation_matrix(direction * -1) * cur_inv_transform;
        }
        else if (command == "rotate") {
            vec3 axis = normalize(vec3{std::stof(tokens[1]), std::stof(tokens[2]), std::stof(tokens[3])});
            float angle = std::stof(tokens[4]);
            while (angle >= 360.0f)
                angle -= 360.0f;
            while (angle <= -360.0f)
                angle += 360.0f;

            cur_transform = cur_transform * get_rotation_matrix(axis, radians(angle));
            cur_inv_transform = get_rotation_matrix(axis, radians(-angle)) * cur_inv_transform;
        }
        else if (command == "scale") {
            vec3 scale_vec = vec3{std::stof(tokens[1]), std::stof(tokens[2]), std::stof(tokens[3])};
            cur_transform = cur_transform * mat4{vec4{scale_vec, 1.0f}};
            cur_inv_transform =
                    mat4{vec4{1.0f / scale_vec[0], 1.0f / scale_vec[1], 1.0f / scale_vec[2], 1.0f}} * cur_inv_transform;
        }
        else if (command == "pushTransform") {
            transform_stack.push(cur_transform);
            inv_transform_stack.push(cur_inv_transform);
        }
        else if (command == "popTransform") {
            cur_transform = transform_stack.top();
            transform_stack.pop();

            cur_inv_transform = inv_transform_stack.top();
            inv_transform_stack.pop();
        }
        else if (command == "directional") {
            if (is_dirn_light_set)
                throw std::runtime_error("multiple directional lights not implemented yet.");

            vec3 direction_to_light = vec3{std::stof(tokens[1]), std::stof(tokens[2]), std::stof(tokens[3])};
            direction_to_light = normalize(cur_transform * direction_to_light);

            vec3 color = vec3{std::stof(tokens[4]), std::stof(tokens[5]), std::stof(tokens[6])};
            dirn_light = DirnLight{color, direction_to_light};

            is_dirn_light_set = true;
        }
        else if (command == "point") {
            vec4 tmp_position =
                    cur_transform * vec4{std::stof(tokens[1]), std::stof(tokens[2]), std::stof(tokens[3]), 1.0f};
            vec3 position = vec3{
                tmp_position.x / tmp_position.w, tmp_position.y / tmp_position.w, tmp_position.z / tmp_position.w
            };
            vec3 color = vec3{std::stof(tokens[4]), std::stof(tokens[5]), std::stof(tokens[6])};
            pointlights.emplace_back(color, position, atten_const, atten_linear, atten_quadratic);
        }
        else if (command == "attenuation") {
            atten_const = std::stof(tokens[1]);
            atten_linear = std::stof(tokens[2]);
            atten_quadratic = std::stof(tokens[3]);
        }
        else if (command == "ambient") {
            ambient.x = std::stof(tokens[1]);
            ambient.y = std::stof(tokens[2]);
            ambient.z = std::stof(tokens[3]);
        }
        else if (command == "diffuse") {
            diffuse.x = std::stof(tokens[1]);
            diffuse.y = std::stof(tokens[2]);
            diffuse.z = std::stof(tokens[3]);
        }
        else if (command == "specular") {
            specular.x = std::stof(tokens[1]);
            specular.y = std::stof(tokens[2]);
            specular.z = std::stof(tokens[3]);
        }
        else if (command == "shininess") {
            shininess = std::stof(tokens[1]);
        }
        else if (command == "emission") {
            emission.x = std::stof(tokens[1]);
            emission.y = std::stof(tokens[2]);
            emission.z = std::stof(tokens[3]);
        }
        else if (command == "use_bvh") {
            if (tokens[1] == "true")
                this->use_bvh = true;
            else if (tokens[1] == "false")
                this->use_bvh = false;
        }
        else if (command == "integrator") {
            this->integrator = tokens[1];
        }
        else if (command == "quadLight") {
            vec4 a_temp = vec4{std::stof(tokens[1]), std::stof(tokens[2]), std::stof(tokens[3]), 1.0f};
            vec4 b_temp = vec4{
                std::stof(tokens[4]) + a_temp.x, std::stof(tokens[5]) + a_temp.y, std::stof(tokens[6]) + a_temp.z, 1.0f
            };
            vec4 c_temp = vec4{
                std::stof(tokens[7]) + a_temp.x, std::stof(tokens[8]) + a_temp.y, std::stof(tokens[9]) + a_temp.z, 1.0f
            };
            vec4 d_temp = vec4{
                -a_temp.x + b_temp.x + c_temp.x, -a_temp.y + b_temp.y + c_temp.y, -a_temp.z + b_temp.z + c_temp.z, 1.0f
            };
            a_temp = cur_transform * a_temp;
            b_temp = cur_transform * b_temp;
            c_temp = cur_transform * c_temp;
            d_temp = cur_transform * d_temp;
            vec3 a = vec3{a_temp.x / a_temp.w, a_temp.y / a_temp.w, a_temp.z / a_temp.w};
            vec3 b = vec3{b_temp.x / b_temp.w, b_temp.y / b_temp.w, b_temp.z / b_temp.w};
            vec3 c = vec3{c_temp.x / c_temp.w, c_temp.y / c_temp.w, c_temp.z / c_temp.w};
            vec3 d = vec3{d_temp.x / d_temp.w, d_temp.y / d_temp.w, d_temp.z / d_temp.w};

            areaLights.emplace_back(vec3{std::stof(tokens[10]), std::stof(tokens[11]), std::stof(tokens[12])}, a, b, c,
                                    d);
        }
        else if (command == "lightsamples") {
            this->light_samples = std::stoi(tokens[1]);
        }
        else if (command == "lightstratify") {
            if (tokens[1] == "on")
                this->light_stratify = true;
            else if (tokens[1] == "off")
                this->light_stratify = false;
        }
        else if (command == "spp") {
            this->spp = std::stoi(tokens[1]);
        }
        else if (command == "nexteventestimation") {
            if (tokens[1] == "on")
                this->next_event_estimation = true;
            else if (tokens[1] == "off")
                this->next_event_estimation = false;
        }
        else if (command == "russianroulette") {
            if (tokens[1] == "on")
                this->russian_roulette = true;
            else if (tokens[1] == "off")
                this->russian_roulette = false;
        }
        else if (command == "parallel_run") {
            if (tokens[1] == "true")
                this->parallel_run = true;
            else if (tokens[1] == "false")
                this->parallel_run = false;
        }
        else if (command == "importancesampling") {
            if (tokens[1] == "hemisphere")
                this->importance_sampling_type = "hemisphere";
            else if (tokens[1] == "cosine")
                this->importance_sampling_type = "cosine";
        }
        else
            throw std::runtime_error(std::format("Invalid command: {}", command));
    }

    file.close();
}
