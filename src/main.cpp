#include <iostream>
#include <chrono>
#include <cstdlib>
#include "Scene.h"
#include "Renderer.h"

#define ARCHIVE_IMAGES true

void run(int argc, char *argv[])
{
    if (argc < 2)
        throw std::runtime_error("no source scene files provided");
    
    double total_render_time = 0.0;
    for (int i = 1; i < argc; ++i)
    {
        Scene main_scene = Scene{};
        const char *input_filename = argv[i];
        std::string output_filename;
        main_scene.parse_scene_file(input_filename, output_filename);
        
        main_scene.init();
        
        Film film{main_scene.width, main_scene.height};
        
        Renderer renderer{&main_scene, &film};
        
        std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
        renderer.render();
        std::chrono::steady_clock::time_point finish_time = std::chrono::steady_clock::now();
        double duration = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds >(finish_time - start_time).count());
        total_render_time += duration;
        std::cout << input_filename << " -> " << output_filename << "  duration: " << duration / 1000.0f << " seconds" << std::endl;
        
        film.write_to_file_png(("../ResultImages/" + output_filename).c_str());
        
        main_scene.terminate();
    }
    std::cout << "\ntotal render time: " << total_render_time / 1000.0 << std::endl << std::endl;
    
    if (ARCHIVE_IMAGES)
    {
        const char* command = "powershell.exe -ExecutionPolicy Bypass -File C:/Users/masou/Desktop/Codes/UCSD_Renderer/libs/compress_images.ps1";
        int result = system(command);
        
        if (result != 0)
            std::cout << "powershell command failed to execute" << std::endl;
    }
}

int main(int argc, char *argv[])
{
    try
    {
        run(argc, argv);
    } catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    
    // std::cout << "\nPress Enter to close the console...";
    // std::cin.get();
    return 0;
}
