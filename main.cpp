#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "geometry.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

void write_ppm(
    std::vector<Vec3f> &framebuffer, 
    const int width, const int height, 
    const char* filename)
{
    std::ofstream ofs;
    ofs.open(filename);
    ofs << "P6\n"
        << width << " " << height << "\n255\n";
    for (size_t i = 0; i < width * height; ++i)
    {
        for (size_t j = 0; j < 3; j++)
        {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
}

struct Sphere {
    Vec3f center;
    float radius;

    Sphere(const Vec3f& c, const float r) : center{c}, radius{r} {}

    bool ray_intersect(const Vec3f& o, const Vec3f& dir, float& t0) {
        bool inside_sphere = (o-center).norm() < radius;
        if (inside_sphere) return false;
        auto oc = center - o;
        float k = oc*dir.normalize();
        bool ahead_sphere = k < 0;
        if (ahead_sphere) return false;
        auto a = oc*dir.normalize();
        auto c = oc.norm();
        auto b_2 = c*c-a*a;
        auto b = sqrt(b_2);
        if (b < radius) {
            //2 intersections
            t0 = a - sqrt(radius*radius-b_2);
            return true;
        }
        else if (b - radius < 1e-4) {
            //1 intersection
            t0 = a;
            return true;
        }
        //no intersection
        return false;
    }
};

void render() {
    const int width = 1024;
    const int height = 768;
    std::vector<Vec3f> framebuffer(width*height);
    std::fill(framebuffer.begin(), framebuffer.end(), Vec3f(0.2, 0.7, 0.8));

    float aspect = float(height) / float(width);
    float fov = M_PI_2;
    float fov_2 = M_PI_4;
    //scene
    Sphere s({0, 0, -15.5f}, 1.5f);

    float tan_fov_2 = tanf(fov_2);
    for (size_t y = 0; y < height; y++) {
        for (size_t x = 0; x < width; x++) {
            float wx = (x - width/2.f)*tan_fov_2/(width/2.0f);
            float wy = (y - height/2.f)*tan_fov_2*aspect/(height/2.0f);
            float t0;
            bool hit = s.ray_intersect({0.0f, 0.0f, 0.0f}, {wx, wy, -1.0f}, t0);
            if (hit) {
                framebuffer[x+y*width] = Vec3f(0.4, 0.4, 0.3);
            }
        }
    }
    write_ppm(framebuffer, width, height, "./out.ppm");
}


int main()
{
    render();
    return 0;
}