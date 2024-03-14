#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "geometry.h"
// #define STB_IMAGE_IMPLEMENTATION
// #include "stb_image.h"
// #define STB_IMAGE_WRITE_IMPLEMENTATION
// #include "stb_image_write.h"

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

struct Material {
    Vec3f diffuse_color;
    Vec2f albedo; //diffuse_coff, spec_coff
    float spec_expo;
    Material() : diffuse_color() {}
    Material(const Vec2f &a, const Vec3f& color, float spec) : diffuse_color{color}, albedo{a}, spec_expo{spec} {}
};

struct Light {
    Vec3f pos;
    float intensity;
    Light(const Vec3f& p, const float i) : pos{p}, intensity{i} {}
};

struct Hit {
    Vec3f o, dir;
    int isHit;
    int numberOfHits;
    float t0; //p + t0*dir = the first hit point
    Hit(const Vec3f& orig, const Vec3f& dir) : o{orig}, dir{dir} {}
    Vec3f hitPoint() const {
        return o + dir.normalize() * t0;
    }
};

struct Sphere {
    Vec3f center;
    float radius;
    Material material;

    Sphere(const Vec3f& c, const float r, const Material& mat) : center{c}, radius{r}, material(mat) {}

    bool ray_intersect(const Vec3f& o, const Vec3f& dir, Hit& hit) {
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
            hit.t0 = a - sqrt(radius*radius-b_2);
            hit.isHit = 1;
            hit.numberOfHits = 2;
            return true;
        }
        else if (b - radius < 1e-4) {
            //1 intersection
            hit.t0 = a;
            hit.isHit = 1;
            hit.numberOfHits = 1;
            return true;
        }
        //no intersection
        hit.isHit = 0;
        return false;
    }
};

Vec3f shade(
    const std::vector<Light>& lights, 
    const Sphere& sphere,
    const Hit& hit) {
    
    Vec3f p = hit.hitPoint();
    Vec3f n = (p - sphere.center).normalize();
    float diffuse_intensity = 0.0f;
    float spec_intensity = 0.0f;
    for (auto& light : lights) {
        Vec3f l = (light.pos - p).normalize();
        Vec3f r = (n*2*(n*l) - l).normalize();
        Vec3f v = (hit.o - p).normalize();
        diffuse_intensity += light.intensity * std::max(0.0f, n * l);
        spec_intensity += light.intensity * powf(std::max(0.0f, v*r), sphere.material.spec_expo);
    }

    Vec3f diffuse = sphere.material.diffuse_color * diffuse_intensity * sphere.material.albedo.x;
    Vec3f specular = Vec3f(1.f, 1.f, 1.f) * spec_intensity * sphere.material.albedo.y;
    return diffuse + specular;
}

void render() {
    const int width = 1024;
    const int height = 768;
    std::vector<Vec3f> framebuffer(width*height);
    std::fill(framebuffer.begin(), framebuffer.end(), Vec3f(0.2, 0.7, 0.8));

    float aspect = float(height) / float(width);
    float fov = M_PI_2;
    float fov_2 = M_PI_4;
    //scene
    std::vector<Light> lights = {
        Light(Vec3f(-20, 20,  20), 1.5),
        Light(Vec3f( 30, 50, -25), 1.8),
        Light(Vec3f( 30, 20,  30), 1.7)
    };

    Material ivory({0.6,  0.3}, Vec3f(0.4, 0.4, 0.3), 50.0f);
    Material red_rubber({0.9,  0.1}, Vec3f(0.3, 0.1, 0.1), 10.0f);
    std::vector<Sphere> spheres = {
        Sphere(Vec3f(-3,    0,   -16), 2,      ivory),
        Sphere(Vec3f(-1.0, -1.5, -12), 2, red_rubber),
        Sphere(Vec3f( 1.5, -0.5, -18), 3, red_rubber),
        Sphere(Vec3f( 7,    5,   -18), 4,      ivory)
    };

    float tan_fov_2 = tanf(fov_2);
    for (size_t y = 0; y < height; y++) {
        for (size_t x = 0; x < width; x++) {
            float wx = (x - width/2.f)*tan_fov_2/(width/2.0f);
            float wy = (y - height/2.f)*tan_fov_2*aspect/(height/2.0f);

            float nearest_factor = std::numeric_limits<float>::max();
            for (Sphere& s : spheres) {
                Hit hit({0.f,0.f,0.f}, {wx, wy, -1.0f});
                bool ishit = s.ray_intersect({0.0f, 0.0f, 0.0f}, {wx, wy, -1.0f}, hit);
                if (ishit && hit.t0 < nearest_factor) {
                    nearest_factor = hit.t0;
                    framebuffer[x+(height-1-y)*width] = shade(lights, s, hit);
                }
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