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

#include "model.h"

int envmap_width, envmap_height;
std::vector<Vec3f> envmap;

Model *mesh = nullptr;

const Vec3f background = Vec3f(0.2, 0.7, 0.8);
const int MAX_TRACE_DEPTH = 4;

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
    Vec4f albedo; //diffuse_coff, spec_coff, reflect_coff, refract_coff
    float spec_expo;
    float ridx;//refrection factor
    Material() : diffuse_color() {}
    Material(const Vec4f &a, const Vec3f& color, float spec, float ridx = 1.0f) 
    : diffuse_color{color}, albedo{a}, spec_expo{spec}, ridx{ridx} {}
};

Material ivory({0.6,  0.3, 0.1, 0.0f}, Vec3f(0.4, 0.4, 0.3), 50.0f);
Material red_rubber({0.9,  0.1, 0.0, 0.0f}, Vec3f(0.3, 0.1, 0.1), 10.0f);
Material mirror({0.0, 10.0, 0.8, 0.0f}, Vec3f(1.0, 1.0, 1.0), 1425.);
Material glass({0.0,  0.5, 0.1, 0.8}, Vec3f(0.6, 0.7, 0.8),  125., 1.5);

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
    Material mat;
    Vec3f normal;
    Hit(const Vec3f& orig, const Vec3f& dir) : o{orig}, dir{dir}, isHit{0}, numberOfHits{0}, t0{0} {}
    Vec3f hitPoint() const {
        return o + dir.normalize() * t0;
    }
};

struct Sphere {
    Vec3f center;
    float radius;
    Material material;

    Sphere(const Vec3f& c, const float r, const Material& mat) : center{c}, radius{r}, material(mat) {}

    bool ray_intersect(const Vec3f& o, const Vec3f& dir, Hit& hit) const {
        bool inside_sphere = (o-center).norm() < radius;
        auto oc = center - o;
        float k = oc*dir.normalize();
        bool ahead_sphere = k < 0;
        if (ahead_sphere && !inside_sphere) return false;
        auto a = oc*dir.normalize();
        auto c = oc.norm();
        auto b_2 = c*c-a*a;
        auto b = sqrt(b_2);
        if (b < radius) {
            //maybe 2 intersections
            auto d = sqrt(radius*radius-b_2);
            hit.t0 = inside_sphere ? a + d : a - d;
            hit.isHit = 1;
            hit.numberOfHits = inside_sphere ? 1 : 2;
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

bool scene_intersect(
    const Vec3f& o, 
    const Vec3f& dir, 
    const std::vector<Sphere>& spheres, 
    Hit& hit, Sphere const ** hitter) {
    
    float nearest_factor = std::numeric_limits<float>::max();
    for (auto& ss : spheres) {
        Hit h(o, dir);
        bool ishit = ss.ray_intersect(o, dir, h);
        if (ishit && h.t0 < nearest_factor) {
            nearest_factor = h.t0;
            *hitter = &ss;
            hit = h;
        }
    }

    auto nfaces = mesh->nfaces();
    for (int fid = 0; fid < nfaces; fid++) {
        float tnear;
        if (mesh->ray_triangle_intersect(fid, o, dir, tnear) && tnear < nearest_factor) {
            nearest_factor = tnear;
            *hitter = nullptr;
            hit.isHit = 1;
            hit.t0 = tnear;
            hit.mat = mirror;
            auto v0 = mesh->point(mesh->vert(fid, 0));
            auto v1 = mesh->point(mesh->vert(fid, 1));
            auto v2 = mesh->point(mesh->vert(fid, 2));
            auto n = cross((v1-v0), (v2-v0)).normalize();
            hit.normal = n;
        }
    }

    //define a plane at y = -4
    float d = (-4 - o.y) / dir.y;
    auto p = o + dir*d;
    if (d > 0.0f && ((hit.isHit && d < nearest_factor) || !hit.isHit) && p.x >= -6 && p.x <=6 && p.z >= -28 && p.z <= -10) {
        hit.isHit = 1;
        hit.t0 = d;
        hit.mat = Material({0.3f, 0.0f, 0.0f, 0.0f}, ((int(p.x*0.5 + 101) + int(p.z*0.5)) & 1) ? Vec3f{0.9, 0.9, 0.9} : Vec3f{0.7, 0.7, 0.7}, 1.0f);
        hit.normal = Vec3f{0.0f, 1.0f, 0.0f};
        *hitter = nullptr;
    }

    return hit.isHit == 1;
}

Vec3f reflect(const Vec3f& l, const Vec3f& n) {
    return n * 2.0f * (l * n) - l;
}

Vec3f refrect(const Vec3f& l, const Vec3f& n, float n1, float n2) {
    auto r = n1 / n2;
    auto c = (-n*l);
    Vec3f nn = n;
    if (c < 0) {
        nn = nn * -1.0f;
        r = n2/n1;
        c = (-nn*l);
    }
    return (l*r+nn*(r*c-sqrt(1.0f - r*r*(1.0f - c*c)))).normalize();
}

Vec3f ray_trace(
    const std::vector<Light>& lights,
    const std::vector<Sphere>& spheres,
    const Vec3f& origin,
    const Vec3f& dir,
    int depth = 1) {

    if (depth > MAX_TRACE_DEPTH) return background;
    
    Hit hit(origin, dir);
    Sphere const * s = nullptr;
    bool ishit = scene_intersect(origin, dir.normalize(), spheres, hit, &s);
    if (ishit) {
        
        Vec3f p = hit.hitPoint();
        Vec3f n = (s != nullptr) ? (p - s->center).normalize() : hit.normal;
        Material m = (s != nullptr) ? s->material : hit.mat;
        float diffuse_intensity = 0.0f;
        float spec_intensity = 0.0f;
        for (auto& light : lights) {
            Vec3f l = (light.pos - p).normalize();
            Vec3f r = (n*2*(n*l) - l).normalize();
            Vec3f v = (hit.o - p).normalize();
        
            //shadow
            //shoot a ray from intersection poitn p to the light source
            //see if there is any thing between the point and the light source.
            Hit h(p + l*0.0001, l);
            Sphere const * hitter = nullptr;
            if (scene_intersect(h.o, h.dir, spheres, h, &hitter)) {
                continue;
            }

            diffuse_intensity += light.intensity * std::max(0.0f, n * l);
            spec_intensity += light.intensity * powf(std::max(0.0f, v*r), m.spec_expo);
        }

        Vec3f diffuse =m.diffuse_color * diffuse_intensity * m.albedo.x;
        Vec3f specular = Vec3f(1.f, 1.f, 1.f) * spec_intensity * m.albedo.y;
        Vec3f reflection = ray_trace(lights, spheres, p+n*0.001, reflect((origin-p), n).normalize(), depth+1) * m.albedo.z;
        //refraction
        Vec3f refrection = {0.0f, 0.0f, 0.0f};
        if (m.albedo[3] > 0.f) {
            auto r1 = refrect(dir.normalize(), n, 1.0f, m.ridx);
            refrection = ray_trace(lights, spheres, p+r1*0.0001, r1, depth+1) * m.albedo[3];
        }

        return diffuse + specular + reflection + refrection;
    }

    auto ndir = dir.normalize();
    float theta = std::atan2(ndir.z, ndir.x);
    float phi = std::acos(ndir.y);
    float u = (theta + M_PI) / (2.0f * M_PI);
    float v = phi / M_PI;
    return envmap[int(u*envmap_width)+int(v*envmap_height)*(envmap_width)];
    
}

void render() {
    const int width = 1024;
    const int height = 768;
    std::vector<Vec3f> framebuffer(width*height);
    std::fill(framebuffer.begin(), framebuffer.end(), Vec3f(0.2, 0.7, 0.8));

    float aspect = float(height) / float(width);
    float fov = M_PI_2;
    float fov_2 = M_PI_4;

    int n;
    unsigned char *pixmap = stbi_load("./envmap.jpg", &envmap_width, &envmap_height, &n, 0);
    if (!pixmap || 3!=n) {
        std::cerr << "Error: can not load the environment map" << std::endl;
        return;
    }
    envmap = std::vector<Vec3f>(envmap_width*envmap_height);
    for (int j = envmap_height-1; j>=0 ; j--) {
        for (int i = 0; i<envmap_width; i++) {
            envmap[i+j*envmap_width] = Vec3f(pixmap[(i+j*envmap_width)*3+0], pixmap[(i+j*envmap_width)*3+1], pixmap[(i+j*envmap_width)*3+2])*(1/255.);
        }
    }
    stbi_image_free(pixmap);

    mesh = new Model("duck.obj");
    //scene
    std::vector<Light> lights = {
        Light(Vec3f(-20, 20,  20), 1.5),
        Light(Vec3f( 30, 50, -25), 1.8),
        Light(Vec3f( 30, 20,  30), 1.7)
    };

    std::vector<Sphere> spheres = {
        Sphere(Vec3f(-3,    0,   -16), 2,      ivory),
        Sphere(Vec3f(-1.0, -1.5, -12), 2,      glass),
        Sphere(Vec3f( 1.5, -0.5, -18), 3,      red_rubber),
        Sphere(Vec3f( 7,    5,   -18), 4,      mirror)
    };

    float tan_fov_2 = tanf(fov_2);
    for (size_t y = 0; y < height; y++) {
        for (size_t x = 0; x < width; x++) {
            float wx = (x - width/2.f)*tan_fov_2/(width/2.0f);
            float wy = (y - height/2.f)*tan_fov_2*aspect/(height/2.0f);

            framebuffer[x+(height-1-y)*width] = ray_trace(
                    lights, 
                    spheres, 
                    {0.0f, 0.0f, 0.0f}, {wx, wy, -1.0f});
        }
    }
    write_ppm(framebuffer, width, height, "./out.ppm");
}


int main()
{
    render();
    return 0;
}