# TinyRaytracer
My raytracer practice
![result](out.png)

technical difficulties:
* screen space to view space transformation
* sphere ray intersection
* triangle ray intersection 
* ray_trace function recursive construct
* refrection calculation(especital when ray is inside the shape)
* environment map coordinates translation
* axis-align plane rendering and checkerboard pattern generation
---
2024.4.18
1) render obj file
2) refactor the refraction calculation

2024.4.17
1) add checkerboard plane
2) add environment map

2024.3.14
1) refactor a little bit and implement basic light;
2) implement phong shading;
3) add shadow;
4) implement reflection;

2024.3.13
1) Implement sphere class with ray intersect funciton;
2) Transform screen coordiates to world coordinates in order to shoot rays;
3) Add material class and the capability to render multiple spheres;
4) And flip the final image vertically.

2024.3.12
1) setting up the project. the ppm viewer extension in vscode is very handy 