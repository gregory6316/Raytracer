#include <iostream>
#include <cstdint>
#include <stdlib.h>
#include <limits>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <omp.h>
#include "Bitmap.h"


template<typename T>
class Vec3
{
public:
    T x, y, z;
    Vec3() : x(T(0)), y(T(0)), z(T(0)) {}
    Vec3(T xx) : x(xx), y(xx), z(xx) {}
    Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}
    Vec3& normalize()
    {
        T nor2 = length2();
        if (nor2 > 0) {
            T invNor = 1 / sqrt(nor2);
            x *= invNor, y *= invNor, z *= invNor;
        }
        return *this;
    }
    Vec3<T> operator * (const T &f) const { return Vec3<T>(x * f, y * f, z * f); }
    Vec3<T> operator / (const T &f) const { return Vec3<T>(x / f, y / f, z / f); }
    Vec3<T> operator * (const Vec3<T> &v) const { return Vec3<T>(x * v.x, y * v.y, z * v.z); }
    T dot(const Vec3<T> &v) const { return x * v.x + y * v.y + z * v.z; }
    Vec3<T> operator - (const Vec3<T> &v) const { return Vec3<T>(x - v.x, y - v.y, z - v.z); }
    Vec3<T> operator + (const Vec3<T> &v) const { return Vec3<T>(x + v.x, y + v.y, z + v.z); }
    Vec3<T>& operator += (const Vec3<T> &v) { x += v.x, y += v.y, z += v.z; return *this; }
    Vec3<T>& operator *= (const Vec3<T> &v) { x *= v.x, y *= v.y, z *= v.z; return *this; }
    Vec3 cross(const Vec3<T> &v) const
    { return Vec3<T>(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x); }
    Vec3<T> operator - () const { return Vec3<T>(-x, -y, -z); }
    T length2() const { return x * x + y * y + z * z; }
    T length() const { return sqrt(length2()); }
    friend std::ostream & operator << (std::ostream &os, const Vec3<T> &v)
    {
        os << "[" << v.x << " " << v.y << " " << v.z << "]";
        return os;
    }
};

typedef Vec3<float> Vec3f;

class Sphere
{
public:
    Vec3f center;                           
    float radius, radius2;                  
    Vec3f surfaceColor, emissionColor;     
    float transparency, reflection;         
    Sphere(
        const Vec3f &c,
        const float &r,
        const Vec3f &sc,
        const float &refl = 0,
        const float &transp = 0,
        const Vec3f &ec = 0) :
        center(c), radius(r), surfaceColor(sc), emissionColor(ec),
        transparency(transp), reflection(refl)
    {  }
    bool intersect(const Vec3f &rayorig, const Vec3f &raydir, float &t0, float &t1) const
    {
        Vec3f l = center - rayorig;
        float tca = l.dot(raydir);
        if (tca < 0) return false;
        float d2 = l.dot(l) - tca * tca;
        if (d2 > radius*radius) return false;
        float thc = sqrt(radius*radius - d2);
        t0 = tca - thc;
        t1 = tca + thc;
        
        return true;
    }
};
class Triangle
{
public:
    Vec3f vec0, vec1, vec2;
    Vec3f surfaceColor, emissionColor;
    float transparency, reflection;
    Triangle(
        const Vec3f &v0,
        const Vec3f &v1,
        const Vec3f &v2,
        const Vec3f &sc) :
        vec0(v0), vec1(v1), vec2(v2), surfaceColor(sc)
    {  }
    bool rayTriangleIntersect(
        const Vec3f &rayorig, const Vec3f &raydir,
        float &t, float &u, float &v) const
    {
        Vec3f v0v1 = vec1 - vec0;
        Vec3f v0v2 = vec2 - vec0;
        Vec3f pvec = raydir.cross(v0v2);
        float det = v0v1.dot(pvec); 
        if (fabs(det) < 1e-8) return false;

        float invDet = 1 / det;

        Vec3f tvec = rayorig - vec0;
        u = tvec.dot(pvec) * invDet;
        if (u < 0 || u > 1) return false;

        Vec3f qvec = tvec.cross(v0v1);
        v = raydir.dot(qvec) * invDet;
        if (v < 0 || u + v > 1) return false;
    
        t = v0v2.dot(qvec) * invDet;
    
        return true;

    }
};





Vec3f trace(
    const Vec3f &rayorig,
    const Vec3f &raydir,
    const std::vector<Sphere> &spheres,
    const std::vector<Triangle> &triangles,
    const int &depth)
{
    float tnear = 1e8;
    const Sphere* sphere = NULL;
    const Triangle* triangle = NULL;
    for (unsigned i = 0; i < spheres.size(); ++i) {
        float t0 = 1e8, t1 = 1e8;
        if (spheres[i].intersect(rayorig, raydir, t0, t1)) {
            if (t0 < 0) t0 = t1;
            if (t0 < tnear) {
                tnear = t0;
                sphere = &spheres[i];
            }
        }
    }
    for (unsigned i = 0; i < triangles.size(); ++i) {
        float t0 = 1e8, t1 = 1e8, t2 = 1e8;
        if (triangles[i].rayTriangleIntersect(rayorig, raydir, t0, t1, t2)) {
            if (t0 < 0) t0 = t1;
            if (t0 < tnear) {
                tnear = t0;
                triangle = &triangles[i];
            }
        }
    }
    int shadow = 0;
    float spheres_dist = std::numeric_limits<float>::max();
    float bias = 1e-4;
    Vec3f surfaceColor = 0;
    if (fabs(raydir.y)>1e-3 && (!sphere || sphere == &spheres[spheres.size()-1]) && !triangle)  {
        float d = -(rayorig.y+4)/raydir.y; 
        Vec3f pt = rayorig + raydir*d;
        if (d>0 && fabs(pt.x)<1000 && pt.z<100 && pt.z>-100 && d<spheres_dist) {
            surfaceColor = (int(.5*pt.x+1000) + int(.5*pt.z)) & 1 ? Vec3f(0.5, 0.3, 0.0) : Vec3f(.9, .9, .9);
            shadow = 1;
        }
     }
    if ((!sphere || sphere == &spheres[spheres.size()-1])&&(shadow == 0) && !triangle){
         return surfaceColor;
    };
    if (shadow == 0 && !triangle){
        Vec3f phit = rayorig + raydir * tnear; 
        Vec3f nhit = phit - sphere->center; 
        nhit.normalize(); 
        bool inside = false;
        if (raydir.dot(nhit) > 0) nhit = -nhit, inside = true;
            if ((sphere->transparency > 0 || sphere->reflection > 0) && depth < 5) {
        
                Vec3f refldir = raydir - nhit * 2 * raydir.dot(nhit);
                refldir.normalize();
                Vec3f reflection = trace(phit + nhit * bias, refldir, spheres, triangles, depth + 1);
                Vec3f refraction = 0;
                unsigned flag = 1;
        
                if (sphere->transparency) {
                    float ior = 1.1, eta = (inside) ? ior : 1 / ior; 
                    float cosi = -nhit.dot(raydir);
                    flag = 0;
                    float k = 1 - eta * eta * (1 - cosi * cosi);
                    Vec3f refrdir = (raydir * eta + nhit * (eta *  cosi - sqrt(k)));
                    refrdir.normalize();
                    refraction = trace(phit - nhit * bias, refrdir, spheres, triangles, depth+1 );
                }
                surfaceColor = (
                reflection  * flag +
                refraction  * sphere->transparency) ;
            }
            else {
                for (unsigned i = 0; i < spheres.size(); ++i) {
                    if (spheres[i].emissionColor.x > 0) {
                        Vec3f transmission = 1;
                        Vec3f lightDirection = spheres[i].center - phit;
                        lightDirection.normalize();
                        for (unsigned j = 0; j < spheres.size(); ++j) {
                            if (i != j) {
                                float t0, t1;
                                if (spheres[j].intersect(phit + nhit * bias, lightDirection, t0, t1)) {
                                    transmission = 0;
                                    break;
                                }
                            }
                        }
                        for (unsigned j = 0; j < triangles.size(); ++j) {
                                    float t0, t1, t2;
                                    if (triangles[j].rayTriangleIntersect(phit + nhit * bias, lightDirection, t0, t1, t2)) {
                                        transmission = 0;
                                        if (sphere->center.z > triangles[j].vec1.z){
                                            transmission = 1;
                                        }
                                        break;
                                    }
                        }
                        surfaceColor += sphere->surfaceColor * transmission *
                        std::max(float(0), nhit.dot(lightDirection)) ;
                    }
                }
            }
    
            return surfaceColor + sphere->emissionColor;
        } 
        Vec3f transmission = 1;
        if (triangle){
            Vec3f phit = rayorig + raydir * tnear; 
            Vec3f nhit = phit - (triangle->vec0+triangle->vec1+triangle->vec2)/3;
            nhit.normalize(); 
            for (unsigned i = 0; i < spheres.size(); ++i) {
                        if (spheres[i].emissionColor.x > 0) {
                            Vec3f lightDirection = spheres[i].center - phit;
                            lightDirection.normalize();
                            surfaceColor = spheres[i].emissionColor;
                            for (unsigned j = 0; j < spheres.size(); ++j) {
                                if (i != j) {
                                    float t0, t1;
                                    if (spheres[j].intersect(phit + nhit * bias, lightDirection, t0, t1)) {
                                        transmission = 0;
                                        break;
                                    }
                                }
                            }
                               for (unsigned j = 0; j < triangles.size(); ++j) {
                                    float t0, t1, t2;
                                    if (triangles[j].rayTriangleIntersect(phit + nhit * bias, lightDirection, t0, t1, t2) && triangle != &triangles[j]) {
                                        transmission = 0;
                                        break;
                                    }
                                }    
                        }
            }
                 return triangle->surfaceColor * transmission;
        }  
        if (sphere == &spheres[spheres.size()-1]){
            Vec3f phit = rayorig + raydir * tnear; 
            Vec3f nhit = phit - sphere->center; 
            nhit.normalize();
            for (unsigned i = 0; i < spheres.size(); ++i) {
                        if (spheres[i].emissionColor.x > 0) {
                            Vec3f lightDirection = spheres[i].center - phit;
                            lightDirection.normalize();
                            for (unsigned j = 0; j < spheres.size(); ++j) {
                                if (i != j) {
                                    float t0, t1;
                                    if (spheres[j].intersect(phit + nhit * bias, lightDirection, t0, t1)) {
                                        transmission = 0;
                                        break;
                                    }
                                
                            }
                            
                                for (unsigned j = 0; j < triangles.size(); ++j) {
                                    float t0, t1, t2;
                                    if (triangles[j].rayTriangleIntersect(phit + nhit * bias, lightDirection, t0, t1, t2)) {
                                        transmission = 0;
                                        break;
                                    }
                                }    
                            

                            }
                    }
            }
        }
        return surfaceColor * transmission;
    }



void render(const std::vector<Sphere> &spheres, const std::vector<Triangle> &triangles, int argc, const char** argv)
{
    std::unordered_map<std::string, std::string> cmdLineParams;
    int threads = 2;
    for(int i=0; i<argc; i++)
    {
      std::string key(argv[i]);

      if(key.size() > 0 && key[0]=='-')
      {
        if(i != argc-1) 
        {
          cmdLineParams[key] = argv[i+1];
          i++;
        }
        else
          cmdLineParams[key] = "";
      }
    }
    std::string outFilePath = "319_Nanyan.bmp";
  if(cmdLineParams.find("-out") != cmdLineParams.end())
    outFilePath = cmdLineParams["-out"];

  int sceneId = 0;
  if(cmdLineParams.find("-scene") != cmdLineParams.end())
    sceneId = atoi(cmdLineParams["-scene"].c_str());
  if(cmdLineParams.find("-threads") != cmdLineParams.end())
    threads = stoi(cmdLineParams["-threads"]);

  
  
  std::vector<uint32_t> images(1024*768); 

  

  std::cout << "end. " << std::endl;
    unsigned width = 1024, height = 768;
    Vec3f *image = new Vec3f[width * height], *pixel = image;
    float invWidth = 1 / float(width), invHeight = 1 / float(height);
    float fov = 30, aspectratio = width / float(height);
    float angle = tan(3.141592653589793 * 0.5 * fov / 180.);
    unsigned x;
    unsigned y;
    Vec3f *pixels = pixel; 
    omp_set_num_threads(threads);
    #pragma omp parallel for private(y, x)  
    for ( y = 0; y < height; ++y) {
        for ( x = 0; x < width; ++x) {
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
            Vec3f raydir(xx, yy, -1);
            raydir.normalize();
            pixel = pixels + x+(y*width);
            *pixel = trace(Vec3f(0), raydir, spheres, triangles, 0);
        }
    }
    
    for (unsigned i = 1; i <= width * height; ++i) {
      if (i / width != height){
        images[(height-i/width)*width-(width - i % width)] = (int(image[i-1].z * 255))*65536 + (int(image[i-1].y*255))*256 + (int(image[i-1].x*255));
      }
      else{
        images[width] = (int(image[width].z * 255))*65536 + (int(image[width].y*255))*256 + (int(image[width].x*255));
      }
    }
    delete [] image;
  if(sceneId == 1 || sceneId == 0)
    SaveBMP(outFilePath.c_str(), images.data(), 1024, 768);
  

}

int main(int argc, const char** argv)
{

  std::vector<Sphere> spheres;
  std::vector<Triangle> triangles;
    // координаты вершин, цвет
    triangles.push_back(Triangle(Vec3f(-9, -4, -30), Vec3f(-5, -4, -30), Vec3f(-7, -1, -30), Vec3f(0.8,0.8,0.8)));
    // позиция, радиус, цвет, отражение, прозрачность, свет
    spheres.push_back(Sphere(Vec3f( 3.0,      -3, -28),     1, Vec3f(0.80, 0.80, 0.80), 1.0, 0.0));//зеркальный
    spheres.push_back(Sphere(Vec3f( -7.0,      -1, -25),     0.5, Vec3f(1.0, 0.80, 0.80), 0.0, 0.0));//матовый
    spheres.push_back(Sphere(Vec3f( 4, -1.3, -18),     1, Vec3f(0.90, 0.76, 0.46), 0, 0.0));//матовый
    spheres.push_back(Sphere(Vec3f( -1.5, -1.5, -12),     1, Vec3f(1.0, 1.0, 1.0), 1.0, 1.0));//стеклянный
    //свет
    spheres.push_back(Sphere(Vec3f( 0.0,     20, 30),     1, Vec3f(0.00, 0.00, 0.00), 0, 0.0, Vec3f(3)));
    spheres.push_back(Sphere(Vec3f( 0.0,     40, -100),     1, Vec3f(0.00, 0.00, 0.00), 0, 0.0, Vec3f(3)));
    //тень
    spheres.push_back(Sphere(Vec3f( 0.0, -10004, -20), 10000, Vec3f(0.20, 0.20, 0.20), 0, 0.0));
    render(spheres, triangles, argc, argv);
    
    return 0;

  

}
