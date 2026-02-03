// TinyDEM
// Copyright (c) 2025 Roman Vetter
// vetterro@ethz.ch
// ETH Zurich

#ifndef TINYDEM_HPP
#define TINYDEM_HPP

// Allow building with older Clang-based compilers against newer MSVC STL.
// This avoids STL1000 when using icx 2023.x with VS2022 17.10+ headers.
#ifndef _ALLOW_COMPILER_AND_STL_VERSION_MISMATCH
#define _ALLOW_COMPILER_AND_STL_VERSION_MISMATCH
#endif

// Work around icx 2023.x handling of "static_assert(false)" in discarded
// if constexpr branches inside the MSVC STL (Clang 16 bug).
// We temporarily replace static_assert with a no-op declaration for STL headers.
#ifndef _ALLOW_KEYWORD_MACROS
#define _ALLOW_KEYWORD_MACROS
#endif
#ifndef TINYDEM_STL_STATIC_ASSERT_HACK
#define TINYDEM_STL_STATIC_ASSERT_HACK
#define TINYDEM_STATIC_ASSERT_GLUE(a, b) a##b
#define TINYDEM_STATIC_ASSERT_JOIN(a, b) TINYDEM_STATIC_ASSERT_GLUE(a, b)
#define static_assert(...) typedef int TINYDEM_STATIC_ASSERT_JOIN(_stl_static_assert_, __COUNTER__)[1]
#endif

#include <vector>
#include <list>
#include <cmath>
#include <random>
#include <limits>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <filesystem>

#ifdef static_assert
#undef static_assert
#endif
#ifdef TINYDEM_STL_STATIC_ASSERT_HACK
#undef TINYDEM_STATIC_ASSERT_JOIN
#undef TINYDEM_STATIC_ASSERT_GLUE
#undef TINYDEM_STL_STATIC_ASSERT_HACK
#endif

typedef double real; // floating point precision

// minimal class providing the required vector algebra in 3D
struct Point
{
  real x, y, z;
  Point operator+(const Point& p) const { return {x + p.x, y + p.y, z + p.z}; } // sum
  Point operator-(const Point& p) const { return {x - p.x, y - p.y, z - p.z}; } // difference
  real operator*(const Point& p) const { return x * p.x + y * p.y + z * p.z; } // dot product
  Point cross(const Point& p) const { return {y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x}; } // cross product
  real length() const { return std::sqrt(x * x + y * y + z * z); } // norm
  void add(const real a, const Point& p) { // add a * p to this point
    #pragma omp atomic
    x += a * p.x;
    #pragma omp atomic
    y += a * p.y;
    #pragma omp atomic
    z += a * p.z;
  }
};
Point operator*(const real a, const Point& p) { return {a * p.x, a * p.y, a * p.z}; } // scalar multiplication

// minimal bounding box struct
struct Box
{
  Point xmin, xmax;
  void include(const Point& p0, const Point& p1) // assuming p1 >= p0
  {
    if (p0.x < xmin.x) xmin.x = p0.x;
    if (p0.y < xmin.y) xmin.y = p0.y;
    if (p0.z < xmin.z) xmin.z = p0.z;
    if (p1.x > xmax.x) xmax.x = p1.x;
    if (p1.y > xmax.y) xmax.y = p1.y;
    if (p1.z > xmax.z) xmax.z = p1.z;
  }
};

// minimal class providing the required quaternion algebra
struct Quaternion
{
  real s; // scalar part
  Point v; // vectorial part
  void mult(const Quaternion& q) // Hamilton product
  {
    const real snew = s * q.s - v * q.v;
    v = s * q.v + q.s * v + v.cross(q.v);
    s = snew;
  }
  void update(const real dt, const Point& w, const Point& alpha) // SPIRAL quaternion update
  {
    const real w_norm = w.length();
    if (w_norm > 0)
    {
      const real phi = dt / 2 * w_norm;
      mult({std::cos(phi), std::sin(phi) / w_norm * w});
    }
    const real alpha_norm = alpha.length();
    if (alpha_norm > 0)
    {
      const real theta = dt * dt / 4 * alpha_norm;
      mult({std::cos(theta), std::sin(theta) / alpha_norm * alpha});
    }
  }
};

// data struct for contacts
struct Contact
{
  std::size_t j; // index of the other object (particle if >0, mesh if overflow to <=0)
  Point us, ur; // sliding and rolling friction spring displacements
  real ut; // twisting spring displacement
  bool psi; // flag used to remove contacts that have ended
};

// minimal particle struct
struct Particle
{
  real R; // radius
  Point x, v, a, w, alpha; // position, velocity, acceleration, angular velocity, angular acceleration
  Quaternion q; // orientation
  std::size_t next; // next particle in same cell
  std::list<Contact> c; // list of contact data for collisions (only particles with greater index than this)
};

// minimal spatial cell index class
struct Index
{
  long int x, y, z; // signed indices in x, y, z directions
  Index(const Point& p = {0,0,0}, real d = 1) : x(p.x / d), y(p.y / d), z(p.z / d) {}
  std::size_t global(const Index& N) const { return x * N.y * N.z + y * N.z + z; } // global index in array of total size N
};

// minimal mesh face struct
struct Face
{
  std::vector<std::size_t> v; // vertices
  real R; // semi-thickness
  Box box; // bounding box
};

// returns the closest point on the line segment (a,b) to point p
Point edge_point(const Point& p, const Point& a, const Point& b)
{
  const Point ab = b - a;
  const Point ap = p - a;
  const real s = ab * ap / (ab * ab);
  if (s <= 0) return a; // closest point is vertex a
  if (s >= 1) return b; // closest point is vertex b
  return a + s * ab; // closest point is between a and b
}

// returns the closest point on the triangle (a,b,c) to point p
Point triangle_point(const Point& p, const Point& a, const Point& b, const Point& c)
{
  const Point ab = b - a;
  const Point ac = c - a;
  const Point ap = p - a;
  const real d1 = ab * ap;
  const real d2 = ac * ap;
  if (d1 <= 0 && d2 <= 0)
    return a; // closest point is vertex a

  const Point bp = p - b;
  const real d3 = ab * bp;
  const real d4 = ac * bp;
  if (d3 >= 0 && d4 <= d3)
    return b; // closest point is vertex b

  const real vc = d1 * d4 - d3 * d2;
  if (vc <= 0 && d1 >= 0 && d3 <= 0)
    return a + (d1 / (d1 - d3)) * ab; // closest point is on ab edge

  const Point cp = p - c;
  const real d5 = ab * cp;
  const real d6 = ac * cp;
  if (d6 >= 0 && d5 <= d6)
    return c; // closest point is vertex c

  const real vb = d5 * d2 - d1 * d6;
  if (vb <= 0 && d2 >= 0 && d6 <= 0)
    return a + (d2 / (d2 - d6)) * ac; // closest point is on ac edge

  const real va = d3 * d6 - d5 * d4;
  const real d43 = d4 - d3;
  const real d56 = d5 - d6;
  if (va <= 0 && d43 >= 0 && d56 >= 0)
    return b + (d43 / (d43 + d56)) * (c - b); // closest point is on bc edge

  const real invv = 1 / (va + vb + vc);
  return a + (vb * invv) * ab + (vc * invv) * ac; // closest point is in interior
}

// returns the closest point on the rectangle (a,b,c,d=a+ab+ac) to point p
Point rectangle_point(const Point& p, const Point& a, const Point& b, const Point& c)
{
  const Point ab = b - a;
  const Point ac = c - a;
  const Point ap = p - a;
  Point q = a;

  // clamp in the ab direction
  real s = ab * ap / (ab * ab);
  if (s > 0)
    q.add(s < 1 ? s : 1, ab);

  // clamp in the ac direction
  s = ac * ap / (ac * ac);
  if (s > 0)
    q.add(s < 1 ? s : 1, ac);

  return q;
}

// minimal class for polygonal meshes in the Open File Format
struct Mesh
{
  std::vector<Point> p; // list of points
  std::vector<Face> f; // list of faces

  // read a mesh from an OFF file
  Mesh(const std::string filename, const real Rmesh [])
  {
    if (filename.empty()) return;

    std::ifstream file(filename.c_str());
    readline(file); // skip header line "OFF"
    std::size_t Np, Nf, Ne, index;
    readline(file) >> Np >> Nf >> Ne; // Ne is unused

    // read points
    p.resize(Np);
    for (auto& point : p)
      readline(file) >> point.x >> point.y >> point.z;

    // read faces
    f.resize(Nf);
    for (std::size_t i = 0; i < Nf; ++i)
    {
      std::istringstream iss = readline(file);
      iss >> Np;
      f[i].v.resize(Np);
      for (std::size_t j = 0; j < Np; ++j)
        iss >> f[i].v[j];

      // read index for Rmesh array, if one is given (defaults to 0)
      index = 0;
      iss >> index;
      f[i].R = Rmesh[index];

      // precompute bounding box
      f[i].box = {p[f[i].v[0]], p[f[i].v[0]]};
      for (auto& v : f[i].v)
        f[i].box.include(p[v], p[v]);
      f[i].box.xmin.add(-f[i].R, {1,1,1});
      f[i].box.xmax.add( f[i].R, {1,1,1});
    }
  }

  // read a line from a file, skipping empty and comment-only lines
  std::istringstream readline(std::ifstream& file)
  {
    std::string line;
    do std::getline(file, line);
    while (line.empty() || line[0] == '#');
    return std::istringstream(line);
  }
  
  // returns the point on face j that is closest to x
  Point closest(const Point& x, std::size_t j) const
  {
    if (f[j].v.size() == 1) // point
      return p[f[j].v[0]];
    else if (f[j].v.size() == 2) // edge
      return edge_point(x, p[f[j].v[0]], p[f[j].v[1]]);
    else if (f[j].v.size() == 3) // triangle
      return triangle_point(x, p[f[j].v[0]], p[f[j].v[1]], p[f[j].v[2]]);
    else if (f[j].v.size() == 4) // rectangle
      return rectangle_point(x, p[f[j].v[0]], p[f[j].v[1]], p[f[j].v[3]]);
    else return x; // mesh polygons with more vertices are not implemented
  }
};

// ensemble of particles
struct Ensemble
{
  std::vector<Particle> p; // list of particles
  Mesh& mesh;
  real t = 0; // time
  static constexpr std::size_t invalid = std::numeric_limits<std::size_t>::max(); // invalid particle index
  std::mt19937 rng; // random number generator
  std::uniform_real_distribution<real> uni_dist; // standard uniform distribution U[0,1)
  real rand = uni_dist(rng); // random number for the next particle

  // read a particle ensemble from a CSV file
  Ensemble(const std::string filename, Mesh& m, const std::size_t N = 0) : mesh(m)
  {
    if (filename.empty()) return;
    p.reserve(N); // optionally preallocate memory
    std::string line;
    std::ifstream file(filename.c_str());
    std::getline(file, line); // skip header line
    real x, y, z, R, vx, vy, vz, q0, q1, q2, q3, wx, wy, wz;
    while (std::getline(file, line))
    {
      std::replace(line.begin(), line.end(), ',', ' ');
      std::istringstream iss(line);
      iss >> x >> y >> z >> R >> vx >> vy >> vz >> q0 >> q1 >> q2 >> q3 >> wx >> wy >> wz;
      p.push_back({R, {x,y,z}, {vx,vy,vz}, {0,0,0}, {wx,wy,wz}, {0,0,0}, {q0, {q1,q2,q3}}, invalid, {}});
    }
  }

  void step(); // defined externally

  // print & write output
  void output(const std::string outdir, const std::size_t f) const
  {
    // print frame number, simulated time, and number of particles to stdout
    std::cout << "frame " << f << ": t=" << t << ", " << p.size() << " particles" << std::endl;

    // create the output directory if it doesn't exist yet
    std::error_code ec;
    std::filesystem::create_directories(outdir, ec);

    // write current particle state to a CSV file
    char filename [20];
    snprintf(filename, 20, "particles%06zu.csv", f);
    std::ofstream file(outdir + "/" + filename);
    file << std::setprecision(std::numeric_limits<real>::max_digits10);
    file << "x,y,z,R,vx,vy,vz,q0,q1,q2,q3,wx,wy,wz\n";
    for (auto& pi : p)
      file << pi.x.x << "," << pi.x.y << "," << pi.x.z << "," << pi.R << ","
           << pi.v.x << "," << pi.v.y << "," << pi.v.z << ","
           << pi.q.s << "," << pi.q.v.x << "," << pi.q.v.y << "," << pi.q.v.z << ","
           << pi.w.x << "," << pi.w.y << "," << pi.w.z << "\n";
  }
};

#endif // TINYDEM_HPP
