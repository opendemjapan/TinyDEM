// TinyDEM
// Copyright (c) 2025 Roman Vetter
// vetterro@ethz.ch
// ETH Zurich

#include "tinydem.hpp"

#if defined(_MSC_VER)
// Avoid duplicate math symbol definitions when linking icx with MSVC's UCRT.
#pragma comment(linker, "/NODEFAULTLIB:libmmt.lib")
#endif

constexpr real Rmin = 0.002; // [L] minimum particle radius
constexpr real Rmax = 0.01; // [L] maximum particle radius
constexpr real Rmesh [] = {0.002,0.02,0.01}; // [L] list of mesh radii (length must match number of indices in mesh file)
constexpr Point smin = {-0.1,-0.1,0.3}; // [L] spawn box lower corner
constexpr Point smax = {0.1,0.1,0.4}; // [L] spawn box upper corner

constexpr real rho = 1e3; // [M/L^3] mass density
constexpr real E = 1e6; // [M/LT^2] Young's modulus
constexpr real nu = 0.3; // [-] Poisson's ratio
constexpr real en = 0.5; // [-] coefficient of normal restitution
constexpr real mus = 0.3; // [-] Coulomb friction coefficient for sliding
constexpr real mur = 0.3; // [-] Coulomb friction coefficient for rolling
constexpr real mut = 2 * mus / 3; // [-] Coulomb friction coefficient for twisting
constexpr real eta = 1e-2; // [M/LT] dynamic viscosity for drag

constexpr Point v0 = {0,0,-1}; // [L/T] initial linear velocity
constexpr Point w0 = {0,0,0}; // [1/T] initial angular velocity
constexpr Point g = {0,0,-9.81}; // [L/T^2] external acceleration (gravity etc.)

constexpr real PI = 3.14159265358979323846;
const real dt = Rmin * std::sqrt(rho * (1 - nu*nu) / E); // [T] time step
constexpr std::size_t Nmax = 1000; // [-] maximum number of spawning particles
constexpr std::size_t Nframes = 100; // [-] number of output frames
constexpr std::size_t Nsteps = 500; // [-] number of time steps between frames
constexpr std::size_t Nspawn = 100; // [-] number of particle spawn attempts per time step

// add linear and angular accelerations from collisions
void collision(const Point& xi, const Point& xj, const Point& vi, const Point& vj, const Point& wi, const Point& wj, const real Ri, real Rj,
               Point& ai, Point& aj, Point& alphai, Point& alphaj, std::list<Contact>& contacts, const std::size_t j, const bool is_mesh)
{
  const Point x = xi - xj; // distance vector from center to center
  const real x2 = x * x; // squared distance
  const real Rsum = Ri + Rj;
  if (x2 < Rsum * Rsum) // do the objects overlap?
  {
    const real d = std::sqrt(x2); // distance from center to center
    const Point n = 1 / d * x; // normal vector
    const real depth = Rsum - d; // depth of overlap

    // mesh elements are treated as flat walls
    if (is_mesh)
      Rj = std::numeric_limits<real>::infinity();

    // radii, masses, moments of inertia
    const real ri = Ri - depth / 2; // radius of object i to contact point
    const real rj = Rj - depth / 2; // radius of object j to contact point
    const real R  = 1 / (1 / Ri + 1 / Rj); // effective radius
    const real mi = 4 * PI / 3 * rho * Ri * Ri * Ri; // mass of object i
    const real mj = 4 * PI / 3 * rho * Rj * Rj * Rj; // mass of object j
    const real m  = 1 / (1 / mi + 1 / mj); // effective mass
    const real Ii = 2 * mi * Ri * Ri / 5; // moment of inertia of object i
    const real Ij = 2 * mj * Rj * Rj / 5; // moment of inertia of object j

    // stiffness, damping and friction coefficients
    const real a  = std::sqrt(R * depth); // contact radius
    const real kn = 4 * E / (6 * (1 - nu*nu)) * (1 + is_mesh) * a; // normal spring coefficient
    const real ks = 2 * E / ((2 - nu) * (1 + nu)) * (1 + is_mesh) * a; // shearing spring coefficient
    const real kr = kn; // rolling spring coefficient
    const real kt = ks; // twisting spring coefficient
    const real cn = std::sqrt(5 * m * kn) / std::sqrt(1 + std::pow(PI / std::log(en), 2)); // normal damping coefficient
    const real cs = std::sqrt((1 - nu) / (10 * (2 - nu))) * cn; // shearing damping coefficient
    const real cr = cn; // rolling damping coefficient
    const real ct = cs; // twisting damping coefficient

    // relative linear and angular velocity
    const Point v  = vi - vj + n.cross(ri * wi + (is_mesh ? 0 : rj) * wj); // relative velocity at contact point
    const Point vn = (v * n) * n; // normal velocity
    const Point vs = v - vn; // sliding velocity
    const Point w  = wi - wj; // relative angular velocity
    const Point vr = 1 / (1 / ri + 1 / rj) * w.cross(n); // rolling velocity
    const real  vt = R * (w * n); // twisting velocity

    // determine the contact data to use
    const auto cend = contacts.end();
    const auto it = std::find_if(contacts.begin(), cend, [j](const Contact& c){ return c.j == j; }); // locate existing contact entry
    Contact& c = it == cend ? *contacts.insert(cend, {j, {0,0,0}, {0,0,0}, 0, true}) : *it; // add a new one if none was found, otherwise use existing entry
    c.psi = true; // flag it as active so it won't be removed

    // project sliding and rolling friction springs onto current tangent plane
    const real us_before = c.us.length();
    const real ur_before = c.ur.length();
    c.us.add(-(c.us * n), n); // subtract normal component from shear displacement
    c.ur.add(-(c.ur * n), n); // subtract normal component from rolling displacement
    const real us_after = c.us.length();
    const real ur_after = c.ur.length();
    if (us_after > 0) c.us = us_before / us_after * c.us; // rescale shear displacement to the previous length
    if (ur_after > 0) c.ur = ur_before / ur_after * c.ur; // rescale rolling displacement to the previous length

    // forces
    const Point fn = kn * depth * n - cn * vn; // normal force
    Point fs = -ks * c.us - cs * vs; // shear/sliding resistance
    Point fr = -kr * c.ur - cr * vr; // rolling resistance
    real  ft = -kt * c.ut - ct * vt; // twisting resistance

    // forward Euler step for static friction springs
    c.us.add(dt, vs);
    c.ur.add(dt, vr);
    c.ut += dt * vt;

    // Coulomb's law
    const real fn_norm = fn.length();
    const real fs_norm = fs.length();
    const real fr_norm = fr.length();
    const real ft_norm = std::abs(ft);
    if (fs_norm > mus * fn_norm) // dynamic sliding?
    {
      fs = mus * fn_norm / fs_norm * fs; // clip sliding force
      c.us = -1 / ks * (fs + cs * vs); // set shear displacement to match sliding resistance
    }
    if (fr_norm > mur * fn_norm) // dynamic rolling?
    {
      fr = mur * fn_norm / fr_norm * fr; // clip rolling force
      c.ur = -1 / kr * (fr + cr * vr); // set rolling displacement to match rolling resistance
    }
    if (ft_norm > mut * fn_norm) // dynamic twisting?
    {
      ft = mut * fn_norm / ft_norm * ft; // clip twisting force
      c.ut = -1 / kt * (ft + ct * vt); // set twisting displacement to match twisting resistance
    }

    // sliding, rolling and twisting torques
    const Point taus = ri * fs.cross(n);
    const Point taur = R * n.cross(fr);
    const Point taut = ft * R * n;

    // apply the forces and torques to both objects except static mesh elements
    ai.add(1 / mi, fn + fs);
    alphai.add(1 / Ii, taus + taur + taut);
    if (!is_mesh)
    {
      aj.add(-1 / mj, fn + fs);
      alphaj.add(-1 / Ij, -(rj/ri) * taus + taur + taut);
    }
  }
}

// make a step forward in time
void Ensemble::step()
{
  // compute global bounding box
  real xmin = std::numeric_limits<real>::max(), ymin = xmin, zmin = xmin;
  real xmax = std::numeric_limits<real>::lowest(), ymax = xmax, zmax = xmax;
  #pragma omp parallel for reduction(min:xmin,ymin,zmin) reduction(max:xmax,ymax,zmax)
  for (std::size_t i = 0; i < p.size(); ++i)
  {
    if (p[i].x.x < xmin) xmin = p[i].x.x;
    if (p[i].x.x > xmax) xmax = p[i].x.x;
    if (p[i].x.y < ymin) ymin = p[i].x.y;
    if (p[i].x.y > ymax) ymax = p[i].x.y;
    if (p[i].x.z < zmin) zmin = p[i].x.z;
    if (p[i].x.z > zmax) zmax = p[i].x.z;
  }
  Box box{{xmin,ymin,zmin}, {xmax,ymax,zmax}};

  // include spawn box if it is still used
  if (Nspawn > 0 && p.size() < Nmax) // are new particles still being spawned?
    box.include(smin, smax);
  else if (p.size() == 0)
    box.xmax = box.xmin; // allows simulation to run without particles

  // build linked cell list
  Index Nc(box.xmax - box.xmin, 2 * Rmax); // number of cells in each direction
  ++Nc.x; ++Nc.y; ++Nc.z; // add one cell in all 3 directions because the cell number was rounded down
  std::vector<std::size_t> first(Nc.x * Nc.y * Nc.z, invalid); // first particle in each cell
  for (std::size_t i = 0; i < p.size(); ++i)
  {
    const auto c = Index(p[i].x - box.xmin, 2 * Rmax).global(Nc); // global cell index
    p[i].next = first[c]; // point from particle i to previously first in cell
    first[c] = i; // make particle i the first in its cell
  }

  // try to spawn new particles at random
  for (std::size_t i = 0; i < Nspawn && p.size() < Nmax; ++i)
  {
    bool overlap = false;
    const real Rnew = Rmin + rand * (Rmax - Rmin); // random radius
    const Point xnew = {smin.x + uni_dist(rng) * (smax.x - smin.x), // random position
                        smin.y + uni_dist(rng) * (smax.y - smin.y),
                        smin.z + uni_dist(rng) * (smax.z - smin.z)};
    const Index ci(xnew - box.xmin, 2 * Rmax); // cell index of the new particle

    // check distance to other particles in local Moore neighborhood
    Index cj; // cell index of the other object (particle or mesh element)
    for (cj.x = std::max(ci.x, 1l) - 1; cj.x <= ci.x + 1 && cj.x < Nc.x; ++cj.x)
      for (cj.y = std::max(ci.y, 1l) - 1; cj.y <= ci.y + 1 && cj.y < Nc.y; ++cj.y)
        for (cj.z = std::max(ci.z, 1l) - 1; cj.z <= ci.z + 1 && cj.z < Nc.z; ++cj.z)
          for (std::size_t j = first[cj.global(Nc)]; j != invalid; j = p[j].next) // loop over all particles in this cell
          {
            const Point x = xnew - p[j].x; // distance vector from center to center
            const real Rsum = Rnew + p[j].R;
            overlap = x * x < Rsum * Rsum; // do the particles overlap?
            if (overlap) goto done; // do not check further particles
          }

    // check distance to mesh elements
    for (std::size_t j = 0; j < mesh.f.size() && !overlap; ++j)
    {
      const Point x = xnew - mesh.closest(xnew, j); // distance vector from center to closest point of approach on mesh element j
      const real Rsum = Rnew + mesh.f[j].R;
      overlap = x * x < Rsum * Rsum; // does the particle clooside with the mesh element?
    }

    done:
    if (!overlap) // no overlap found with existing particles of mesh elements?
    {
      const auto c = ci.global(Nc); // global cell index
      p.push_back({Rnew, xnew, v0, {0,0,0}, w0, {0,0,0}, {1, {0,0,0}}, first[c], {}}); // create new particle
      first[c] = p.size() - 1; // make new particle the first in its cells
      rand = uni_dist(rng); // draw new random number for next particle to spawn
    }
  }

  // compute all accelerations other than from collisions
  #pragma omp parallel for
  for (std::size_t i = 0; i < p.size(); ++i)
  {
    const real c = 9 * eta / (2 * rho * p[i].R * p[i].R); // linear Stokes drag rate
    const real gam = 10 * c / 3; // angular Stokes drag rate
    p[i].a = g - c * p[i].v; // external acceleration (gravity etc.) and linear viscous damping
    p[i].alpha = -gam * p[i].w; // angular viscous damping
  }

  // particle-mesh collisions
  #pragma omp parallel for
  for (std::size_t j = 0; j < mesh.f.size(); ++j)
  {
    // check all particles in cells overlapping with bounding box of this mesh element, +/- 1 cell
    Index cmin(mesh.f[j].box.xmin - box.xmin, 2 * Rmax); // cell index of lower box corner
    Index cmax(mesh.f[j].box.xmax - box.xmin, 2 * Rmax); // cell index of upper box corner
    Index ci; // cell index of particle i
    for (ci.x = std::max(cmin.x, 1l) - 1; ci.x <= cmax.x + 1 && ci.x < Nc.x; ++ci.x)
      for (ci.y = std::max(cmin.y, 1l) - 1; ci.y <= cmax.y + 1 && ci.y < Nc.y; ++ci.y)
        for (ci.z = std::max(cmin.z, 1l) - 1; ci.z <= cmax.z + 1 && ci.z < Nc.z; ++ci.z)
          for (std::size_t i = first[ci.global(Nc)]; i != invalid; i = p[i].next) // loop over all particles in this cell
            collision(p[i].x, mesh.closest(p[i].x, j), p[i].v, {0,0,0}, p[i].w, {0,0,0}, p[i].R, mesh.f[j].R, p[i].a, p[i].a, p[i].alpha, p[i].alpha, p[i].c, -j, true);
  }

  // particle-particle collisions
  #pragma omp parallel for
  for (std::size_t i = 0; i < p.size(); ++i)
  {
    const Index ci(p[i].x - box.xmin, 2 * Rmax); // cell index of particle i
    Index cj; // cell index of the other particle j
    for (cj.x = std::max(ci.x, 1l) - 1; cj.x <= ci.x + 1 && cj.x < Nc.x; ++cj.x) // loop over local Moore neighborhood (up to 27 cells)
      for (cj.y = std::max(ci.y, 1l) - 1; cj.y <= ci.y + 1 && cj.y < Nc.y; ++cj.y)
        for (cj.z = std::max(ci.z, 1l) - 1; cj.z <= ci.z + 1 && cj.z < Nc.z; ++cj.z)
          for (std::size_t j = first[cj.global(Nc)]; j != invalid; j = p[j].next) // loop over all particles in this cell
            if (i < j) // test each particle pair only once
              collision(p[i].x, p[j].x, p[i].v, p[j].v, p[i].w, p[j].w, p[i].R, p[j].R, p[i].a, p[j].a, p[i].alpha, p[j].alpha, p[i].c, j, false);

    // remove contact data for collisions that have ended and reset the flag of all others for the next iteration
    p[i].c.remove_if([](Contact& c){ return c.psi = !c.psi; });
  }

  // time integration (semi-implicity Euler method for translation, SPIRAL for rotation)
  #pragma omp parallel for
  for (std::size_t i = 0; i < p.size(); ++i)
  {
    p[i].v.add(dt, p[i].a); // new velocity
    p[i].x.add(dt, p[i].v); // new position
    p[i].q.update(dt, p[i].w, p[i].alpha); // update quaternion before angular velocity
    p[i].w.add(dt, p[i].alpha); // new angular velocity
  }
  t += dt; // advance the time
}

int main(int argc, char** argv)
{
  std::string mfile  = (argc > 1 ? argv[1] : ""); // optional mesh file (Object File Format)
  std::string pfile  = (argc > 2 ? argv[2] : ""); // optional particle input file (CSV file in same format as output)
  std::string outdir = (argc > 3 ? argv[3] : "output"); // output directory

  Mesh mesh(mfile, Rmesh); // generate mesh
  Ensemble ensemble(pfile, mesh); // generate particle ensemble

  ensemble.output(outdir, 0); // print initial state
  for (std::size_t f = 1; f <= Nframes; ++f)
  {
    for (std::size_t s = 0; s < Nsteps; ++s)
      ensemble.step(); // make a step forward in time
    ensemble.output(outdir, f); // print a frame
  }
}
