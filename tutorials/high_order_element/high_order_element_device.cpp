// ======================================================================== //
// Copyright 2009-2015 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#include "../common/tutorial/tutorial_device.h"
#include <math.h>
#include <iostream>

/* render function to use */
renderPixelFunc renderPixel;

/* error reporting function */
void error_handler(const RTCError code, const char* str)
{
  printf("Embree: ");
  switch (code) {
  case RTC_UNKNOWN_ERROR    : printf("RTC_UNKNOWN_ERROR"); break;
  case RTC_INVALID_ARGUMENT : printf("RTC_INVALID_ARGUMENT"); break;
  case RTC_INVALID_OPERATION: printf("RTC_INVALID_OPERATION"); break;
  case RTC_OUT_OF_MEMORY    : printf("RTC_OUT_OF_MEMORY"); break;
  case RTC_UNSUPPORTED_CPU  : printf("RTC_UNSUPPORTED_CPU"); break;
  case RTC_CANCELLED        : printf("RTC_CANCELLED"); break;
  default                   : printf("invalid error code"); break;
  }
  if (str) { 
    printf(" ("); 
    while (*str) putchar(*str++); 
    printf(")\n"); 
  }
  exit(1);
}

// ======================================================================== //
//                     User defined patch geometry                          //
// ======================================================================== //

struct Patch
{
  ALIGNED_STRUCT
  Vec3fa v[8];
  Vec3fa center;
  float radius;
  unsigned int geomID;
};

Vec3fa patchSurfaceFunc(const Patch& patch, const float u, const float v)
{
  return 0.25*(1.0 - u)*(1.0 - v)*(-u - v - 1)*patch.v[0] + 
         0.25*(1.0 + u)*(1.0 - v)*( u - v - 1)*patch.v[1] + 
         0.25*(1.0 + u)*(1.0 + v)*( u + v - 1)*patch.v[2] + 
         0.25*(1.0 - u)*(1.0 + v)*(-u + v - 1)*patch.v[3] + 
         0.5*(1 - u)*(1 - v*v)*patch.v[4] + 
         0.5*(1 - u*u)*(1 - v)*patch.v[5] + 
         0.5*(1 + u)*(1 - v*v)*patch.v[6] + 
         0.5*(1 - u*u)*(1 + v)*patch.v[7];
}

Vec3fa patchSurfaceDerivU(const Patch& patch, const float u, const float v)
{
  return (-0.25*(v - 1.0)*(u + v + 1) - 0.25*(u - 1.0)*(v - 1.0))*patch.v[0] + 
         (-0.25*(v - 1.0)*(u - v - 1) - 0.25*(u + 1.0)*(v - 1.0))*patch.v[1] + 
         ( 0.25*(v + 1.0)*(u + v - 1) + 0.25*(u + 1.0)*(v + 1.0))*patch.v[2] + 
         ( 0.25*(v + 1.0)*(u - v + 1) + 0.25*(u - 1.0)*(v + 1.0))*patch.v[3] + 
         0.5*(v*v - 1.0)*patch.v[4] + u*(v - 1.0)*patch.v[5] - 
         0.5*(v*v - 1.0)*patch.v[6] - u*(v + 1.0)*patch.v[7];
}

Vec3fa patchSurfaceDerivV(const Patch& patch, const float u, const float v)
{
  return (-0.25*(u - 1.0)*(u + v + 1) - 0.25*(u - 1.0)*(v - 1.0))*patch.v[0] + 
         (-0.25*(u + 1.0)*(u - v - 1) + 0.25*(u + 1.0)*(v - 1.0))*patch.v[1] + 
         ( 0.25*(u + 1.0)*(u + v - 1) + 0.25*(u + 1.0)*(v + 1.0))*patch.v[2] + 
         ( 0.25*(u - 1.0)*(u - v + 1) - 0.25*(u - 1.0)*(v + 1.0))*patch.v[3] + 
         0.5*(u*u - 1.0)*patch.v[5] + v*(u - 1.0)*patch.v[4] - 
         0.5*(u*u - 1.0)*patch.v[7] - v*(u + 1.0)*patch.v[6];
}

void patchBoundsFunc(const Patch* patches, size_t item, RTCBounds* bounds_o)
{
  const Patch& patch = patches[item];

  float lo_x = 1.0e300;
  float lo_y = 1.0e300;
  float lo_z = 1.0e300;

  float hi_x = -1.0e300;
  float hi_y = -1.0e300;
  float hi_z = -1.0e300;

  for (int i = 0; i < 8; ++i)
    {
      lo_x = min(lo_x, patch.v[i].x);
      lo_y = min(lo_y, patch.v[i].y);
      lo_z = min(lo_z, patch.v[i].z);
      hi_x = max(hi_x, patch.v[i].x);
      hi_y = max(hi_y, patch.v[i].y);
      hi_z = max(hi_z, patch.v[i].z);
    }

  bounds_o->lower_x = lo_x;
  bounds_o->lower_y = lo_y;
  bounds_o->lower_z = lo_z;
  bounds_o->upper_x = hi_x;
  bounds_o->upper_y = hi_y;
  bounds_o->upper_z = hi_z;
}


void patchIntersectFunc(const Patch* patches, RTCRay& ray, size_t item)
{
  const Patch& patch = patches[item];

  // return if the ray does not intersect the bounding sphere
  const Vec3fa vec = ray.org - patch.center;
  const float A = dot(ray.dir,ray.dir);
  const float B = 2.0f*dot(vec,ray.dir);
  const float C = dot(vec, vec) - sqr(patch.radius);
  const float D = B*B - 4.0f*A*C;
  if (D < 0.0f) return;

  // otherwise, iterate to the get the true hit position
  const Vec3fa dir = ray.dir / A;
  Vec3fa N1, N2;
  
  if ((fabs(dir.x) > fabs(dir.y)) && (fabs(dir.x) > fabs(dir.z))) {
      N1 = Vec3fa(dir.y, -dir.x, 0.0);
    }
  else {
    N1 = Vec3fa(0.0, dir.z, -dir.y);
  }
  N2 = cross(N1, dir);

  float d1 = - dot(N1, ray.org);
  float d2 = - dot(N2, ray.org);

  float u = 0.0;
  float v = 0.0;
  
  Vec3fa S = patchSurfaceFunc(patch, u, v);
  float fu = dot(N1, S) + d1;
  float fv = dot(N2, S) + d2;

  float err = max(fabs(fu), fabs(fv));
  float tol = 1.0e-6;
  int iterations = 0;
  int max_iter = 5;
  bool converged = false;
  Vec3fa Su, Sv;
  while ((err > tol) and (iterations < max_iter))
    {
      // compute the Jacobian
      Su = patchSurfaceDerivU(patch, u, v);
      Sv = patchSurfaceDerivV(patch, u, v);  
      float J11 = dot(N1, Su);
      float J12 = dot(N1, Sv);
      float J21 = dot(N2, Su);
      float J22 = dot(N2, Sv);
      float det = (J11*J22 - J12*J21);

      // update the u, v values
      u -= ( J22*fu - J12*fv) / det;
      v -= (-J21*fu + J11*fv) / det;

      S = patchSurfaceFunc(patch, u, v);
      fu = dot(N1, S) + d1;
      fv = dot(N2, S) + d2;
      
      err = max(fabs(fu), fabs(fv));
      iterations += 1;
    }

  if (err < tol) {
    converged = true;
  }

  float t = distance(S, ray.org);
  if ((t < ray.tnear) or (t > ray.tfar)) {
    return;
  }

  if (fabs(u) <= 1.0 and fabs(v) <= 1.0 and iterations < max_iter) {
    ray.u = u;
    ray.v = v;
    ray.tfar = t;
    ray.geomID = patch.geomID;
    ray.primID = item;
    ray.Ng = cross(Su, Sv);
  }

  return;
}

void patchOccludedFunc(const Patch* patches, RTCRay& ray, size_t item)
{
  const Patch& patch = patches[item];
  const float A = dot(ray.dir, ray.dir);
  const Vec3fa dir = ray.dir / A;
  Vec3fa N1, N2;
  
  if ((fabs(dir.x) > fabs(dir.y)) && (fabs(dir.x) > fabs(dir.z))) {
      N1 = Vec3fa(dir.y, -dir.x, 0.0);
    }
  else {
    N1 = Vec3fa(0.0, dir.z, -dir.y);
  }
  N2 = cross(N1, dir);

  float d1 = - dot(N1, ray.org);
  float d2 = - dot(N2, ray.org);

  float u = 0.0;
  float v = 0.0;
  
  Vec3fa S = patchSurfaceFunc(patch, u, v);
  float fu = dot(N1, S) + d1;
  float fv = dot(N2, S) + d2;

  float err = max(fabs(fu), fabs(fv));
  float tol = 1.0e-6;
  int iterations = 0;
  int max_iter = 5;
  bool converged = false;
  while ((err > tol) and (iterations < max_iter))
    {
      // compute the Jacobian
      Vec3fa Su = patchSurfaceDerivU(patch, u, v);
      Vec3fa Sv = patchSurfaceDerivV(patch, u, v);  
      float J11 = dot(N1, Su);
      float J12 = dot(N1, Sv);
      float J21 = dot(N2, Su);
      float J22 = dot(N2, Sv);
      float det = (J11*J22 - J12*J21);

      // update the u, v values
      u -= ( J22*fu - J12*fv) / det;
      v -= (-J21*fu + J11*fv) / det;

      S = patchSurfaceFunc(patch, u, v);
      fu = dot(N1, S) + d1;
      fv = dot(N2, S) + d2;
      
      err = max(fabs(fu), fabs(fv));
      iterations += 1;
    }

  if (err < tol) {
    converged = true;
  }

  float t = distance(S, ray.org);
  if ((t < ray.tnear) or (t > ray.tfar)) {
    return;
  }

  if (fabs(u) <= 1.0 and fabs(v) <= 1.0 and iterations < max_iter) {
    ray.geomID = 0;
  }

  return;
}

void setBoundingSphere(Patch& patch) 
{

  Vec3fa c = Vec3fa(0.0, 0.0, 0.0);
  for (int i = 0; i < 8; ++i) {
    c += patch.v[i];
  }
  c /= 8;
  patch.center = c;

  float r = 0.0;
  for (int i = 0; i < 8; ++i) {
    r = max(r, distance(patch.v[i], patch.center));
  }
  patch.radius = 1.05*r;
}

/* scene data */
RTCScene g_scene  = nullptr;

Vec3fa colors[6];

/* called by the C++ code for initialization */
extern "C" void device_init (char* cfg)
{
  /* initialize ray tracing core */
  rtcInit(cfg);

  /* set error handler */
  rtcSetErrorFunction(error_handler);

  /* create scene */
  g_scene = rtcNewScene(RTC_SCENE_STATIC,RTC_INTERSECT1);

  unsigned int geomID = rtcNewUserGeometry(g_scene, 6);
  Patch* patches = new Patch[6];

  patches[0].geomID = geomID;
  patches[0].v[0] = Vec3fa( 2.01600338,  0.30931164, -0.77566128);
  patches[0].v[1] = Vec3fa( 0.21646419,  2.16927057,  0.0948798);
  patches[0].v[2] = Vec3fa(-0.79038795,  0.33802214,  1.92612823);
  patches[0].v[3] = Vec3fa( 1.00915124, -1.52193679,  1.05558716);
  patches[0].v[4] = Vec3fa( 1.51257731, -0.60631257,  0.13996294);
  patches[0].v[5] = Vec3fa( 1.11623378,  1.23929111, -0.34039074);
  patches[0].v[6] = Vec3fa(-0.28696188,  1.25364635,  1.01050402);
  patches[0].v[7] = Vec3fa( 0.10938164, -0.59195733,  1.49085769);
  
  patches[1].geomID = geomID;
  patches[1].v[0] = Vec3fa(0.66487957, -0.40704697, -2.03809756);
  patches[1].v[1] = Vec3fa(2.01600338,  0.30931164, -0.77566128);
  patches[1].v[2] = Vec3fa(1.00915124, -1.52193679,  1.05558716);
  patches[1].v[3] = Vec3fa(-0.09433292, -2.10699819,  0.0245355);
  patches[1].v[4] = Vec3fa(0.28527332, -1.25702258, -1.00678103);
  patches[1].v[5] = Vec3fa(1.34044147, -0.04886766, -1.40687942);
  patches[1].v[6] = Vec3fa(1.51257731, -0.60631257,  0.13996294);
  patches[1].v[7] = Vec3fa(0.46585194, -1.79758637,  0.52144638);

  patches[2].geomID = geomID;
  patches[2].v[0] = Vec3fa(-1.13465961,  1.45291196, -1.16755648);
  patches[2].v[1] = Vec3fa( 0.21646419,  2.16927057,  0.0948798 );
  patches[2].v[2] = Vec3fa(-0.79038795,  0.33802214,  1.92612823);
  patches[2].v[3] = Vec3fa(-1.89387211, -0.24703926,  0.89507658);
  patches[2].v[4] = Vec3fa(-1.51426586,  0.60293635, -0.13623995);
  patches[2].v[5] = Vec3fa(-0.45909771,  1.81109126, -0.53633834);
  patches[2].v[6] = Vec3fa(-0.28696188,  1.25364635,  1.01050402);
  patches[2].v[7] = Vec3fa(-1.33368725,  0.06237256,  1.39198746);

  patches[3].geomID = geomID;
  patches[3].v[0] = Vec3fa( 0.66487957, -0.40704697, -2.03809756);
  patches[3].v[1] = Vec3fa(-1.13465961,  1.45291196, -1.16755648);
  patches[3].v[2] = Vec3fa(-1.89387211, -0.24703926,  0.89507658);
  patches[3].v[3] = Vec3fa(-0.09433292, -2.10699819,  0.0245355 );
  patches[3].v[4] = Vec3fa(  0.28527332, -1.25702258, -1.00678103);
  patches[3].v[5] = Vec3fa(-0.23489002,  0.52293249, -1.60282702);
  patches[3].v[6] = Vec3fa(-1.51426586,  0.60293635, -0.13623995);
  patches[3].v[7] = Vec3fa(-0.99410252, -1.17701872,  0.45980604);

  patches[4].geomID = geomID;
  patches[4].v[0] = Vec3fa(-0.09433292, -2.10699819,  0.0245355 );
  patches[4].v[1] = Vec3fa( 1.00915124, -1.52193679,  1.05558716);
  patches[4].v[2] = Vec3fa(-0.79038795,  0.33802214,  1.92612823);
  patches[4].v[3] = Vec3fa(-1.89387211, -0.24703926,  0.89507658);
  patches[4].v[4] = Vec3fa(-0.99410252, -1.17701872,  0.45980604);
  patches[4].v[5] = Vec3fa( 0.46585194, -1.79758637,  0.52144638);
  patches[4].v[6] = Vec3fa( 0.10938164, -0.59195733,  1.49085769);
  patches[4].v[7] = Vec3fa(-1.33368725,  0.06237256,  1.39198746);

  patches[5].geomID = geomID;
  patches[5].v[0] = Vec3fa(0.66487957, -0.40704697, -2.03809756);
  patches[5].v[1] = Vec3fa( 2.01600338,  0.30931164, -0.77566128);
  patches[5].v[2] = Vec3fa( 0.21646419,  2.16927057,  0.0948798 );
  patches[5].v[3] = Vec3fa(-1.13465961,  1.45291196, -1.16755648);
  patches[5].v[4] = Vec3fa(-0.23489002,  0.52293249, -1.60282702);
  patches[5].v[5] = Vec3fa( 1.34044147, -0.04886766, -1.40687942);
  patches[5].v[6] = Vec3fa( 1.11623378,  1.23929111, -0.34039074);
  patches[5].v[7] = Vec3fa(-0.45909771,  1.81109126, -0.53633834);

  for (int i = 0; i < 6; ++i) {
    setBoundingSphere(patches[i]);
  }

  rtcSetUserData(g_scene, geomID, patches);
  rtcSetBoundsFunction(g_scene, geomID, (RTCBoundsFunc)&patchBoundsFunc);
  rtcSetIntersectFunction(g_scene, geomID, (RTCIntersectFunc)&patchIntersectFunc);
  rtcSetOccludedFunction (g_scene, geomID, (RTCOccludedFunc )&patchOccludedFunc);
  rtcCommit(g_scene);

  /* set all colors */
  colors[0] = Vec3fa(1.0,0.0,0);
  colors[1] = Vec3fa(1.0,0.5,0);
  colors[2] = Vec3fa(1.0,1.0,0);
  colors[3] = Vec3fa(0.5,1.0,0.0);
  colors[4] = Vec3fa(1.0,0.0,1.0);
  colors[5] = Vec3fa(1.0,0.5,1.0);
  colors[6] = Vec3fa(0.5,1.0,1.0);

  /* set start render mode */
  renderPixel = renderPixelStandard;
}

/* task that renders a single screen tile */
Vec3fa renderPixelStandard(float x, float y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
{
  /* initialize ray */
  RTCRay ray;
  ray.org = p;
  ray.dir = normalize(x*vx + y*vy + vz);
  ray.tnear = 0.0f;
  ray.tfar = inf;
  ray.geomID = RTC_INVALID_GEOMETRY_ID;
  ray.primID = RTC_INVALID_GEOMETRY_ID;
  ray.instID = 4; // set default instance ID
  ray.mask = -1;
  ray.time = 0;
  
  /* intersect ray with scene */
  rtcIntersect(g_scene,ray);
  
  /* shade pixels */
  Vec3fa color = Vec3fa(0.0f);
  if (ray.geomID != RTC_INVALID_GEOMETRY_ID) 
  {
    Vec3fa diffuse = Vec3fa(0.0f);
    diffuse = colors[ray.primID];
    color = color + diffuse;

    // color the element edges white
    if ((fabs(fabs(ray.u) - 1.0) < 2.0e-2) or (fabs(fabs(ray.v) - 1.0) < 2.0e-2)) {
      color = Vec3fa(1.0, 1.0, 1.0);
    }
  }
  return color;
}

/* task that renders a single screen tile */
void renderTile(int taskIndex, int* pixels,
                     const int width,
                     const int height, 
                     const float time,
                     const Vec3fa& vx, 
                     const Vec3fa& vy, 
                     const Vec3fa& vz, 
                     const Vec3fa& p,
                     const int numTilesX, 
                     const int numTilesY)
{
  const int tileY = taskIndex / numTilesX;
  const int tileX = taskIndex - tileY * numTilesX;
  const int x0 = tileX * TILE_SIZE_X;
  const int x1 = min(x0+TILE_SIZE_X,width);
  const int y0 = tileY * TILE_SIZE_Y;
  const int y1 = min(y0+TILE_SIZE_Y,height);

  for (int y = y0; y<y1; y++) for (int x = x0; x<x1; x++)
  {
    /* calculate pixel color */
    Vec3fa color = renderPixel(x,y,vx,vy,vz,p);

    /* write color to framebuffer */
    unsigned int r = (unsigned int) (255.0f * clamp(color.x,0.0f,1.0f));
    unsigned int g = (unsigned int) (255.0f * clamp(color.y,0.0f,1.0f));
    unsigned int b = (unsigned int) (255.0f * clamp(color.z,0.0f,1.0f));
    pixels[y*width+x] = (b << 16) + (g << 8) + r;
  }
}

/* called by the C++ code to render */
extern "C" void device_render (int* pixels,
                           const int width,
                           const int height, 
                           const float time,
                           const Vec3fa& vx, 
                           const Vec3fa& vy, 
                           const Vec3fa& vz, 
                           const Vec3fa& p)
{
  /* render all pixels */
  const int numTilesX = (width +TILE_SIZE_X-1)/TILE_SIZE_X;
  const int numTilesY = (height+TILE_SIZE_Y-1)/TILE_SIZE_Y;
  launch_renderTile(numTilesX*numTilesY,pixels,width,height,time,vx,vy,vz,p,numTilesX,numTilesY); 
}

/* called by the C++ code for cleanup */
extern "C" void device_cleanup ()
{
  rtcDeleteScene (g_scene);
  rtcExit();
}
