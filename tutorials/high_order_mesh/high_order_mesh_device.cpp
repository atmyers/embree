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
  float d[8];
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

float patchInterpolationFunc(const Patch& patch, const float u, const float v)
{
  return 0.25*(1.0 - u)*(1.0 - v)*(-u - v - 1)*patch.d[0] +
    0.25*(1.0 + u)*(1.0 - v)*( u - v - 1)*patch.d[1] +
    0.25*(1.0 + u)*(1.0 + v)*( u + v - 1)*patch.d[2] +
    0.25*(1.0 - u)*(1.0 + v)*(-u + v - 1)*patch.d[3] +
    0.5*(1 - u)*(1 - v*v)*patch.d[4] +
    0.5*(1 - u*u)*(1 - v)*patch.d[5] +
    0.5*(1 + u)*(1 - v*v)*patch.d[6] +
    0.5*(1 - u*u)*(1 + v)*patch.d[7];
}

Vec3fa mapToColormap(float v,
		     float vmin,
		     float vmax)
{
  Vec3fa c;
  c.x = 1.0;
  c.y = 1.0;
  c.z = 1.0;

  float dv;

  if (v < vmin)
    v = vmin;
  if (v > vmax)
    v = vmax;
  dv = vmax - vmin;

  if (v < (vmin + 0.25 * dv)) {
    c.x = 0;
    c.y = 4 * (v - vmin) / dv;
  } else if (v < (vmin + 0.5 * dv)) {
    c.x = 0;
    c.z = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
  } else if (v < (vmin + 0.75 * dv)) {
    c.x = 4 * (v - vmin - 0.5 * dv) / dv;
    c.z = 0;
  } else {
    c.y = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
    c.z = 0;
  }

  return(c);
}

void patchIntersectFunc(const Patch* patches, RTCRay& ray, size_t item)
{
  const Patch& patch = patches[item];

  // otherwise, iterate to the get the true hit position
  const Vec3fa dir = ray.dir / dot(ray.dir, ray.dir);
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
  float tol = 1.0e-5;
  int iterations = 0;
  int max_iter = 10;
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
    ray.time = patchInterpolationFunc(patch, u, v);
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

Patch* createMesh(RTCScene scene, size_t N)
{
  unsigned int geomID = rtcNewUserGeometry(scene, N);
  Patch* patches = new Patch[N];

  int index;
  float f1, f2, f3, f4, f5, f6, f7, f8;
  float f9, f10, f11, f12, f13, f14, f15, f16;
  float f17, f18, f19, f20, f21, f22, f23, f24;
  std::ifstream meshfile ("/Users/atmyers/embree/tutorials/high_order_mesh/hex20_out.txt");
  if (meshfile.is_open())
    {
      while (meshfile >> index >> f1 >> f2 >> f3 >> f4 >> f5 >> f6 >> f7 >> f8 >>
	     f9 >> f10 >> f11 >> f12 >> f13 >> f14 >> f15 >> f16 >> f17 >> f18 >>
	     f19 >> f20 >> f21 >> f22 >> f23 >> f24)
	{

	  patches[index].geomID = geomID;
	  
	  patches[index].v[0] = Vec3fa( f1,  f2,  f3);
	  patches[index].v[1] = Vec3fa( f4,  f5,  f6);
	  patches[index].v[2] = Vec3fa( f7,  f8,  f9);
	  patches[index].v[3] = Vec3fa( f10, f11, f12);
	  patches[index].v[4] = Vec3fa( f13, f14, f15);
	  patches[index].v[5] = Vec3fa( f16, f17, f18);
	  patches[index].v[6] = Vec3fa( f19, f20, f21);
	  patches[index].v[7] = Vec3fa( f22, f23, f24);
	}
      meshfile.close();
    }

  else std::cout << "Unable to open file"; 

  std::ifstream datafile ("/Users/atmyers/embree/tutorials/high_order_mesh/hex20_temp.txt");
  if (datafile.is_open())
    {
      while (datafile >> index >> f1 >> f2 >> f3 >> f4 >> f5 >> f6 >> f7 >> f8)
	{	  
	  patches[index].d[0] = f1;
	  patches[index].d[1] = f2;
	  patches[index].d[2] = f3;
	  patches[index].d[3] = f4;
	  patches[index].d[4] = f5;
	  patches[index].d[5] = f6;
	  patches[index].d[6] = f7;
	  patches[index].d[7] = f8;
	}
      datafile.close();
    }

  else std::cout << "Unable to open file"; 

  rtcSetUserData(scene, geomID, patches);
  rtcSetBoundsFunction(scene, geomID, (RTCBoundsFunc)&patchBoundsFunc);
  rtcSetIntersectFunction(scene, geomID, (RTCIntersectFunc)&patchIntersectFunc);
  rtcSetOccludedFunction (scene, geomID, (RTCOccludedFunc )&patchOccludedFunc);
  return patches;
}

/* scene data */
RTCScene g_scene  = nullptr;

/* called by the C++ code for initialization */
extern "C" void device_init (char* cfg)
{

  /* initialize ray tracing core */
  rtcInit(cfg);

  /* set error handler */
  rtcSetErrorFunction(error_handler);

  /* create scene */
  g_scene = rtcNewScene(RTC_SCENE_STATIC,RTC_INTERSECT1);

  Patch* patches = createMesh(g_scene, 86940);

  rtcCommit(g_scene);

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
  rtcIntersect(g_scene, ray);
  
  /* shade pixels */
  Vec3fa color = Vec3fa(0.0f);
  if (ray.geomID != RTC_INVALID_GEOMETRY_ID) 
  {
    color = mapToColormap(ray.time, 200.0, 1500.0);
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
