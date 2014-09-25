// ======================================================================== //
// Copyright 2009-2014 Intel Corporation                                    //
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

#pragma once

#include "common/default.h"
#include "common/scene_subdivision.h"
#include <iomanip>

namespace embree
{
  //namespace isa
  //{
    template<typename T>
    class Array2D
    {
    public:
      Array2D () : array(NULL), size_x(0), size_y(0) {}
      ~Array2D () { delete[] array; }

      void init(size_t width, size_t height) {
        size_x = width; size_y = height;
        delete[] array; array = new T[width*height];
      }

      __forceinline size_t width() const { 
        return size_x;
      }

      __forceinline size_t height() const {
        return size_y;
      }

      __forceinline Vec3fa& operator() (size_t x, size_t y) {
        assert(x<size_x);
        assert(y<size_y);
        return array[y*size_x+x];
      }

      __forceinline const Vec3fa& operator() (size_t x, size_t y) const {
        assert(x<size_x);
        assert(y<size_y);
        return array[y*size_x+x];
      }

    private:
      T* array;
      size_t size_x, size_y;
    };

    class __aligned(64) IrregularSubdividedCatmullClarkPatch
    {
      struct Ring
      {
        Vec3fa vtx;
        Vec3fa ring[2*MAX_VALENCE]; //!< two vertices per face
        unsigned int N;             //!< number of vertices

        __forceinline Vec3fa& get(size_t i) {
          assert(i<N);
          return ring[i];
        }

        __forceinline const Vec3fa& get(size_t i) const {
          assert(i<N);
          return ring[i];
        }
        
        __forceinline const Vec3fa& first() const {
          assert(N>=4);
          return ring[4];
        }

        __forceinline const Vec3fa& last() const {
          assert(N>=2);
          return ring[N-2];
        }

        __forceinline void init(const SubdivMesh::HalfEdge* const h, const Vec3fa *const vertices)
        {
          size_t i=0;
          vtx = vertices[ h->getStartVertexIndex() ];
          SubdivMesh::HalfEdge* p = (SubdivMesh::HalfEdge*)h;
          do {
            //p = p->opposite();
            p = p->next();
            assert( i < 2*MAX_VALENCE );
            ring[i++] = vertices[ p->getStartVertexIndex() ];
            p = p->next();
            assert( i < 2*MAX_VALENCE );
            ring[i++] = vertices[ p->getStartVertexIndex() ];
            //ring[i++] = vertices[ p->prev()->getStartVertexIndex() ];
            
            /*! continue with next adjacent edge. */
            p = p->next();
            p = p->opposite();
          } while( p != h);
          N = i;
          assert( N < MAX_VALENCE );
        }

        __forceinline void subdivide (Ring& dest) const
        {
          dest.N = N;

          for (size_t i=0; i<N; i++) dest.ring[i] = Vec3fa(nan);

          /* compute face points */
          Vec3fa F(zero);
          for (size_t i=1; i<N; i+=2) {
            const Vec3fa f = 0.25f*(vtx+get(i-1)+get(i)+get((i+1)%N)); F += f; // FIXME: optimize %N
            dest.get(i) = f;
          }
          
          /* compute edge points */
          Vec3fa R(zero);
          //const Vec3fa r = 0.5f*(vtx+ring[0]); R += r;
          //const Vec3fa f = 0.5f*(dest.ring[N-1] + dest.ring[1]);
          //dest.ring[0] = 0.5f*(f+r); 
          for (size_t i=0; i<N; i+=2) {
            const Vec3fa r = 0.5f*(vtx+get(i)); R += r;
            const Vec3fa f = 0.5f*(dest.get((i-1)%N) + dest.get(i+1)); // FIXME: optimize %N
            dest.get(i) = 0.5f*(f+r); 
          }

          //PRINT("second");
          //for (size_t i=0; i<N; i++) PRINT2(i,dest.get(i));

          /* compute new point */
          const size_t valence = N/2;
          F /= (float)valence;
          R /= (float)valence; 
          dest.vtx = (F + 2.0f * R + (float)(valence-3)*vtx) / valence;
        }

        friend __forceinline std::ostream &operator<<(std::ostream &o, const Ring& h)
        {
          o << "vtx = " << h.vtx << ", ring = [";
          for (size_t i=0; i<h.N; i++) o << h.get(i) << " ";
          o << "], first = " << h.first() << ", last = " << h.last();
          return o;
        } 
    };

    public:
      Ring ring00,ring01,ring10,ring11;
      Array2D<Vec3fa> points;

      __forceinline Vec3fa& face (size_t x, size_t y) { return points(2*x,2*y); }
      __forceinline Vec3fa& vedge(size_t x, size_t y) { return points(2*x-1,2*y); }
      __forceinline Vec3fa& hedge(size_t x, size_t y) { return points(2*x,2*y-1); }
      __forceinline Vec3fa& point(size_t x, size_t y) { return points(2*x-1,2*y-1); }
      __forceinline       Vec3fa& get  (size_t x, size_t y)       { return points(x,y); }
      __forceinline const Vec3fa& get  (size_t x, size_t y) const { return points(x,y); }

      IrregularSubdividedCatmullClarkPatch () {}

      IrregularSubdividedCatmullClarkPatch (size_t width, size_t height) 
      {
        points.init(width,height);
      }

      void init (size_t width, size_t height) {
        points.init(width,height);
      }

      IrregularSubdividedCatmullClarkPatch (const SubdivMesh::HalfEdge* h, const Vec3fa *const vertices)
      {
        points.init(4,4);
        ring00.init(h,vertices); h = h->next();
        ring01.init(h,vertices); h = h->next();
        ring11.init(h,vertices); h = h->next();
        ring10.init(h,vertices); h = h->next();
        handle_corners();
      }

      void handle_corners()
      {
        size_t w = points.width();
        size_t h = points.height();
        
        points(0,0) = Vec3fa(zero);
        points(1,1) = ring00.vtx;
        points(1,0) = ring00.first();
        points(0,1) = ring00.last();

        points(0,h-1) = Vec3fa(zero);
        points(1,h-2) = ring01.vtx;
        points(1,h-1) = ring01.last();
        points(0,h-2) = ring01.first();

        points(w-1,h-1) = Vec3fa(zero);
        points(w-2,h-2) = ring11.vtx;
        points(w-2,h-1) = ring11.first();
        points(w-1,h-2) = ring11.last();

        points(w-1,0) = Vec3fa(zero);
        points(w-2,1) = ring10.vtx;
        points(w-2,0) = ring10.last();
        points(w-1,1) = ring10.first();
      }
      
      __forceinline void subdivide(IrregularSubdividedCatmullClarkPatch* out)
      {
        /* greate properly sized output patch */
        const size_t width  = points.width();
        const size_t height = points.height();
        out->init(2*width-3, 2*height-3);
        for (size_t y=0; y<2*height-3; y++)
          for (size_t x=0; x<2*width-3; x++)
            out->points(x,y) = Vec3fa(1E10);

        /* subdivide corner rings first */
        ring00.subdivide(out->ring00);
        ring01.subdivide(out->ring01);
        ring10.subdivide(out->ring10);
        ring11.subdivide(out->ring11);

        /* calculate face points */
        for (size_t y=0; y<height-1; y++) {
          for (size_t x=0; x<width-1; x++) {
            out->face(x,y) = 0.25f*(get(x+0,y+0) + get(x+0,y+1) + get(x+1,y+0) + get(x+1,y+1));
          }
        }

        /* calculate vertical edge centers */
        for (size_t y=0; y<height-1; y++) {
          for (size_t x=1; x<width-1; x++) {
            out->vedge(x,y) = 0.5f*(get(x,y)+get(x,y+1));
          }
        }

        /* calculate horizontal edge centers */
        for (size_t y=1; y<height-1; y++) {
          for (size_t x=0; x<width-1; x++) {
            out->hedge(x,y) = 0.5f*(get(x,y)+get(x+1,y));
          }
        }

        /* calculate base points */
        for (size_t y=1; y<height-1; y++) {
          for (size_t x=1; x<width-1; x++) {
            const Vec3fa F = 0.25f*(out->face(x-1,y-1)+out->face(x-1,y)+out->face(x,y-1)+out->face(x,y));
            const Vec3fa R = 0.25f*(out->hedge(x-1,y)+out->hedge(x,y)+out->vedge(x,y-1)+out->vedge(x,y));
            const Vec3fa P = get(x,y);
            out->point(x,y) = 0.25*F + 0.5*R + 0.25*P;
          }
        }

        /* calculate vertical edge points */
        for (size_t y=0; y<height-1; y++) {
          for (size_t x=1; x<width-1; x++) {
            out->vedge(x,y) = 0.5f*(out->vedge(x,y) + 0.5f*(out->face(x-1,y)+out->face(x,y)));
          }
        }

        /* calculate horizontal edge points */
        for (size_t y=1; y<height-1; y++) {
          for (size_t x=0; x<width-1; x++) {
            out->hedge(x,y) = 0.5f*(out->hedge(x,y) + 0.5f*(out->face(x,y-1)+out->face(x,y)));
          }
        }

        /* copy invalid corner points from corner rings */
        out->handle_corners();
      }

      friend __forceinline std::ostream &operator<<(std::ostream &o, const IrregularSubdividedCatmullClarkPatch& patch)
      {
        o << "ring00 = " << patch.ring00 << std::endl;
        o << "ring01 = " << patch.ring01 << std::endl;
        o << "ring10 = " << patch.ring10 << std::endl;
        o << "ring11 = " << patch.ring11 << std::endl;
        for (size_t y=0; y<patch.points.height(); y++) {
          for (size_t x=0; x<patch.points.width(); x++) {
            //o << patch.points(x,y) << " ";
            o << std::setw(10) << patch.points(x,y).z << " ";
          }
          o << std::endl;
        }
        return o;
      } 
    };
//}
}