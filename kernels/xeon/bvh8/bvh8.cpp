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

#include "bvh8.h"
#include "bvh8_statistics.h"
#include "geometry/triangle4.h"
#include "geometry/triangle8.h"
#include "common/accelinstance.h"

namespace embree
{
  DECLARE_SYMBOL(Accel::Intersector1,BVH8Triangle4Intersector1Moeller);
  DECLARE_SYMBOL(Accel::Intersector4,BVH8Triangle4Intersector4HybridMoeller);
  DECLARE_SYMBOL(Accel::Intersector4,BVH8Triangle4Intersector4HybridMoellerNoFilter);
  DECLARE_SYMBOL(Accel::Intersector8,BVH8Triangle4Intersector8ChunkMoeller);
  DECLARE_SYMBOL(Accel::Intersector8,BVH8Triangle4Intersector8HybridMoeller);
  DECLARE_SYMBOL(Accel::Intersector8,BVH8Triangle4Intersector8HybridMoellerNoFilter);

  DECLARE_SYMBOL(Accel::Intersector1,BVH8Triangle8Intersector1Moeller);
  DECLARE_SYMBOL(Accel::Intersector4,BVH8Triangle8Intersector4HybridMoeller);
  DECLARE_SYMBOL(Accel::Intersector4,BVH8Triangle8Intersector4HybridMoellerNoFilter);
  DECLARE_SYMBOL(Accel::Intersector8,BVH8Triangle8Intersector8ChunkMoeller);
  DECLARE_SYMBOL(Accel::Intersector8,BVH8Triangle8Intersector8HybridMoeller);
  DECLARE_SYMBOL(Accel::Intersector8,BVH8Triangle8Intersector8HybridMoellerNoFilter);

  DECLARE_BUILDER(void,Scene,size_t,BVH8Triangle4SceneBuilderSAH);
  DECLARE_BUILDER(void,Scene,size_t,BVH8Triangle8SceneBuilderSAH);
	DECLARE_BUILDER(void,Scene,size_t,BVH8Triangle4SceneBuilderBonsai);


  DECLARE_BUILDER(void,Scene,size_t,BVH8Triangle4SceneBuilderSpatialSAH);
  DECLARE_BUILDER(void,Scene,size_t,BVH8Triangle8SceneBuilderSpatialSAH);

  void BVH8Register ()
  {
    int features = getCPUFeatures();

    /* select builders */
    SELECT_SYMBOL_AVX(features,BVH8Triangle4SceneBuilderSAH);
    SELECT_SYMBOL_AVX(features,BVH8Triangle8SceneBuilderSAH);

	  SELECT_SYMBOL_AVX2(features,BVH8Triangle4SceneBuilderBonsai);

    SELECT_SYMBOL_AVX(features,BVH8Triangle4SceneBuilderSpatialSAH);
    SELECT_SYMBOL_AVX(features,BVH8Triangle8SceneBuilderSpatialSAH);

    /* select intersectors1 */
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle4Intersector1Moeller);
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle8Intersector1Moeller);

    /* select intersectors4 */
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle4Intersector4HybridMoeller);
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle4Intersector4HybridMoellerNoFilter);
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle8Intersector4HybridMoeller);
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle8Intersector4HybridMoellerNoFilter);

    /* select intersectors8 */
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle4Intersector8ChunkMoeller);
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle4Intersector8HybridMoeller);
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle4Intersector8HybridMoellerNoFilter);
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle8Intersector8ChunkMoeller);
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle8Intersector8HybridMoeller);
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle8Intersector8HybridMoellerNoFilter);
  }

  BVH8::BVH8 (const PrimitiveType& primTy, Scene* scene)
    : primTy(primTy), scene(scene), root(emptyNode),
	numPrimitives(0), numVertices(0) {}

  BVH8::~BVH8 () {
    for (size_t i=0; i<objects.size(); i++)
      delete objects[i];
  }

#if 0 // FIXME: remove
  void BVH8::init(size_t nodeSize, size_t numPrimitives, size_t numThreads)
  {
    /* allocate as much memory as likely needed and reserve conservative amounts of memory */
    size_t blockSize = LinearAllocatorPerThread::allocBlockSize;

    size_t numPrimBlocks = primTy.blocks(numPrimitives);
    size_t numAllocatedNodes      = min(size_t(0.6*numPrimBlocks),numPrimitives);
    size_t numAllocatedPrimitives = min(size_t(1.2*numPrimBlocks),numPrimitives);
#if defined(__X86_64__)
    size_t numReservedNodes = 2*numPrimitives;
    size_t numReservedPrimitives = 2*numPrimitives;
#else
    size_t numReservedNodes = 1.5*numAllocatedNodes;
    size_t numReservedPrimitives = 1.5*numAllocatedPrimitives;
#endif

    size_t bytesAllocated = numAllocatedNodes * nodeSize + numAllocatedPrimitives * primTy.bytes; // required also for parallel split stage in BVH4BuilderFast
    size_t bytesReserved  = numReservedNodes  * nodeSize + numReservedPrimitives  * primTy.bytes;
    if (numPrimitives) bytesReserved = (bytesReserved+blockSize-1)/blockSize*blockSize + numThreads*blockSize*2;

    root = emptyNode;
    bounds = empty;
    alloc.init(bytesAllocated,bytesReserved);
  }
#endif

  void BVH8::clear()
  {
    set(BVH8::emptyNode,empty,0);
    alloc2.clear();
  }

  void BVH8::set (NodeRef root, const BBox3fa& bounds, size_t numPrimitives)
  {
    this->root = root;
    this->bounds = bounds;
    this->numPrimitives = numPrimitives;
  }

   void BVH8::printStatistics()
   {
     std::cout << BVH8Statistics(this).str();
     std::cout << "  "; alloc2.print_statistics();
   }

  void BVH8::clearBarrier(NodeRef& node)
  {
    if (node.isBarrier())
      node.clearBarrier();
    else if (!node.isLeaf()) {
      Node* n = node.node();
      for (size_t c=0; c<N; c++)
        clearBarrier(n->child(c));
    }
  }

  void BVH8::layoutLargeNodes(size_t N)
  {
    struct NodeArea
    {
      __forceinline NodeArea() {}

      __forceinline NodeArea(NodeRef& node, const BBox3fa& bounds)
        : node(&node), A(node.isLeaf() ? float(neg_inf) : area(bounds)) {}

      __forceinline bool operator< (const NodeArea& other) const {
        return this->A < other.A;
      }

      NodeRef* node;
      float A;
    };
    std::vector<NodeArea> lst;
    lst.reserve(N);
    lst.push_back(NodeArea(root,empty));

    while (lst.size() < N)
    {
      std::pop_heap(lst.begin(), lst.end());
      NodeArea n = lst.back(); lst.pop_back();
      if (!n.node->isNode()) break;
      Node* node = n.node->node();
      for (size_t i=0; i<BVH8::N; i++) {
        if (node->child(i) == BVH8::emptyNode) continue;
        lst.push_back(NodeArea(node->child(i),node->bounds(i)));
        std::push_heap(lst.begin(), lst.end());
      }
    }

    for (size_t i=0; i<lst.size(); i++)
      lst[i].node->setBarrier();

    root = layoutLargeNodesRecursion(root);
  }

  BVH8::NodeRef BVH8::layoutLargeNodesRecursion(NodeRef& node)
  {
    if (node.isBarrier()) {
      node.clearBarrier();
      return node;
    }
    else if (node.isNode())
    {
      Node* oldnode = node.node();
      Node* newnode = (BVH8::Node*) alloc2.threadLocal2()->alloc0.malloc(sizeof(BVH8::Node)); // FIXME: optimize access to threadLocal2
      *newnode = *oldnode;
      for (size_t c=0; c<BVH8::N; c++)
        newnode->child(c) = layoutLargeNodesRecursion(oldnode->child(c));
      return encodeNode(newnode);
    }
    else return node;
  }

  Accel::Intersectors BVH8Triangle4Intersectors(BVH8* bvh)
  {
    Accel::Intersectors intersectors;
    intersectors.ptr = bvh;
    intersectors.intersector1 = BVH8Triangle4Intersector1Moeller;
    intersectors.intersector4_filter = BVH8Triangle4Intersector4HybridMoeller;
    intersectors.intersector4_nofilter = BVH8Triangle4Intersector4HybridMoellerNoFilter;
    intersectors.intersector8_filter = BVH8Triangle4Intersector8HybridMoeller;
    intersectors.intersector8_nofilter = BVH8Triangle4Intersector8HybridMoellerNoFilter;
    intersectors.intersector16 = nullptr;
    return intersectors;
  }

  Accel::Intersectors BVH8Triangle8Intersectors(BVH8* bvh)
  {
    Accel::Intersectors intersectors;
    intersectors.ptr = bvh;
    intersectors.intersector1 = BVH8Triangle8Intersector1Moeller;
    intersectors.intersector4_filter = BVH8Triangle8Intersector4HybridMoeller;
    intersectors.intersector4_nofilter = BVH8Triangle8Intersector4HybridMoellerNoFilter;
    intersectors.intersector8_filter = BVH8Triangle8Intersector8HybridMoeller;
    intersectors.intersector8_nofilter = BVH8Triangle8Intersector8HybridMoellerNoFilter;
    intersectors.intersector16 = nullptr;
    return intersectors;
  }

  Accel* BVH8::BVH8Triangle4(Scene* scene)
  {
    BVH8* accel = new BVH8(Triangle4::type,scene);
    Accel::Intersectors intersectors= BVH8Triangle4Intersectors(accel);

    Builder* builder = nullptr;
    if      (State::instance()->tri_builder == "default"     ) builder = BVH8Triangle4SceneBuilderSAH(accel,scene,0);
    else if (State::instance()->tri_builder == "binned_sah2" ) builder = BVH8Triangle4SceneBuilderSAH(accel,scene,0);
    else if (State::instance()->tri_builder == "binned_sah2_spatial" ) builder = BVH8Triangle4SceneBuilderSpatialSAH(accel,scene,0);
    else if (State::instance()->tri_builder == "binned_sah2_presplit" ) builder = BVH8Triangle4SceneBuilderSAH(accel,scene,MODE_HIGH_QUALITY);
    else THROW_RUNTIME_ERROR("unknown builder "+State::instance()->tri_builder+" for BVH8<Triangle4>");

    return new AccelInstance(accel,builder,intersectors);
  }

  Accel* BVH8::BVH8Triangle4ObjectSplit(Scene* scene)
  {
    BVH8* accel = new BVH8(Triangle4::type,scene);
    Accel::Intersectors intersectors= BVH8Triangle4Intersectors(accel);
    Builder* builder = BVH8Triangle4SceneBuilderSAH(accel,scene,0);
    return new AccelInstance(accel,builder,intersectors);
  }

	Accel* BVH8::BVH8Triangle4Bonsai(Scene* scene)
	{
		BVH8* accel = new BVH8(Triangle4::type,scene);
		Accel::Intersectors intersectors = BVH8Triangle4Intersectors(accel);
		Builder* builder = BVH8Triangle4SceneBuilderBonsai(accel,scene,0);
		return new AccelInstance(accel,builder,intersectors);
	}

  Accel* BVH8::BVH8Triangle4SpatialSplit(Scene* scene)
  {
    BVH8* accel = new BVH8(Triangle4::type,scene);
    Accel::Intersectors intersectors= BVH8Triangle4Intersectors(accel);
    Builder* builder = BVH8Triangle4SceneBuilderSpatialSAH(accel,scene,0);
    return new AccelInstance(accel,builder,intersectors);
  }

  Accel* BVH8::BVH8Triangle8(Scene* scene)
  {
    BVH8* accel = new BVH8(Triangle8::type,scene);
    Accel::Intersectors intersectors= BVH8Triangle8Intersectors(accel);

    Builder* builder = nullptr;
    if      (State::instance()->tri_builder == "default"     ) builder = BVH8Triangle8SceneBuilderSAH(accel,scene,0);
    else if (State::instance()->tri_builder == "binned_sah2" ) builder = BVH8Triangle8SceneBuilderSAH(accel,scene,0);
    else if (State::instance()->tri_builder == "binned_sah2_spatial" ) builder = BVH8Triangle8SceneBuilderSpatialSAH(accel,scene,0);
    else if (State::instance()->tri_builder == "binned_sah2_presplit" ) builder = BVH8Triangle8SceneBuilderSAH(accel,scene,MODE_HIGH_QUALITY);
    else THROW_RUNTIME_ERROR("unknown builder "+State::instance()->tri_builder+" for BVH8<Triangle8>");

    return new AccelInstance(accel,builder,intersectors);
  }

  Accel* BVH8::BVH8Triangle8ObjectSplit(Scene* scene)
  {
    BVH8* accel = new BVH8(Triangle8::type,scene);
    Accel::Intersectors intersectors= BVH8Triangle8Intersectors(accel);
    Builder* builder = BVH8Triangle8SceneBuilderSAH(accel,scene,0);
    return new AccelInstance(accel,builder,intersectors);
  }

  Accel* BVH8::BVH8Triangle8SpatialSplit(Scene* scene)
  {
    BVH8* accel = new BVH8(Triangle8::type,scene);
    Accel::Intersectors intersectors= BVH8Triangle8Intersectors(accel);
    Builder* builder = BVH8Triangle8SceneBuilderSpatialSAH(accel,scene,0);
    return new AccelInstance(accel,builder,intersectors);
  }
}
