#pragma once

#include <d3d11.h>
#include "bvh.h"

#define D3D_RELEASE(Ptr) if (Ptr != nullptr) { Ptr->Release(); } 


class Timer
{
public:

   void start()
   {

   }

   void stop()
   {

   }

   float getSeconds()
   {
      return 0.0f;
   }

   float getMiliseconds()
   {
      return 0.0f;
   }

   float getMicroseconds()
   {
      return 0.0f;
   }

};

// Setups the global data for debug lines
void DebugDraw_OnInit(ID3D11Device* device, ID3D11DeviceContext* context);
void DebugDraw_OnShutdown();

class DebugDraw
{
public:
   DebugDraw(ID3D11Device* device, ID3D11DeviceContext* context);
   ~DebugDraw();

   void AddLine(math::Float4 a, math::Float4 b, unsigned int color);
   void AddAABB(math::Float4 min, math::Float4 max, unsigned int color);

   void Render(math::Matrix4 viewMatrix, math::Matrix4 projMatrix);
   
private:
   struct Vertex
   {
      float position[3];
      unsigned int color;
   };

   std::vector<Vertex> mVertices;

   ID3D11Buffer* mVertexBuffer = nullptr;
   size_t mVertexBufferSize = 0;
};
