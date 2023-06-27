#pragma once

#include <d3d11.h>
#include <DirectXMath.h>
#include <chrono>

#include "bvh.h"
#include "imgui.h"

#define D3D_RELEASE(Ptr) if (Ptr != nullptr) { Ptr->Release(); } 

// Setups the global data for debug lines
void DebugDraw_OnInit(ID3D11Device* device, ID3D11DeviceContext* context);
void DebugDraw_OnShutdown();

void ImGui_DebugDrawLine_SLOW(const DirectX::XMMATRIX& projViewMatrix, math::Float4 a, math::Float4 b, ImU32 col);
void ImGui_DebugDrawAABB_SLOW(const DirectX::XMMATRIX& projViewMatrix, math::Float4 min, math::Float4 max, ImU32 col);

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

   ID3D11Buffer* mConstantBuffer = nullptr;
};
