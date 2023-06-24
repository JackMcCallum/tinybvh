#pragma once
#include "utils.h"
#include "imgui.h"

#include <vector>
#include <string>
#include <d3dcompiler.h> // shader compiler

ID3DBlob* CompileShaderBlob(const char* code, const char* entry, const char* target)
{
   ID3DBlob* codeBlob = nullptr;
   ID3DBlob* errorBlob = nullptr;
   HRESULT hr = D3DCompile(code, strlen(code), nullptr, nullptr, nullptr, entry, target, 0, 0, &codeBlob, &errorBlob);

   if (errorBlob)
   {
      std::string errorString;
      errorString.resize(errorBlob->GetBufferSize());
      memcpy(&errorString[0], errorBlob->GetBufferPointer(), errorBlob->GetBufferSize());
      printf("Errors/Warnings while compiling shader\n%s", errorString.c_str());
      errorBlob->Release();
   }

   assert(SUCCEEDED(hr));

   return codeBlob;
}

ID3D11VertexShader* CompileVertexShader(ID3D11Device* device, const char* code, const char* entry)
{
   ID3D11VertexShader* shader = nullptr;

   ID3DBlob* blob = CompileShaderBlob(code, entry, "vs_5_0");
   assert(blob);

   if (blob)
   {
      HRESULT hr = device->CreateVertexShader(blob->GetBufferPointer(), blob->GetBufferSize(), nullptr, &shader);
      assert(SUCCEEDED(hr));

      blob->Release();
   }

   return shader;
}

ID3D11PixelShader* CompilePixelShader(ID3D11Device* device, const char* code, const char* entry)
{
   ID3D11PixelShader* shader = nullptr;

   ID3DBlob* blob = CompileShaderBlob(code, entry, "ps_5_0");
   assert(blob);

   if (blob)
   {
      HRESULT hr = device->CreatePixelShader(blob->GetBufferPointer(), blob->GetBufferSize(), nullptr, &shader);
      assert(SUCCEEDED(hr));
   }

   return shader;
}

ID3D11Buffer* CreateConstantBuffer(ID3D11Device* device, UINT size)
{
   D3D11_BUFFER_DESC desc = {};
   desc.ByteWidth = size;
   desc.ByteWidth = (desc.ByteWidth + 255) & ~255;
   desc.Usage = D3D11_USAGE_DYNAMIC;
   desc.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
   desc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;

   ID3D11Buffer* buffer = nullptr;
   HRESULT hr = device->CreateBuffer(&desc, nullptr, &buffer);
   assert(SUCCEEDED(hr));

   return buffer;
}

ID3D11Buffer* CreateVertexBuffer(ID3D11Device* device, UINT size)
{
   D3D11_BUFFER_DESC desc = {};
   desc.ByteWidth = size;
   desc.ByteWidth = (desc.ByteWidth + 255) & ~255;
   desc.Usage = D3D11_USAGE_DYNAMIC;
   desc.BindFlags = D3D11_BIND_VERTEX_BUFFER;
   desc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;

   ID3D11Buffer* buffer = nullptr;
   HRESULT hr = device->CreateBuffer(&desc, nullptr, &buffer);
   assert(SUCCEEDED(hr));

   return buffer;
}

void UpdateBuffer(ID3D11DeviceContext* context, ID3D11Buffer* buffer, void* data, UINT size)
{
   D3D11_MAPPED_SUBRESOURCE mapped = {};
   HRESULT hr = context->Map(buffer, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped);
   assert(SUCCEEDED(hr));




}



const char* debugLinesShader = R"(

struct VSInput
{
   float3 position : POSITION;
   float4 color : COLOR;
};

struct VSOutput
{
   float4 position : SV_POSITION;
   float4 color : TEXCOORD0;
};

cbuffer Camera : register(b0)
{
   float4x4 viewProjMatrix;
};

VSOutput VSMain(VSInput IN)
{
   VSOutput OUT;
   OUT.position = mul(float4(IN.position, 1), viewProjMatrix);
   OUT.color = IN.color;
   return OUT;
}

float4 PSMain(VSOutput IN) : SV_Target0
{
   return float4(IN.color.rgb, IN.color.a);
}

)";

ID3D11VertexShader* mDebugLinesVS = nullptr;
ID3D11PixelShader* mDebugLinesPS = nullptr;

void DebugDraw_OnInit(ID3D11Device* device, ID3D11DeviceContext* context)
{
   mDebugLinesVS = CompileVertexShader(device, debugLinesShader, "VSMain");
   mDebugLinesPS = CompilePixelShader(device, debugLinesShader, "PSMain");
}

void DebugDraw_OnShutdown()
{
   D3D_RELEASE(mDebugLinesVS);
   D3D_RELEASE(mDebugLinesPS);
}

DebugDraw::DebugDraw(ID3D11Device* device, ID3D11DeviceContext* context)
{
   mConstantBuffer = CreateBuf
}

DebugDraw::~DebugDraw()
{
}

void DebugDraw::AddLine(math::Float4 a, math::Float4 b, unsigned int color)
{
   auto ac = a.SplitComponents();
   auto bc = b.SplitComponents();

   Vertex vertex;
   vertex.color = color;

   vertex.position[0] = ac.x;
   vertex.position[1] = ac.y;
   vertex.position[2] = ac.z;
   mVertices.push_back(vertex);

   vertex.position[0] = bc.x;
   vertex.position[1] = bc.y;
   vertex.position[2] = bc.z;
   mVertices.push_back(vertex);
}

void DebugDraw::AddAABB(math::Float4 min, math::Float4 max, unsigned int color)
{
}

void DebugDraw::Render(math::Matrix4 viewMatrix, math::Matrix4 projMatrix)
{

}
