
#include "utils.h"
#include "bvh.h"
#include "imgui.h"
#include "camera.h"

#include <DirectXMath.h>

struct MyUser
{
   int user = 0;
};

math::Matrix4 BVHCameraFrom;
math::Matrix4 BVHCameraTo;

DebugDraw* gDebugLines = nullptr;

void BVH_OnInit(ID3D11Device* device, ID3D11DeviceContext* context)
{
   gDebugLines = new DebugDraw(device, context);

   BVHCameraFrom.m[3] = math::Float4(10, 10, 10, 1);
   BVHCameraTo.m[3] = math::Float4(0, 0, 0, 1);
}

void BVH_OnShutdown()
{
   delete gDebugLines;
   gDebugLines = nullptr;
}

static bool BVHEdit = 0;
static bool BVHRebuild = 0;
static bool BVHOptimize = 0;
static bool BVHDraw = 0;

static int BVHDepth = -1;
static int BVHCount = 128;

static int BVHForceRefresh = -1;

static float BVHRebuildTime = 0;
static float BVHOptimizeTime = 0;
static float BVHInsertTime = 0;

math::Matrix4 QueryProjMatrix;
math::Matrix4 QueryViewMatrix;

bvh::BVH3D<MyUser, unsigned int> sBVHTest(BVHCount);

FreeCamera gFreeCamera;

void BVH_OnTick()
{
   class BVHDrawInterface : public bvh::DrawInterface
   {
   public:
      virtual void DrawBounds(const bvh::AABB& bounds, const bvh::AABB& parentBounds, int id, int depth, bool isLeaf) const
      {
         gDebugLines->AddAABB(bounds.min, bounds.max, 0xFFFFFFFF);
      }
   };

   gFreeCamera.Update();

   CameraMatrices cameraMatrices = gFreeCamera.ComputeMatrices(60.0f, 1.666f, 0.1, 10.0f);

   ImU32 gridCol = ImGui::GetColorU32(ImVec4(0.2f, 0.2f, 0.2f, 1));

   for (int i = -10; i <= 10; i++)
   {
      if (i == 0)
      {
         ImGui_DebugDrawLine_SLOW(cameraMatrices.projViewMatrix, math::Float4(-10, 0, i, 1), math::Float4(10, 0, i, 1), ImGui::GetColorU32(ImVec4(1, 0, 0, 1)));
         ImGui_DebugDrawLine_SLOW(cameraMatrices.projViewMatrix, math::Float4(i, 0, -10, 1), math::Float4(i, 0, 10, 1), ImGui::GetColorU32(ImVec4(0, 0, 1, 1)));
      }
      else
      {
         ImGui_DebugDrawLine_SLOW(cameraMatrices.projViewMatrix, math::Float4(-10, 0, i, 1), math::Float4(10, 0, i, 1), gridCol);
         ImGui_DebugDrawLine_SLOW(cameraMatrices.projViewMatrix, math::Float4(i, 0, -10, 1), math::Float4(i, 0, 10, 1), gridCol);
      }
   }

   ImGui_DebugDrawLine_SLOW(cameraMatrices.projViewMatrix, math::Float4(0, 0, 0, 1), math::Float4(0, 1, 0, 1), ImGui::GetColorU32(ImVec4(0, 1, 0, 1)));

   if (BVHForceRefresh)
   {
      BVHForceRefresh = false;

      sBVHTest = bvh::BVH3D<MyUser, unsigned int>(BVHCount);

      struct Data
      {
         bvh::AABB aabb;
         bvh::Handle handle;
      };

      std::vector<Data> data;
      data.resize(BVHCount);

      srand(0);

      for (int i = 0; i < BVHCount; i++)
      {
         auto random = [](float min, float max) -> float
         {
            float r = (rand() / (float)(RAND_MAX));
            return r * (max - min) + min;
         };

         math::Float4 center(random(-128, 128), random(-4, 4), random(-128, 128), 0);

         math::Float4 halfextents(random(0.1f, 8), random(0.1f, 8), random(0.1f, 8), 0);
         math::Float4 extentsmin = center.Sub(halfextents);
         math::Float4 extentsmax = center.Add(halfextents);

         data[i].aabb = bvh::AABB(extentsmin, extentsmax);
      }

      Timer insertTimer;

      int numSamples = 0;
      bool hasStarted = false;

      for (int i = 0; i < BVHCount; i++)
      {
         if (i == BVHCount - 100)
         {
            insertTimer.start();
            hasStarted = true;
         }

         data[i].handle = sBVHTest.Insert(data[i].aabb);
         
         if (hasStarted)
         {
            numSamples += 1;
         }
      }

      if (hasStarted)
      {
         insertTimer.stop();
         BVHInsertTime = (insertTimer.getMiliseconds()) / numSamples;
      }
      else
      {
         BVHInsertTime = 0;
      }

      if (BVHRebuild)
      {
         Timer rebuildTimer;
         rebuildTimer.start();
         sBVHTest.Rebuild();
         rebuildTimer.stop();
         BVHRebuildTime = rebuildTimer.getMiliseconds();
      }
      else
      {
         BVHRebuildTime = 0;
      }

      if (BVHOptimize)
      {
         Timer optimizeTimer;
         optimizeTimer.start();
         sBVHTest.Optimize();
         optimizeTimer.stop();
         BVHOptimizeTime = optimizeTimer.getMiliseconds();
      }
      else
      {
         BVHOptimizeTime = 0;
      }
   }
   
   bvh::BVH3D<MyUser, unsigned int>::QueryStats stats;
   Timer timer;


   {
      auto drawLine = [](math::Float4 a, math::Float4 b, unsigned int c)
      {
         auto aC = a.SplitComponents();
         auto bC = b.SplitComponents();

         gDebugLines->AddLine(a, b, c);
      };

      drawLine(BVHCameraFrom.m[3], BVHCameraTo.m[3], 0xFFFFFFFF);

      auto fromS = BVHCameraFrom.m[3].SplitComponents();
      auto toS = BVHCameraTo.m[3].SplitComponents();
      auto upS = math::Float4(0, 1, 0, 1);

      DirectX::XMMATRIX projMatrix = DirectX::XMMatrixPerspectiveFovLH(DirectX::XMConvertToRadians(45.0f), 1.0f, 0.1f, 10.0f);
      DirectX::XMMATRIX viewMatrix = DirectX::XMMatrixLookAtLH(BVHCameraFrom.m[3].m, BVHCameraTo.m[3].m, math::Float4(0, 1, 0, 1).m);

      DirectX::XMMATRIX projViewMatrix = DirectX::XMMatrixMultiply(viewMatrix, projMatrix);
      DirectX::XMMATRIX invProjViewMatrix = DirectX::XMMatrixInverse(nullptr, projViewMatrix);

      math::Matrix4 invProjViewMatrix2 = *(math::Matrix4*)&invProjViewMatrix;

      bvh::FrustumConvexHull<true> frustum;
      frustum.InitializeFromProjectionMatrix(invProjViewMatrix2);

      unsigned int color = 0xFFFFFFFF;

      drawLine(frustum.cornerPositions[0], frustum.cornerPositions[1], color);
      drawLine(frustum.cornerPositions[1], frustum.cornerPositions[3], color);
      drawLine(frustum.cornerPositions[3], frustum.cornerPositions[2], color);
      drawLine(frustum.cornerPositions[2], frustum.cornerPositions[0], color);

      drawLine(frustum.cornerPositions[4], frustum.cornerPositions[5], color);
      drawLine(frustum.cornerPositions[5], frustum.cornerPositions[7], color);
      drawLine(frustum.cornerPositions[7], frustum.cornerPositions[6], color);
      drawLine(frustum.cornerPositions[6], frustum.cornerPositions[4], color);

      drawLine(frustum.cornerPositions[0], frustum.cornerPositions[4], color);
      drawLine(frustum.cornerPositions[1], frustum.cornerPositions[5], color);
      drawLine(frustum.cornerPositions[3], frustum.cornerPositions[7], color);
      drawLine(frustum.cornerPositions[2], frustum.cornerPositions[6], color);

      timer.start();

      sBVHTest.Query(frustum.GetView(), [](const MyUser& data, const bvh::AABB& bounds, void* user)
      {
      }, nullptr, &stats);

      timer.stop();
      // 1,000,000 objects dynamically added
      // QueryTimer 13.344900ms f:27181 s:27900 l:360
      

      sBVHTest.Query(frustum.GetView(), [](const MyUser& data, const bvh::AABB& bounds, void* user)
      {
         BVHDrawInterface drawInterface;
         drawInterface.DrawBounds(bounds, bounds, 0, 0, true);
      }, nullptr, & stats);

   }


   BVHDrawInterface drawInterface;


#if 0
   auto view = (float*)&gRender->getViewD3DMatrix();
   auto proj = (float*)&gRender->getProjectionD3DMatrix();

   ImGuizmo::SetID(2);
   ImGuizmo::Manipulate(view, proj, ImGuizmo::TRANSLATE, ImGuizmo::WORLD, (float*)&BVHCameraFrom);
   ImGuizmo::SetID(1);
   ImGuizmo::Manipulate(view, proj, ImGuizmo::TRANSLATE, ImGuizmo::WORLD, (float*)&BVHCameraTo);
   ImGuizmo::SetID(0);
#endif

   ImGui::Begin("BVH");

   if (ImGui::Button("Reset Camera"))
   {
      BVHCameraFrom.m[3] = math::Float4(0, 0, 0, 1);
      BVHCameraTo.m[3] = math::Float4(10, 10, 10, 1);
   }

   if (ImGui::Button("Edit"))
   {
      BVHEdit = !BVHEdit;
   }

   //ImGuizmo::SetID(2);
   //ImGuizmo::Manipulate(view, proj, ImGuizmo::TRANSLATE, ImGuizmo::WORLD, (float*)&BVHCameraTo);
   //ImGuizmo::SetID(0);

   if (ImGui::SliderInt("Depth", &BVHDepth, -1, 32))
   {
      BVHForceRefresh = true;
   }

   if (ImGui::SliderInt("Count", &BVHCount, 0, 1024))
   {
      BVHForceRefresh = true;
   }

   if (ImGui::Checkbox("Rebuild", &BVHRebuild))
   {
      BVHForceRefresh = true;
   }

   if (ImGui::Checkbox("Optimize", &BVHOptimize))
   {
      BVHForceRefresh = true;
   }

   if (ImGui::Checkbox("Draw", &BVHDraw))
   {
   }

   static float QueryTime = 0;

   QueryTime = QueryTime * 0.5f + timer.getMiliseconds() * 0.05f;

   ImGui::Text("Cost: %f", sBVHTest.ComputeTotalCost());
   ImGui::Text("Query Timer %fms", QueryTime);
   ImGui::Text("Rebuild Timer %fms", BVHRebuildTime);
   ImGui::Text("Optimize Timer %fms", BVHOptimizeTime);
   ImGui::Text("Insert Timer %fus", BVHInsertTime * 1000);
   ImGui::Text("Failed Intersections: %i", stats.failedIntersections);
   ImGui::Text("Successful Intersections: %i", stats.successfulIntersections);
   ImGui::Text("Leaf Count: %i", stats.leafCount);

   ImGui::End();

   if (BVHDraw)
   {
      sBVHTest.Draw(&drawInterface);
   }


}
