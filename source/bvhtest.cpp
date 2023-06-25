
#include "utils.h"
#include "bvh.h"
#include "imgui.h"
#include "camera.h"

#include <DirectXMath.h>

struct MyUser
{
   int user = 0;
};

DebugDraw* gDebugLines = nullptr;

static bool BVHEdit = 0;
static bool BVHRebuild = 0;
static bool BVHOptimize = 0;
static bool BVHDraw = 0;

static int BVHDepth = -1;
static int BVHNumStaticObjects = 128;
static int BVHNumDynamicObjects = 128;

static int BVHForceRefresh = -1;

static float BVHRebuildTime = 0;
static float BVHOptimizeTime = 0;
static float BVHInsertTime = 0;
static float BVHQueryTime = 0;
static float BVHUpdateTime = 0;

static float BVHAreaSize = 128;
static float BVHAreaHeight = 4;

math::Matrix4 QueryProjMatrix;
math::Matrix4 QueryViewMatrix;

struct Object
{
   math::Float4 position;
   math::Float4 halfextents;
   bvh::Handle handle;

   math::Float4 velocity;
   bool isStatic;
};

std::vector<Object> Objects;

bvh::BVH3D<MyUser, unsigned int> sBVHTest;
bvh::BVH3D<MyUser, unsigned int>::QueryStats sBVHStats;

FreeCamera gFreeCamera;
CameraMatrices gCameraMatrices;
CameraMatrices gQueryMatrices;

void UpdateCamera()
{
   gFreeCamera.Update();
   gCameraMatrices = gFreeCamera.ComputeMatrices(60.0f, 1.666f, 0.1, 100.0f);
}


void BVH_OnInit(ID3D11Device* device, ID3D11DeviceContext* context)
{
   gDebugLines = new DebugDraw(device, context);

   UpdateCamera();
   gQueryMatrices = gCameraMatrices;
}

void BVH_OnShutdown()
{
   delete gDebugLines;
   gDebugLines = nullptr;
}

void DrawGrid()
{
   ImU32 gridCol = ImGui::GetColorU32(ImVec4(0.2f, 0.2f, 0.2f, 1));

   for (int i = -10; i <= 10; i++)
   {
      if (i == 0)
      {
         ImGui_DebugDrawLine_SLOW(gCameraMatrices.projViewMatrix, math::Float4(-10, 0, i, 1), math::Float4(10, 0, i, 1), ImGui::GetColorU32(ImVec4(1, 0, 0, 1)));
         ImGui_DebugDrawLine_SLOW(gCameraMatrices.projViewMatrix, math::Float4(i, 0, -10, 1), math::Float4(i, 0, 10, 1), ImGui::GetColorU32(ImVec4(0, 0, 1, 1)));
      }
      else
      {
         ImGui_DebugDrawLine_SLOW(gCameraMatrices.projViewMatrix, math::Float4(-10, 0, i, 1), math::Float4(10, 0, i, 1), gridCol);
         ImGui_DebugDrawLine_SLOW(gCameraMatrices.projViewMatrix, math::Float4(i, 0, -10, 1), math::Float4(i, 0, 10, 1), gridCol);
      }
   }

   ImGui_DebugDrawLine_SLOW(gCameraMatrices.projViewMatrix, math::Float4(0, 0, 0, 1), math::Float4(0, 1, 0, 1), ImGui::GetColorU32(ImVec4(0, 1, 0, 1)));
}

void OnGUI()
{
   ImGui::Begin("BVH");

   if (ImGui::Button("Reset Camera"))
   {
      gQueryMatrices = gCameraMatrices;
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

   if (ImGui::SliderInt("Static Objects", &BVHNumStaticObjects, 0, 1024))
   {
      BVHForceRefresh = true;
   }

   if (ImGui::SliderInt("Dynamic Objects", &BVHNumDynamicObjects, 0, 1024))
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

   ImGui::Text("Cost: %.3f", sBVHTest.ComputeTotalCost());
   ImGui::Text("Query Timer %.3fms", BVHQueryTime);
   ImGui::Text("Rebuild Timer %.3fms", BVHRebuildTime);
   ImGui::Text("Optimize Timer %.3fms", BVHOptimizeTime);
   ImGui::Text("Update Timer %.3fms %.3fus/Object", BVHUpdateTime, (BVHUpdateTime / BVHNumDynamicObjects) * 1000);
   ImGui::Text("Insert Timer %.3fus", BVHInsertTime * 1000);

   ImGui::Text("Failed Intersections: %i", sBVHStats.failedIntersections);
   ImGui::Text("Successful Intersections: %i", sBVHStats.successfulIntersections);
   ImGui::Text("Leaf Count: %i", sBVHStats.leafCount);

   ImGui::End();

}

void RebuildBVH()
{
   sBVHTest = bvh::BVH3D<MyUser, unsigned int>(BVHNumStaticObjects + BVHNumDynamicObjects);

   Objects.resize(BVHNumStaticObjects + BVHNumDynamicObjects);

   srand(0);

   for (int i = 0; i < Objects.size(); i++)
   {
      auto random = [](float min, float max) -> float
      {
         float r = (rand() / (float)(RAND_MAX));
         return r * (max - min) + min;
      };

      math::Float4 position(random(-BVHAreaSize, BVHAreaSize), random(-BVHAreaHeight, BVHAreaHeight), random(-BVHAreaSize, BVHAreaSize), 0);
      math::Float4 halfextents(random(0.1f, 8), random(0.1f, 8), random(0.1f, 8), 0);

      Objects[i].position = position;
      Objects[i].halfextents = halfextents;

      if (i >= BVHNumStaticObjects)
      {
         float v = 0.1f;
         math::Float4 velocity(random(-v, v), random(-v, v), random(-v, v), 0);
         Objects[i].velocity = velocity;
         Objects[i].isStatic = false;
      }
      else
      {
         Objects[i].isStatic = true;
      }
   }

   Timer insertTimer;

   int numSamples = 0;
   bool hasStarted = false;

   for (int i = 0; i < Objects.size(); i++)
   {
      if (i == Objects.size() - 100)
      {
         insertTimer.start();
         hasStarted = true;
      }

      bvh::AABB bounds(
         Objects[i].position.Sub(Objects[i].halfextents),
         Objects[i].position.Add(Objects[i].halfextents));

      Objects[i].handle = sBVHTest.Insert(bounds);

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


void BVH_OnTick()
{
   auto drawLine = [](math::Float4 a, math::Float4 b, unsigned int c)
   {
      ImGui_DebugDrawLine_SLOW(gCameraMatrices.projViewMatrix, a, b, c);
   };

   class BVHDrawInterface : public bvh::DrawInterface
   {
   public:
      virtual void DrawBounds(const bvh::AABB& bounds, const bvh::AABB& parentBounds, int id, int depth, bool isLeaf) const
      {
         ImGui_DebugDrawAABB_SLOW(gCameraMatrices.projViewMatrix, bounds.min, bounds.max, ImGui::GetColorU32(ImVec4(0, 0, 0, 0.1f)));
      }
   };

   UpdateCamera();

   if (BVHForceRefresh)
   {
      BVHForceRefresh = false;
      RebuildBVH();
   }

   float deltaTime = 60.0f / 1000.0f;

   Timer updateObjectsTimer;

   for (int i = 0; i < Objects.size(); i++)
   {
      if (Objects[i].isStatic == false)
      {
         auto pos = Objects[i].position;
         auto vel = Objects[i].velocity;

         auto spos = pos.SplitComponents();
         auto svel = vel.SplitComponents();

         if ((spos.x + svel.x) < -BVHAreaSize || (spos.x + svel.x) > BVHAreaSize)
         {
            svel.x *= -1;
         }

         if ((spos.y + svel.y) < -BVHAreaHeight || (spos.y + svel.y) > BVHAreaHeight)
         {
            svel.y *= -1;
         }

         if ((spos.z + svel.z) < -BVHAreaSize || (spos.z + svel.z) > BVHAreaSize)
         {
            svel.z *= -1;
         }

         Objects[i].velocity = svel;
         Objects[i].position = pos.Add(svel);

         bvh::AABB bounds(
            Objects[i].position.Sub(Objects[i].halfextents),
            Objects[i].position.Add(Objects[i].halfextents));

         sBVHTest.Update(Objects[i].handle, bounds);
      }
   }

   BVHUpdateTime = BVHUpdateTime * 0.5f + updateObjectsTimer.getMiliseconds() * 0.05f;

   // Draw the grid
   {
      ImU32 gridCol = ImGui::GetColorU32(ImVec4(0.2f, 0.2f, 0.2f, 1));

      for (int i = -10; i <= 10; i++)
      {
         if (i == 0)
         {
            drawLine(math::Float4(-10, 0, i, 1), math::Float4(10, 0, i, 1), ImGui::GetColorU32(ImVec4(1, 0, 0, 1)));
            drawLine(math::Float4(i, 0, -10, 1), math::Float4(i, 0, 10, 1), ImGui::GetColorU32(ImVec4(0, 0, 1, 1)));
         }
         else
         {
            drawLine(math::Float4(-10, 0, i, 1), math::Float4(10, 0, i, 1), gridCol);
            drawLine(math::Float4(i, 0, -10, 1), math::Float4(i, 0, 10, 1), gridCol);
         }
      }

      drawLine(math::Float4(0, 0, 0, 1), math::Float4(0, 1, 0, 1), ImGui::GetColorU32(ImVec4(0, 1, 0, 1)));
   }


   Timer timer;

   {

      auto upS = math::Float4(0, 1, 0, 1);

      math::Matrix4 invProjViewMatrix2 = *(math::Matrix4*)&gQueryMatrices.invProjViewMatrix;

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

      // Query and do nothing to get a measurement of just the query
      sBVHTest.Query(frustum.GetView(), [](const MyUser& data, const bvh::AABB& bounds, void* user)
      {
      }, nullptr, &sBVHStats);

      timer.stop();
      
      BVHQueryTime = BVHQueryTime * 0.5f + timer.getMiliseconds() * 0.05f;

      // Query again but this time draw
      sBVHTest.Query(frustum.GetView(), [](const MyUser& data, const bvh::AABB& bounds, void* user)
      {
         ImGui_DebugDrawAABB_SLOW(gCameraMatrices.projViewMatrix, bounds.min, bounds.max, ImGui::GetColorU32(ImVec4(0, 1, 0, 1.0f)));
      }, nullptr, &sBVHStats);
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

   OnGUI();

   if (BVHDraw)
   {
      sBVHTest.Draw(&drawInterface);
   }


}
