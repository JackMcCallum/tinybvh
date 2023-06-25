
#include "utils.h"
#include "bvh.h"
#include "imgui.h"

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

struct CameraMatrices
{
   // The matrix of the camera
   // AKA the matrix you would use to draw an object in worldspace where the camera is positioned
   DirectX::XMMATRIX cameraMatrix;

   // The view matrix, inverse of the camera matrix
   DirectX::XMMATRIX viewMatrix;

   // The projection matrix
   DirectX::XMMATRIX projMatrix;

   // Premultiplied matrices
   DirectX::XMMATRIX projViewMatrix;
   DirectX::XMMATRIX invProjViewMatrix;
};

class FreeCamera
{
public:
   FreeCamera()
   {
      
   }

   ~FreeCamera()
   {

   }

   void Update()
   {
      if (!ImGui::GetIO().WantCaptureMouse)
      {
         if (ImGui::GetIO().MouseDown[ImGuiMouseButton_Right])
         {
            // Dragging the mouse UP leads to a negative delta, need to look up which is positive pitch so we flip it
            mCameraPitch -= ImGui::GetIO().MouseDelta.y * 0.1f;
            mCameraYaw += ImGui::GetIO().MouseDelta.x * 0.1f;

            float speed = 0.1f;

            DirectX::XMMATRIX cameraRot = ComputeCameraRotationMatrix();

            if (ImGui::GetIO().KeysDown[ImGuiKey_W])
            {
               mCameraPos.m = DirectX::XMVectorAdd(mCameraPos.m, DirectX::XMVector3Transform(math::Float4(0, 0, -speed, 0).m, cameraRot));
            }

            if (ImGui::GetIO().KeysDown[ImGuiKey_S])
            {
               mCameraPos.m = DirectX::XMVectorAdd(mCameraPos.m, DirectX::XMVector3Transform(math::Float4(0, 0, speed, 0).m, cameraRot));
            }

            if (ImGui::GetIO().KeysDown[ImGuiKey_D])
            {
               mCameraPos.m = DirectX::XMVectorAdd(mCameraPos.m, DirectX::XMVector3Transform(math::Float4(speed, 0, 0, 0).m, cameraRot));
            }

            if (ImGui::GetIO().KeysDown[ImGuiKey_A])
            {
               mCameraPos.m = DirectX::XMVectorAdd(mCameraPos.m, DirectX::XMVector3Transform(math::Float4(-speed, 0, 0, 0).m, cameraRot));
            }

            if (ImGui::GetIO().KeysDown[ImGuiKey_Q])
            {
               mCameraPos.m = DirectX::XMVectorAdd(mCameraPos.m, DirectX::XMVector3Transform(math::Float4(0, speed, 0, 0).m, cameraRot));
            }

            if (ImGui::GetIO().KeysDown[ImGuiKey_E])
            {
               mCameraPos.m = DirectX::XMVectorAdd(mCameraPos.m, DirectX::XMVector3Transform(math::Float4(0, -speed, 0, 0).m, cameraRot));
            }
         }
      }
   }

   DirectX::XMMATRIX ComputeCameraRotationMatrix()
   {
      DirectX::XMMATRIX cameraRot = DirectX::XMMatrixRotationRollPitchYaw(
         DirectX::XMConvertToRadians(mCameraPitch),
         DirectX::XMConvertToRadians(mCameraYaw),
         DirectX::XMConvertToRadians(0.0f));

      float yaw = DirectX::XMConvertToRadians(mCameraYaw);
      float pitch = DirectX::XMConvertToRadians(mCameraPitch);

      float x = cos(yaw) * cos(pitch);
      float y = -sin(pitch); // Negative pitch looks down
      float z = sin(yaw) * cos(pitch);

      return DirectX::XMMatrixInverse(nullptr, DirectX::XMMatrixLookToLH(math::ZERO.m, math::Float4(x, y, z, 1).m, math::Float4(0, 1, 0, 1).m));
   }

   CameraMatrices ComputeMatrices(float fovDegrees, float aspect, float nearPlane, float farPlane)
   {

      auto pos = mCameraPos.SplitComponents();
      DirectX::XMMATRIX cameraPos = DirectX::XMMatrixTranslation(pos.x, pos.y, pos.z);

      CameraMatrices out;
      out.cameraMatrix = DirectX::XMMatrixMultiply(ComputeCameraRotationMatrix(), cameraPos);
      out.viewMatrix = DirectX::XMMatrixInverse(nullptr, out.cameraMatrix);
      out.projMatrix = DirectX::XMMatrixPerspectiveFovRH(DirectX::XMConvertToRadians(fovDegrees), aspect, nearPlane, farPlane);
      out.projViewMatrix = DirectX::XMMatrixMultiply(out.viewMatrix, out.projMatrix);
      out.invProjViewMatrix = DirectX::XMMatrixInverse(nullptr, out.projViewMatrix);
      return out;
   }

   math::Float4 mCameraPos = math::Float4(10, 10, 10, 1);
   float mCameraPitch = -45;
   float mCameraYaw = 45;
};

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

   auto imguiDrawLine3D = [&cameraMatrices](math::Float4 a, math::Float4 b, ImU32 col)
   {
      ImDrawList* drawList = ImGui::GetBackgroundDrawList();

      math::Float4 at = DirectX::XMVector4Transform(a.m, cameraMatrices.projViewMatrix);
      math::Float4 bt = DirectX::XMVector4Transform(b.m, cameraMatrices.projViewMatrix);

      auto ac = at.SplitComponents();
      auto bc = bt.SplitComponents();

      // Both behind the camera
      if (ac.z < 0 && bc.z < 0)
      {
         return;
      }

      if (ac.z < 0)
      {
         // ac is behind the camera
         float s = (bc.z / (ac.z - bc.z)) * -0.999f;
         at = bt.Add(at.Sub(bt).Mul(s));
      }
      else if (bc.z < 0)
      {
         // bc is behind the camera
         float s = (ac.z / (bc.z - ac.z)) * -0.999f;
         bt = at.Add(bt.Sub(at).Mul(s));
      }

      // Perspective divide, then scale from (-1 to 1) -> (0 to 1)
      at = at.Div(at.Shuffle<3>()).Mul(math::Float4(0.5f, -0.5f, 0.5f, 1.0f)).Add(0.5f);
      bt = bt.Div(bt.Shuffle<3>()).Mul(math::Float4(0.5f, -0.5f, 0.5f, 1.0f)).Add(0.5f);

      math::Float4 displaySize(ImGui::GetIO().DisplaySize.x, ImGui::GetIO().DisplaySize.y, 1, 1);

      ac = at.Mul(displaySize).SplitComponents();
      bc = bt.Mul(displaySize).SplitComponents();

      drawList->AddLine(ImVec2(ac.x, ac.y), ImVec2(bc.x, bc.y), col);
   };

   ImU32 gridCol = ImGui::GetColorU32(ImVec4(0.2f, 0.2f, 0.2f, 1));

   for (int i = -10; i <= 10; i++)
   {
      if (i == 0)
      {
         imguiDrawLine3D(math::Float4(-10, 0, i, 1), math::Float4(10, 0, i, 1), ImGui::GetColorU32(ImVec4(1, 0, 0, 1)));
         imguiDrawLine3D(math::Float4(i, 0, -10, 1), math::Float4(i, 0, 10, 1), ImGui::GetColorU32(ImVec4(0, 0, 1, 1)));
      }
      else
      {
         imguiDrawLine3D(math::Float4(-10, 0, i, 1), math::Float4(10, 0, i, 1), gridCol);
         imguiDrawLine3D(math::Float4(i, 0, -10, 1), math::Float4(i, 0, 10, 1), gridCol);
      }
   }

   imguiDrawLine3D(math::Float4(0, 0, 0, 1), math::Float4(0, 1, 0, 1), ImGui::GetColorU32(ImVec4(0, 1, 0, 1)));

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
