#include "camera.h"

FreeCamera::FreeCamera()
{

}

FreeCamera::~FreeCamera()
{

}

void FreeCamera::Update()
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

DirectX::XMMATRIX FreeCamera::ComputeCameraRotationMatrix()
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

CameraMatrices FreeCamera::ComputeMatrices(float fovDegrees, float aspect, float nearPlane, float farPlane)
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
