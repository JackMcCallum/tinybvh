#pragma once

#include "utils.h"

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
   FreeCamera();
   ~FreeCamera();

   void Update();

   DirectX::XMMATRIX ComputeCameraRotationMatrix();

   CameraMatrices ComputeMatrices(float fovDegrees, float aspect, float nearPlane, float farPlane);

   math::Float4 mCameraPos = math::Float4(10, 10, 10, 1);
   float mCameraPitch = -45;
   float mCameraYaw = 45;
   float mSpeed = 1.0f;
};
