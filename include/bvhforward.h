#pragma once

namespace bvh
{
   template<typename UserDataType>
   class BVH3D;

   struct Handle
   {
      int id = -1;

      bool IsValid()
      {
         return id != -1;
      }

      bool operator==(const Handle& other) const
      {
         return other.id == id;
      }

      bool operator!=(const Handle& other) const
      {
         return other.id != id;
      }
   };

   static const Handle INVALID_HANDLE = { -1 };
}
