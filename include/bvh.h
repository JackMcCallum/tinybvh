/*
TODO LIST:
Implement ray casts
Box queries
Orthographic queries
*/

#pragma once
#include <stdio.h>
#include <vector>
#include <array>
#include <algorithm>
#include <queue>
#include <xmmintrin.h>

#ifndef BVH_ASSERT
#define BVH_ASSERT(cond) do { if (!(cond)) { printf("ASSERT FAILED: %s", #cond); __debugbreak(); } } while (false)
#endif

#ifndef BVH_ASSERTF
#define BVH_ASSERTF(cond, ...) do { if (!(cond)) { printf("ASSERT FAILED: %s - %s", #cond, __VA_ARGS__); __debugbreak(); } } while (false)
#endif

#ifndef BVH_FAIL
#define BVH_FAIL(...) do { printf("ASSERT FAILED: %s", __VA_ARGS__); __debugbreak(); } while (false)
#endif

#ifndef USE_VALIDATION
#define USE_VALIDATION 0
#endif

#define BVH_FORCEINLINE __forceinline

#undef min
#undef max

namespace math
{
	struct Float4
	{
		struct Components
		{
			float x;
			float y;
			float z;
			float w;

			const float& operator[](int i) const
			{
				BVH_ASSERT(i >= 0 && i < 4);
				return (&x)[i];
			}

			float& operator[](int i)
			{
				BVH_ASSERT(i >= 0 && i < 4);
				return (&x)[i];
			}
		};

		BVH_FORCEINLINE Float4()
		{
			__declspec(align(16)) float tmp[4] = { 0, 0, 0, 0 };
			m = _mm_load_ps(tmp);
		}

		BVH_FORCEINLINE Float4(const Components& c)
		{
			__declspec(align(16)) float tmp[4] = { c.x, c.y, c.z, c.w };
			m = _mm_load_ps(tmp);
		}

		BVH_FORCEINLINE Float4(__m128 _m)
		{
			m = _m;
		}

		BVH_FORCEINLINE Float4(float x, float y, float z, float w)
		{
			__declspec(align(16)) float tmp[4] = { x, y, z, w };
			m = _mm_load_ps(tmp);
		}

		BVH_FORCEINLINE Float4 Add(float rhs)
		{
			return _mm_add_ps(m, Float4(rhs, rhs, rhs, rhs).m);
		}

		BVH_FORCEINLINE Float4 Add(Float4 rhs) const
		{
			return _mm_add_ps(m, rhs.m);
		}

		BVH_FORCEINLINE Float4 Sub(float rhs) const
		{
			return _mm_sub_ps(m, Float4(rhs, rhs, rhs, rhs).m);
		}

		BVH_FORCEINLINE Float4 Sub(Float4 rhs) const
		{
			return _mm_sub_ps(m, rhs.m);
		}

		BVH_FORCEINLINE Float4 Mul(float rhs) const
		{
			return _mm_mul_ps(m, Float4(rhs, rhs, rhs, rhs).m);
		}

		BVH_FORCEINLINE Float4 Mul(Float4 rhs) const
		{
			return _mm_mul_ps(m, rhs.m);
		}

		BVH_FORCEINLINE Float4 Div(float rhs) const
		{
			return _mm_div_ps(m, Float4(rhs, rhs, rhs, rhs).m);
		}

		BVH_FORCEINLINE Float4 Div(Float4 rhs) const
		{
			return _mm_div_ps(m, rhs.m);
		}

		BVH_FORCEINLINE Float4 Min(Float4 rhs) const
		{
			return _mm_min_ps(m, rhs.m);
		}

		BVH_FORCEINLINE Float4 Max(Float4 rhs) const
		{
			return _mm_max_ps(m, rhs.m);
		}

		// Compute the dot product on XYZ components
		BVH_FORCEINLINE float Dot3(Float4 rhs) const
		{
			__m128 mult = _mm_mul_ps(m, rhs.m);

			float result[4];
			_mm_store_ps(result, mult);

			return result[0] + result[1] + result[2];
		}

		// Compute the dot product on XYZW components
		BVH_FORCEINLINE float Dot4(Float4 rhs) const
		{
			__m128 mult = _mm_mul_ps(m, rhs.m);

			float result[4];
			_mm_store_ps(result, mult);

			return result[0] + result[1] + result[2] + result[3];
		}

		// Compute the dot product on XYZW components
		BVH_FORCEINLINE Float4 Cross(Float4 rhs) const
		{
			const Float4& lhs = *this;

			Float4 a = lhs.Shuffle<1, 2, 0, 3>();
			Float4 b = rhs.Shuffle<2, 0, 1, 3>();
			Float4 c = lhs.Shuffle<2, 0, 1, 3>();
			Float4 d = rhs.Shuffle<1, 2, 0, 3>();

			Float4 l = a.Mul(b);
			Float4 r = c.Mul(d);

			return l.Sub(r);
		}

		// Compute the dot product on XYZW components
		BVH_FORCEINLINE Components SplitComponents() const
		{
			Components components;
			_mm_store_ps((float*)&components, m);
			return components;
		}

		template<int X, int Y, int Z, int W>
		BVH_FORCEINLINE Float4 Shuffle() const
		{
			// Note, flipped for personal sanity
			// I dont understand why this macro is flipped this way
			// Something to do with the order of components vs the order of bits i think?
			// It makes sense to flip it, so Shuffle<0, 0, 3, 2>() would result in xxwz and !zwxx
			return _mm_shuffle_ps(m, m, _MM_SHUFFLE(W, Z, Y, X));
		}

		template<int N>
		BVH_FORCEINLINE Float4 Shuffle() const
		{
			return _mm_shuffle_ps(m, m, _MM_SHUFFLE(N, N, N, N));
		}

		__m128 m;
	};

	struct Matrix4
	{
		Matrix4()
		{
			// Init to identity
			m[0] = Float4(1, 0, 0, 0);
			m[1] = Float4(0, 1, 0, 0);
			m[2] = Float4(0, 0, 1, 0);
			m[3] = Float4(0, 0, 0, 1);
		}

		Float4 m[4];

		// Multiply vector by matrix
		Float4 Multiply(const Float4& in) const
		{
			Float4 xP = in.Shuffle<0>();
			Float4 yP = in.Shuffle<1>();
			Float4 zP = in.Shuffle<2>();
			Float4 wP = in.Shuffle<3>();

			Float4 xC = xP.Mul(m[0]);
			Float4 yC = yP.Mul(m[1]);
			Float4 zC = zP.Mul(m[2]);
			Float4 wC = wP.Mul(m[3]);

			return xC.Add(yC).Add(zC).Add(wC);
		}
	};

	static const Float4 HALF = Float4(0.5f, 0.5f, 0.5f, 0.5f);
	static const Float4 ZERO = Float4(0, 0, 0, 0);
	static const Float4 ONE = Float4(1, 1, 1, 1);
	static const Float4 LARGEST_POSITIVE = Float4(FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX);
	static const Float4 LARGEST_NEGATIVE = Float4(-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX);
	static const Float4 SMALLEST_POSITIVE = Float4(FLT_MIN, FLT_MIN, FLT_MIN, FLT_MIN);
	static const Float4 SMALLEST_NEGATIVE = Float4(-FLT_MIN, -FLT_MIN, -FLT_MIN, -FLT_MIN);

	static const Float4 AXES_VECTORS[4] = {
	   Float4(1, 0, 0, 0),
	   Float4(0, 1, 0, 0),
	   Float4(0, 0, 1, 0),
	   Float4(0, 0, 0, 1),
	};
}

namespace bvh
{
	struct Handle
	{
		int id = -1;
	};

	class Timer
	{
	public:
		static double GetTimeNowAsSeconds()
		{
			static LARGE_INTEGER frequency;

			static bool init = false;
			if (!init)
			{
				QueryPerformanceFrequency(&frequency);
				init = true;
			}

			LARGE_INTEGER currentTime;
			QueryPerformanceCounter(&currentTime);

			return ((double)currentTime.QuadPart / (double)frequency.QuadPart);
		}

		Timer()
		{
			start();
		}

		void start()
		{
			mStartTime = GetTimeNowAsSeconds();
		}

		void stop()
		{
			mStopTime = GetTimeNowAsSeconds();
		}

		double getSeconds()
		{
			double endTime = 0;
			if (mStopTime == 0)
			{
				endTime = GetTimeNowAsSeconds();
			}
			else
			{
				endTime = mStopTime;
			}

			return endTime - mStartTime;
		}

		double getMiliseconds()
		{
			return getSeconds() * 1000;
		}

		double getMicroseconds()
		{
			return getSeconds() * 1000000;
		}

		double mStartTime = 0;
		double mStopTime = 0;

	};

	// Adds a layer of safety to numeric type casting
	// Protects against accidental data loss use this instead of cstyle casting for basic numeric types
	template<typename Dst, typename Src>
	inline Dst NumericCast(Src value)
	{
#if USE_ASSERT
		typedef std::numeric_limits<Dst> DstLim;
		typedef std::numeric_limits<Src> SrcLim;

		// Disable unsigned signed mismatch warning, this is intention in this function
#pragma warning( push )
#pragma warning( disable: 4018 )

		const bool positiveOverflowPossible = DstLim::max() < SrcLim::max();
		const bool negativeOverflowPossible = SrcLim::is_signed or (DstLim::lowest() > SrcLim::lowest());

		if ((!DstLim::is_signed) and (!SrcLim::is_signed))
		{
			if (positiveOverflowPossible and (value > DstLim::max()))
			{
				FAIL("positive overflow");
			}
		}
		else if ((!DstLim::is_signed) and SrcLim::is_signed)
		{
			if (positiveOverflowPossible and (value > DstLim::max()))
			{
				FAIL("positive overflow");
			}
			else if (negativeOverflowPossible and (value < 0))
			{
				FAIL("negative overflow");
			}
		}
		else if (DstLim::is_signed and (!SrcLim::is_signed))
		{
			if (positiveOverflowPossible and (value > DstLim::max()))
			{
				FAIL("positive overflow");
			}
		}
		else if (DstLim::is_signed and SrcLim::is_signed)
		{
			if (positiveOverflowPossible and (value > DstLim::max()))
			{
				FAIL("positive overflow");
			}
			else if (negativeOverflowPossible and (value < DstLim::lowest()))
			{
				FAIL("negative overflow");
			}
		}

#pragma warning( pop )
#endif

		// limits have been checked, therefore safe to cast
		return static_cast<Dst>(value);
	}

	static const Handle INVALID_HANDLE = { -1 };

	template<typename T>
	class MemoryView
	{
	public:
		MemoryView(const T* memory, int count)
		{
			mMemory = memory;
			mCount = count;
		}

		const T& operator[](int idx) const
		{
			assert(idx >= 0 && idx < mCount);
			return mMemory[idx];
		}

		T& operator[](int idx)
		{
			assert(idx >= 0 && idx < mCount);
			return mMemory[idx];
		}

		int Count() const { return mCount; }

		const T* begin() const { return &mMemory[0]; }
		T* begin() { return &mMemory[0]; }

		const T* end() const { return &mCount[0]; }
		T* end() { return &mCount[0]; }

	private:
		const T* mMemory;
		int mCount;
	};

	template<typename T>
	class ConstMemoryView
	{
	public:
		ConstMemoryView(const T* memory, int count)
		{
			mMemory = memory;
			mCount = count;
		}

		const T& operator[](int idx) const
		{
			BVH_ASSERT(idx >= 0 && idx < mCount);
			return mMemory[idx];
		}

		int Count() const { return mCount; }

		const T* begin() const { return &mMemory[0]; }
		const T* end() const { return &mMemory[mCount]; }

	private:
		const T* mMemory = nullptr;
		int mCount = 0;
	};

	// A view into a convex hull, allowing statically or dynamically allocated, varying sizes to be used
	class ConvexHullView
	{
	public:
		ConvexHullView(
			ConstMemoryView<math::Float4> _edgeAxisVectors,
			ConstMemoryView<math::Float4> _faceAxisVectors,
			ConstMemoryView<math::Float4> _facePlanes,
			ConstMemoryView<math::Float4> _permutedVertices) :
			edgeAxisVectors(_edgeAxisVectors),
			faceAxisVectors(_faceAxisVectors),
			facePlanes(_facePlanes),
			permutedVertices(_permutedVertices)
		{}

		// Edge vectors, stored (XYZ_, XYZ_, XYZ_, ...)
		// These are stored to generate new edges from the other hull when a test needs to be performed
		// They need to be in the most optimized state for generating cross products quickly
		ConstMemoryView<math::Float4> edgeAxisVectors;

		// Face axis vectors, stored (XYZD, XYZD, XYZD, ...)
		// these are the vectors that get tested first since they can be generated ahead of time independant of other hulls
		ConstMemoryView<math::Float4> faceAxisVectors;

		// Face planes, not needed for separating axis theorem but can be very handy for performing some early optimizations
		// We could test if only some of the vertices from one hull are inside the other hull and bypass the SAT test
		// We could also test if all the vertices of one hull are inside another hull, allowing us to iterate all children nodes and bypass their SAT tests
		ConstMemoryView<math::Float4> facePlanes;

		// Permuted vertices, stored (xxxx, yyyy, zzzz, xxxx, yyyy, zzzz, etc...)
		// These are tested against all the face axes of the two convex hulls
		// If the narrow phase is to be performed, its tested against the generated edge axes
		ConstMemoryView<math::Float4> permutedVertices;
	};

	// Check if a sphere is fully iside, fully outside, or intersecting with the hull
	// This is not true sphere/hull accuracy as there will be inaccuracies around edges of the hull
	inline int HullSphereIntersection(const ConvexHullView& hull, const math::Float4& center, float radius)
	{
		int numPlanes = hull.facePlanes.Count();

		int numFullyInside = 0;

		for (int p = 0; p < numPlanes; p++)
		{
			float signedDistance = hull.facePlanes[p].Dot3(center) - hull.facePlanes[p].SplitComponents().w;

			if (signedDistance < -radius)
			{
				numFullyInside++;
			}
			else if (signedDistance > radius)
			{
				// Fully outside any plane means we are fully outside, and can early out
				return 1;
			}
		}

		// Fully outside?
		if (numFullyInside == numPlanes)
		{
			return -1;
		}
		
		return 0;
	}

	template<int MAX_VERTICES = 8>
	inline int CountVerticesInsideHull(const ConvexHullView& hull, const ConstMemoryView<math::Float4>& vertexList)
	{
		int numIterations = (vertexList.Count() + 2) / 3;
		int numPlanes = hull.facePlanes.Count();
		int numInside = 0;

		std::array<__m128i, MAX_VERTICES / 4> numPlanesInside;

		__declspec(align(16)) int int1V[4] = { 1, 1, 1, 1 };
		__m128i int1 = _mm_load_si128((__m128i*)int1V);

		__declspec(align(16)) int numPlanesV[4] = { numPlanes, numPlanes, numPlanes, numPlanes };
		__m128i compare = _mm_load_si128((__m128i*)numPlanesV);

		__declspec(align(16)) int zeroV[4] = { 0, 0, 0, 0 };
		__m128i zero = _mm_load_si128((__m128i*)zeroV);

		for (int i = 0; i < numIterations; i++)
		{
			numPlanesInside[i] = zero;
		}

		for (int p = 0; p < numPlanes; p++)
		{
			// Permute the planes
			math::Float4 planeX = hull.facePlanes[p].Shuffle<0>();
			math::Float4 planeY = hull.facePlanes[p].Shuffle<1>();
			math::Float4 planeZ = hull.facePlanes[p].Shuffle<2>();
			math::Float4 planeD = hull.facePlanes[p].Shuffle<3>();

			for (int i = 0; i < numIterations; i++)
			{
				// XYZ components permuted into vectors
				const math::Float4& x = vertexList[i * 3 + 0];
				const math::Float4& y = vertexList[i * 3 + 1];
				const math::Float4& z = vertexList[i * 3 + 2];

				math::Float4 xx = x.Mul(planeX);
				math::Float4 yy = y.Mul(planeY);
				math::Float4 zz = z.Mul(planeZ);
				math::Float4 dots = xx.Add(yy).Add(zz);
				math::Float4 signedDistances = dots.Sub(planeD);

				// Compare signedDistances less than zero
				__m128i lt = _mm_castps_si128(_mm_cmplt_ps(signedDistances.m, _mm_setzero_ps()));

				__m128i result = _mm_and_si128(lt, int1);
				numPlanesInside[i] = _mm_add_epi32(numPlanesInside[i], result);
			}
		}

		__m128i cmpResult = zero;

		for (int i = 0; i < numIterations; i++)
		{
			cmpResult = _mm_add_epi32(cmpResult, _mm_and_si128(int1, _mm_cmpeq_epi32(compare, numPlanesInside[i])));
		}

		__declspec(align(16)) int cmpResultV[4];
		_mm_store_si128((__m128i*)cmpResultV, cmpResult);

		numInside += cmpResultV[0];
		numInside += cmpResultV[1];
		numInside += cmpResultV[2];
		numInside += cmpResultV[3];

		if (numInside == 0)
		{
			return 1; // 1 indicates all vertices are ouside
		}
		else if (numInside == numIterations * 4)
		{
			return -1; // -1 indicates all vertices are inside
		}

		return 0; // Zero indicates some are inside and outside
	}

	inline void FindMinMaxVertices(const math::Float4& axis, const ConstMemoryView<math::Float4>& vertexList, float& min, float& max)
	{
		min = FLT_MAX;
		max = -FLT_MAX;

		int numIterations = (vertexList.Count() + 2) / 3;

		math::Float4 axisX = axis.Shuffle<0>();
		math::Float4 axisY = axis.Shuffle<1>();
		math::Float4 axisZ = axis.Shuffle<2>();

		for (int i = 0; i < numIterations; i++)
		{
			// XYZ components permuted into vectors
			const math::Float4& x = vertexList[i * 3 + 0];
			const math::Float4& y = vertexList[i * 3 + 1];
			const math::Float4& z = vertexList[i * 3 + 2];

			// Lets compute the dot products for all 4 at once
			math::Float4 xx = x.Mul(axisX);
			math::Float4 yy = y.Mul(axisY);
			math::Float4 zz = z.Mul(axisZ);

			// Sum them up
			math::Float4 sum = xx.Add(yy).Add(zz);

			math::Float4::Components components = sum.SplitComponents();

			//int numVertices = std::min(4, vertexList.Count() - i * 4);
			for (int j = 0; j < 4; j++)
			{
				min = std::min(components[j], min);
				max = std::max(components[j], max);
			}
		}
	}

	inline bool IntersectsWith(const ConvexHullView& hullA, const ConvexHullView& hullB, bool performEdgeTests)
	{
		// First test all of the plane axes since these planes are already precalculated
		// If we can conclude separation here then we get to early out

		// Plane axis separation
		for (int i = 0; i < 2; i++)
		{
			const auto& axisList = (i == 0) ? hullA.faceAxisVectors : hullB.faceAxisVectors;
			for (const auto& axis : axisList)
			{
				float minA, maxA;
				FindMinMaxVertices(axis, hullA.permutedVertices, minA, maxA);

				float minB, maxB;
				FindMinMaxVertices(axis, hullB.permutedVertices, minB, maxB);

				// Test for separation
				if (minB > maxA || minA > maxB)
				{
					return false;
				}
			}
		}

		if (performEdgeTests)
		{
			// Cross each axis against each axis, generate a new axis then add them to the list
			// Axes are quantized before they are hashed so we can eliminate 'close enough' duplicates
			for (int i = 0; i < hullA.edgeAxisVectors.Count(); i++)
			{
				for (int j = 0; j < hullB.edgeAxisVectors.Count(); j++)
				{
					math::Float4 axis = hullA.edgeAxisVectors[i].Cross(hullB.edgeAxisVectors[j]);

					float minA, maxA;
					FindMinMaxVertices(axis, hullA.permutedVertices, minA, maxA);

					float minB, maxB;
					FindMinMaxVertices(axis, hullB.permutedVertices, minB, maxB);

					// Test for separation
					if (minB > maxA || minA > maxB)
					{
						return false;
					}
				}
			}
		}

		return true;
	}

	// USE_NEAR_AND_FAR_PLANES allows for an optimization, see bellow
	template<bool USE_NEAR_AND_FAR_PLANES>
	struct FrustumConvexHull
	{
		ConvexHullView GetView() const
		{
			return ConvexHullView(
				ConstMemoryView<math::Float4>(edgeAxisVectors, 4 + (int)USE_NEAR_AND_FAR_PLANES * 2),
				ConstMemoryView<math::Float4>(faceAxisVectors, 4 + (int)USE_NEAR_AND_FAR_PLANES),
				ConstMemoryView<math::Float4>(facePlanes, 4 + (int)USE_NEAR_AND_FAR_PLANES * 2),
				ConstMemoryView<math::Float4>(permutedVertices, 6)
			);
		}

		// User must provide a valid inverse projection view matrix
		void InitializeFromProjectionMatrix(const math::Matrix4& inverseViewProjectionMatrix4x4)
		{
			// Compute the 8 corners of the frustum
			math::Float4 corners[8]
			{
			   inverseViewProjectionMatrix4x4.Multiply(math::Float4(-1, -1, +0.01f, 1)),
			   inverseViewProjectionMatrix4x4.Multiply(math::Float4(+1, -1, +0.01f, 1)),
			   inverseViewProjectionMatrix4x4.Multiply(math::Float4(-1, +1, +0.01f, 1)),
			   inverseViewProjectionMatrix4x4.Multiply(math::Float4(+1, +1, +0.01f, 1)),
			   inverseViewProjectionMatrix4x4.Multiply(math::Float4(-1, -1, +1, 1)),
			   inverseViewProjectionMatrix4x4.Multiply(math::Float4(+1, -1, +1, 1)),
			   inverseViewProjectionMatrix4x4.Multiply(math::Float4(-1, +1, +1, 1)),
			   inverseViewProjectionMatrix4x4.Multiply(math::Float4(+1, +1, +1, 1)),
			};

			math::Float4::Components sc[8];

			// Perspective divide, (XYZW / W)
			for (int i = 0; i < 8; i++)
			{
				math::Float4 W = corners[i].Shuffle<3>();
				corners[i] = corners[i].Div(W);
				cornerPositions[i] = corners[i];

				sc[i] = corners[i].SplitComponents();
			}

			// Compute the edge axis vectors along the frustum Z axis
			edgeAxisVectors[0] = corners[4].Sub(corners[0]);
			edgeAxisVectors[1] = corners[5].Sub(corners[1]);
			edgeAxisVectors[2] = corners[6].Sub(corners[2]);
			edgeAxisVectors[3] = corners[7].Sub(corners[3]);

			math::Float4 nearX = corners[1].Sub(corners[0]);
			math::Float4 nearY = corners[2].Sub(corners[0]);

			// Stash the near XY vectors if we need them later
			if (USE_NEAR_AND_FAR_PLANES)
			{
				// Compute the optional edge axis vectors for the near/far plane X/Y axis
				edgeAxisVectors[4] = nearX;
				edgeAxisVectors[5] = nearY;
			}

			// The 4 normals for the 
			faceAxisVectors[0] = edgeAxisVectors[0].Cross(nearX);
			faceAxisVectors[1] = nearY.Cross(edgeAxisVectors[0]);
			faceAxisVectors[2] = nearX.Cross(edgeAxisVectors[3]);
			faceAxisVectors[3] = edgeAxisVectors[3].Cross(nearY);
			
			facePlanes[0] = faceAxisVectors[0].Div(sqrt(faceAxisVectors[0].Dot3(faceAxisVectors[0])));
			facePlanes[1] = faceAxisVectors[1].Div(sqrt(faceAxisVectors[1].Dot3(faceAxisVectors[1])));
			facePlanes[2] = faceAxisVectors[2].Div(sqrt(faceAxisVectors[2].Dot3(faceAxisVectors[2])));
			facePlanes[3] = faceAxisVectors[3].Div(sqrt(faceAxisVectors[3].Dot3(faceAxisVectors[3])));

			// Compute optional near/far plane axis
			if (USE_NEAR_AND_FAR_PLANES)
			{
				faceAxisVectors[4] = nearX.Cross(nearY);

				facePlanes[4] = faceAxisVectors[4].Div(sqrt(faceAxisVectors[4].Dot3(faceAxisVectors[4])));
				facePlanes[5] = facePlanes[4].Mul(-1.0f); // Axis of plane 4 and 5 will be identical except flipped
			}

			math::Float4::Components c;

			c = facePlanes[0].SplitComponents();
			facePlanes[0] = math::Float4(c.x, c.y, c.z, facePlanes[0].Dot3(corners[0]));

			c = facePlanes[1].SplitComponents();
			facePlanes[1] = math::Float4(c.x, c.y, c.z, facePlanes[1].Dot3(corners[0]));

			c = facePlanes[2].SplitComponents();
			facePlanes[2] = math::Float4(c.x, c.y, c.z, facePlanes[2].Dot3(corners[3]));

			c = facePlanes[3].SplitComponents();
			facePlanes[3] = math::Float4(c.x, c.y, c.z, facePlanes[3].Dot3(corners[3]));

			// Compute optional near/far plane axis
			if (USE_NEAR_AND_FAR_PLANES)
			{
				c = facePlanes[4].SplitComponents();
				facePlanes[4] = math::Float4(c.x, c.y, c.z, facePlanes[4].Dot3(corners[0]));

				c = facePlanes[5].SplitComponents();
				facePlanes[5] = math::Float4(c.x, c.y, c.z, facePlanes[5].Dot3(corners[0]));
			}

			// Set 1 (Near vertices)
			permutedVertices[0] = math::Float4(sc[0].x, sc[1].x, sc[2].x, sc[3].x); // XXXX
			permutedVertices[1] = math::Float4(sc[0].y, sc[1].y, sc[2].y, sc[3].y); // YYYY
			permutedVertices[2] = math::Float4(sc[0].z, sc[1].z, sc[2].z, sc[3].z); // ZZZZ

			// Set 2 (Far vertices)
			permutedVertices[3] = math::Float4(sc[4].x, sc[5].x, sc[6].x, sc[7].x); // XXXX
			permutedVertices[4] = math::Float4(sc[4].y, sc[5].y, sc[6].y, sc[7].y); // YYYY
			permutedVertices[5] = math::Float4(sc[4].z, sc[5].z, sc[6].z, sc[7].z); // ZZZZ
		}

		// Cached frustum corners, useful for debug drawing the frustum
		math::Float4 cornerPositions[8];

		// For a regular frustum, 4 edge axes are needed for the edges of the view along the Z axis
		// Another optional 2 are needed for the near and far plane X axis and Y axis since they will be the same
		math::Float4 edgeAxisVectors[4 + (int)USE_NEAR_AND_FAR_PLANES * 2];

		// Frustum will need up to 5 face axes, this is because near and far planes will be the same, just opposites
		// however as an optimization we can omit the near and far planes giving us a total of 4 planes
		// This will result in less accurate culling for things that are far away, 
		// but if we know that the far plane is further away than anything we can see then this is safe to do
		// For this to work correctly the matrix supplied to initialize must have a far plane further than any leaf node
		math::Float4 faceAxisVectors[4 + (int)USE_NEAR_AND_FAR_PLANES];

		// Face planes, the planes normals XYZ and distance W 
		math::Float4 facePlanes[4 + (int)USE_NEAR_AND_FAR_PLANES * 2];

		// 8 vertices, 4 per register and 3 components per vector
		// (8 / 4) * 3 == 6
		math::Float4 permutedVertices[6];
	};

	// Specialization for AABB convex hulls, we can abuse certain properties and cutdown on storage requirements significantly
	struct AABBConvexHullStorage
	{
		ConvexHullView GetView() const
		{
			// There is a total of 3 edges with unique axes, and a total of 3 face axes, and they happen to be the same
			// And since they are axis aligned they will always be the same, so we can cheat here and collapse them all into the same static memory
			return ConvexHullView(
				ConstMemoryView<math::Float4>(math::AXES_VECTORS, 3),
				ConstMemoryView<math::Float4>(math::AXES_VECTORS, 3),
				ConstMemoryView<math::Float4>(nullptr, 0), // No faces provided by AABB
				ConstMemoryView<math::Float4>(permutedVertices, 6)
			);
		}

		// 8 vertices, 4 per register and 3 components per vector
		// (8 / 4) * 3 == 6
		math::Float4 permutedVertices[6];
	};

	struct AABB
	{
		AABB() :
			min(math::Float4(0, 0, 0, 1)),
			max(math::Float4(0, 0, 0, 1))
		{}

		AABB(const math::Float4& _min, const math::Float4& _max) :
			min(_min),
			max(_max)
		{}

		AABB ComputeMerged(const AABB& other) const
		{
			AABB out;
			out.min = min.Min(other.min);
			out.max = max.Max(other.max);
			return out;
		}

		// Returns true if this fully contains containee
		bool Contains(const AABB& containee)
		{
			__m128i lt = _mm_castps_si128(_mm_cmplt_ps(containee.min.m, min.m));
			__m128i gt = _mm_castps_si128(_mm_cmpgt_ps(containee.max.m, max.m));
			__m128i one = _mm_set1_epi32(0xFFFFFFFF);

			int r0 = _mm_test_all_zeros(lt, one);
			int r1 = _mm_test_all_zeros(gt, one);

			return r0 != 0 && r1 != 0;
		}

		float ComputeVolume() const
		{
			math::Float4 extents = max.Sub(min);
			math::Float4::Components components = extents.SplitComponents();
			return components.x * components.y * components.z;
		}

		float ComputeSurfaceArea() const
		{
			math::Float4 extents = max.Sub(min);
			math::Float4 rotated = extents.Shuffle<1, 2, 0, 0>();

			math::Float4::Components c = extents.Mul(rotated).SplitComponents();
			return 2 * (c.x + c.y + c.z);
		}

		float ComputeCostHeuristic() const
		{
			return ComputeSurfaceArea();
		}

		math::Float4 ComputeCenter() const
		{
			return min.Add(max).Mul(math::HALF);
		}

		math::Float4 ComputeHalfExtents() const
		{
			return max.Sub(min).Mul(math::HALF);
		}

		// Precompute convex hull optimized storage
		AABBConvexHullStorage ComputeConvexHullStorage() const
		{
			auto smin = min.SplitComponents();
			auto smax = max.SplitComponents();

			AABBConvexHullStorage ret;

			// Set 1 (top 4 vertices)
			ret.permutedVertices[0] = math::Float4(smin.x, smax.x, smin.x, smax.x); // XXXX
			ret.permutedVertices[1] = math::Float4(smin.y, smin.y, smin.y, smin.y); // YYYY
			ret.permutedVertices[2] = math::Float4(smin.z, smin.z, smax.z, smax.z); // ZZZZ

			// Set 2 (bottom 4 vertices)
			ret.permutedVertices[3] = math::Float4(smin.x, smax.x, smin.x, smax.x); // XXXX
			ret.permutedVertices[4] = math::Float4(smax.y, smax.y, smax.y, smax.y); // YYYY
			ret.permutedVertices[5] = math::Float4(smin.z, smin.z, smax.z, smax.z); // ZZZZ

			return ret;
		}

		math::Float4 min;
		math::Float4 max;
	};

	class DrawInterface
	{
	public:
		virtual void DrawBounds(const AABB& bounds, const math::Float4& color) const = 0;
	};

	template<typename UserDataType = uintptr_t>
	class BVH3D
	{
	public:
		typedef unsigned int IndexType;

		static_assert(std::is_integral<IndexType>::value, "Index type must be an integer");
		static_assert(std::is_unsigned<IndexType>::value, "Index type must be signed");

		static const IndexType INVALID_INDEX_VALUE = std::numeric_limits<IndexType>::max();

		// The maximum number of elements this tree can support is index value divided by 1
		// for 16 bit tree that would be 32,000 elements which should be more than plenty for 99% of cases
		static const IndexType MAX_NODES = INVALID_INDEX_VALUE - 1;

		typedef void(*QueryCallbackFunc)(const UserDataType&, const AABB& bounds, void* user);

		struct QueryStats
		{
			int failedSATIntersections = 0;
			int successfulSATIntersections = 0;
			int numVIHTests = 0;
			int numHSITests = 0;
			int separatingAxisTheoremSkipped = 0;
			int intersectionsBypassed = 0;
			int totalNodesVisited = 0;

			double totalHSITime = 0;
			double totalVIHTime = 0;
			double totalSATTime = 0;

			int leafCount = 0;
		};

		struct HandleData
		{
			IndexType node;
			AABB originalBounds;
			UserDataType userData;
		};

		struct Node
		{
			IndexType parent;
			IndexType child0;
			IndexType child1;
			IndexType handleIndex;
			AABB bounds;

			// Cached precomputed state
			float cost;

			bool isLeaf() const
			{
				if (child0 == INVALID_INDEX_VALUE)
				{
					BVH_ASSERT(child1 == INVALID_INDEX_VALUE);
					return true;
				}

				return false;
			}

#if USE_VALIDATION
			bool isDead;
#endif
		};

		BVH3D(int maxExpectedEntries = 0)
		{
			// Need double the number of maximum expected entries to support the branch nodes
			mNodes.reserve(maxExpectedEntries * 2);
			mHandles.reserve(maxExpectedEntries);
		}

		~BVH3D()
		{
		}

		// Insert user data with a bounding volume into the structure, returns a handle
		Handle Insert(const AABB& bounds)
		{
			IndexType index = InsertIntoBVH(bounds);
			if (index == INVALID_INDEX_VALUE)
			{
				return INVALID_HANDLE;
			}

			Handle handle = AllocateHandle();

			// Doubly linked handle
			mHandles[handle.id].originalBounds = bounds;

			mHandles[handle.id].node = index;
			mNodes[index].handleIndex = handle.id;

			return handle;
		}

		void Update(Handle handle, const AABB& newBounds)
		{
			BVH_ASSERT(handle.id != INVALID_HANDLE.id);

			// Reallocate the internal index
			IndexType index = mHandles[handle.id].node;

			auto& node = mNodes[index];

			// No change needed
			if (node.bounds.Contains(newBounds))
			{
				// Update the bounds
				mHandles[handle.id].originalBounds = newBounds;
				return;
			}

			math::Float4 minVel = newBounds.min.Sub(mHandles[handle.id].originalBounds.min);
			math::Float4 maxVel = newBounds.max.Sub(mHandles[handle.id].originalBounds.max);

			static int cycle = 0;
			float advance = 8.0f + (cycle++ % 4);
			AABB predictedAABB = AABB(
				newBounds.min.Add(minVel.Mul(advance)),
				newBounds.max.Add(maxVel.Mul(advance)));

			AABB expandedBounds = predictedAABB.ComputeMerged(newBounds);

			// Grab the user data and store it locally
			IndexType handleIndex = node.handleIndex;

			// Remove the old index
			RemoveFromBVH(index);

			// Reinsert the user data and get the new index
			IndexType newIndex = InsertIntoBVH(expandedBounds);

			mNodes[newIndex].handleIndex = handleIndex;

			// Remap the handle to the new indexs
			mHandles[handle.id].originalBounds = newBounds;
			mHandles[handle.id].node = newIndex;
		}

		// Remove leaf entry
		void Remove(Handle handle)
		{
			// Unlink the handle from the node
			mNodes[mHandles[handle.id].node].handleIndex = INVALID_INDEX_VALUE;
			mHandles[handle.id].node = INVALID_INDEX_VALUE;

			RemoveFromBVH(mHandles[handle.id].node);
			ReleaseHandle(handle);
		}

		// Access the bounding volume of a leaf
		const AABB& GetBounds(Handle handle) const
		{
			return mNodes[handle.id].bounds;
		}

		// Access the user data of a leaf
		UserDataType& GetUserData(Handle handle)
		{
			return mNodes[handle.id].userData;
		}

		// Access the user data a leaf
		const UserDataType& GetUserData(Handle handle) const
		{
			return mNodes[handle.id].userData;
		}



		// Reorganizes the internal memory layout to optimize for query traversal
		void Optimize()
		{
			OptimizeInternal();
		}


		// Rebuild the structure of the BVH to optimize for queries
		// This may be important after loading a level or something
		void Rebuild()
		{
			std::vector<IndexType> leafNodes;
			leafNodes.reserve(mNodes.size());

			struct Range
			{
				IndexType rangeBegin;
				IndexType rangeEnd;
				IndexType node;
			};

			std::vector<Range> pendingRanges;
			pendingRanges.reserve(mNodes.size());

			// Release all non leaf nodes and find the indices for all leaf nodes
			for (IndexType i = 0; i < mNodes.size(); i++)
			{
				if (!mNodes[i].isLeaf())
				{
					ReleaseNode(i);
				}
				else
				{
					leafNodes.push_back(i);
				}
			}

			// TODO handle the case when there is only 1 node?

			auto computeRangeBounds = [&](IndexType rangeBegin, IndexType rangeEnd) -> AABB
			{
				math::Float4 min = math::LARGEST_POSITIVE;
				math::Float4 max = math::LARGEST_NEGATIVE;

				for (IndexType i = rangeBegin; i < rangeEnd; i++)
				{
					min = mNodes[leafNodes[i]].bounds.min.Min(min);
					max = mNodes[leafNodes[i]].bounds.max.Max(max);
				}

				return AABB(min, max);
			};

			Range rootRange;
			rootRange.rangeBegin = (IndexType)0;
			rootRange.rangeEnd = (IndexType)leafNodes.size();
			rootRange.node = AllocateNode(computeRangeBounds(rootRange.rangeBegin, rootRange.rangeEnd));
			mRootIndex = rootRange.node;

			// Starting range is the whole thing
			pendingRanges.push_back(rootRange);

			while (pendingRanges.size())
			{
				Range range = pendingRanges.back();
				pendingRanges.pop_back();

				math::Float4 min = mNodes[range.node].bounds.min;
				math::Float4 max = mNodes[range.node].bounds.max;

				// Measure largest axis
				math::Float4::Components extents = max.Sub(min).SplitComponents();

				int largestAxis = 0;
				if (extents.x > extents.y && extents.x > extents.z)
				{
					largestAxis = 0;
				}
				else if (extents.y > extents.x && extents.y > extents.z)
				{
					largestAxis = 1;
				}
				else
				{
					largestAxis = 2;
				}

				// Sort by that axis
				std::sort(leafNodes.begin() + range.rangeBegin, leafNodes.begin() + range.rangeEnd,
					[&](IndexType a, IndexType b)
					{
						math::Float4 ap = mNodes[a].bounds.max.Add(mNodes[a].bounds.min).Mul(0.5f);
						math::Float4 bp = mNodes[b].bounds.max.Add(mNodes[b].bounds.min).Mul(0.5f);
						return ap.SplitComponents()[largestAxis] > bp.SplitComponents()[largestAxis] ? 1 : 0;
					});

				// Split the leaf nodes in the middle and create two new ranges
				// If we only have 1 item in each range, that means its a leaf
				IndexType rangeS = range.rangeBegin;
				IndexType rangeC = range.rangeBegin + (range.rangeEnd - range.rangeBegin) / 2;
				IndexType rangeE = range.rangeEnd;

				BVH_ASSERT((rangeC - rangeS) >= 1);
				BVH_ASSERT((rangeE - rangeC) >= 1);

				// Relationship links
				IndexType parent = range.node;
				IndexType child0 = INVALID_INDEX_VALUE;
				IndexType child1 = INVALID_INDEX_VALUE;

				if ((rangeC - rangeS) > 1)
				{
					Range newRange;
					newRange.rangeBegin = rangeS;
					newRange.rangeEnd = rangeC;
					newRange.node = AllocateNode(computeRangeBounds(newRange.rangeBegin, newRange.rangeEnd));
					child0 = newRange.node;

					pendingRanges.push_back(newRange);
				}
				else
				{
					child0 = leafNodes[rangeS];
				}

				if ((rangeE - rangeC) > 1)
				{
					Range newRange;
					newRange.rangeBegin = rangeC;
					newRange.rangeEnd = rangeE;
					newRange.node = AllocateNode(computeRangeBounds(newRange.rangeBegin, newRange.rangeEnd));
					child1 = newRange.node;

					pendingRanges.push_back(newRange);
				}
				else
				{
					child1 = leafNodes[rangeC];
				}

				BVH_ASSERT(parent != INVALID_INDEX_VALUE);
				BVH_ASSERT(child0 != INVALID_INDEX_VALUE);
				BVH_ASSERT(child1 != INVALID_INDEX_VALUE);

				// Establish the links between parent and child
				mNodes[parent].child0 = child0;
				mNodes[child0].parent = parent;

				mNodes[parent].child1 = child1;
				mNodes[child1].parent = parent;
			}
		}

		void Query(const ConvexHullView& hull, QueryCallbackFunc callback, void* user, QueryStats* stats, DrawInterface* drawInterface) const
		{
			if (stats)
			{
				*stats = {};
			}

			ConvexHullQueryNode(hull, mRootIndex, callback, user, stats, drawInterface, false);
		}

		void Draw(const DrawInterface* callback)
		{
			DrawTree(callback, mRootIndex, 0);
		}

		float ComputeTotalCost()
		{
			float totalCost = 0;

			for (const auto& node : mNodes)
			{
				if (node.parent != INVALID_INDEX_VALUE &&
					node.child0 != INVALID_INDEX_VALUE &&
					node.child1 != INVALID_INDEX_VALUE)
				{
					totalCost += node.bounds.ComputeSurfaceArea();
				}
			}

			return totalCost;
		}

	private:
		Handle AllocateHandle()
		{
			Handle handle = INVALID_HANDLE;
			if (mFreeHandles.size() > 0)
			{
				handle = mFreeHandles[mFreeHandles.size() - 1];
				mFreeHandles.pop_back();
			}
			else
			{
				handle.id = NumericCast<int>(mHandles.size());
				mHandles.resize(mHandles.size() + 1);
			}

			return handle;
		}

		void ReleaseHandle(Handle handle)
		{
			mFreeHandles.push_back(handle);
		}

		void PrecomputeNodeData(Node& node)
		{
			node.cost = node.bounds.ComputeCostHeuristic();

			math::Float4 center = node.bounds.ComputeCenter();





		}

		IndexType AllocateNode(const AABB& bounds)
		{
			IndexType index = INVALID_INDEX_VALUE;

			if (mFreeNodes.size() > 0)
			{
				index = mFreeNodes[mFreeNodes.size() - 1];
				mFreeNodes.pop_back();
			}
			else
			{
				if (mNodes.size() == MAX_NODES / 2)
				{
					BVH_FAIL("Reached half the maximum node count, \
this is a warning that the internal type might need to be changed to higher bit count. \
The BVH will continue to operate normally until the limit is reached");
				}
				else if (mNodes.size() >= MAX_NODES)
				{
					BVH_FAIL("Maximum node count reached, the BVH behaviour is now undefined");
					return INVALID_INDEX_VALUE;
				}

				index = NumericCast<IndexType>(mNodes.size());
				mNodes.resize(mNodes.size() + 1);
			}

			Node& node = mNodes[index];
			node.child0 = INVALID_INDEX_VALUE;
			node.child1 = INVALID_INDEX_VALUE;
			node.parent = INVALID_INDEX_VALUE;
#if USE_VALIDATION
			node.isDead = false;
#endif
			node.bounds = bounds;
			PrecomputeNodeData(node);

			return index;
		}

		void ReleaseNode(IndexType index)
		{
			Node& node = mNodes[index];
			node.parent = INVALID_INDEX_VALUE;
			node.child0 = INVALID_INDEX_VALUE;
			node.child1 = INVALID_INDEX_VALUE;
#if USE_VALIDATION
			node.isDead = true;
#endif

			mFreeNodes.push_back(index);
		}

		float ComputeCostDelta(const AABB& bounds, IndexType nodeIndex, float directCost) const
		{
			float costBefore = mNodes[nodeIndex].bounds.ComputeCostHeuristic();
			return directCost - costBefore;
		}

		float ComputeDirectCost(const AABB& bounds, IndexType nodeIndex) const
		{
			return mNodes[nodeIndex].bounds.ComputeMerged(bounds).ComputeCostHeuristic();
		}

		IndexType FindBestSibling(const AABB& bounds) const
		{
			float bestCost = FLT_MAX;
			IndexType bestNode = INVALID_INDEX_VALUE;

			float boundsCost = bounds.ComputeCostHeuristic();

			if (mRootIndex == INVALID_INDEX_VALUE)
			{
				return INVALID_INDEX_VALUE;
			}

			struct Helper
			{
				IndexType node;
				float inheritedCost;
			};

			static const int MAX_HELPER_STACK_SIZE = 128;
			Helper helperStack[MAX_HELPER_STACK_SIZE];
			int helperCount = 0;

			auto push = [&](const Helper& helper)
			{
				helperStack[helperCount] = helper;
				helperCount += 1;
				BVH_ASSERT(helperCount < MAX_HELPER_STACK_SIZE);
			};

			auto pop = [&]()
			{
				helperCount -= 1;
				return helperStack[helperCount];
			};

			push({ mRootIndex, 0.0f });

			int numIterations = 0;
			while (helperCount != 0)
			{
				numIterations += 1;
				const auto& entry = pop();

				IndexType nodeIndex = entry.node;
				float inheritedCost = entry.inheritedCost;
				float directCost = ComputeDirectCost(bounds, nodeIndex);

				float cost = inheritedCost + directCost;

				if (cost < bestCost)
				{
					bestNode = nodeIndex;
					bestCost = cost;
				}

				if (!mNodes[nodeIndex].isLeaf())
				{
					float costDelta = (directCost - mNodes[nodeIndex].cost);
					float newCost = inheritedCost + costDelta;

					float lowerBoundChildCost = newCost + boundsCost;
					if (lowerBoundChildCost < bestCost)
					{
						push({ mNodes[nodeIndex].child0, newCost });
						push({ mNodes[nodeIndex].child1, newCost });
					}
				}
			}

			volatile static int maxQueueSize = 0;
			if (numIterations > maxQueueSize)
			{
				maxQueueSize = numIterations;

				char buffer[64];
				snprintf(buffer, sizeof(buffer), "%i\n", maxQueueSize);
				OutputDebugStringA(buffer);
			}

			return bestNode;
		}

		void Rotate(IndexType index)
		{
			IndexType parentIndex = mNodes[index].parent;
			if (parentIndex == INVALID_INDEX_VALUE)
			{
				return;
			}

			IndexType grandParentIndex = mNodes[parentIndex].parent;
			if (grandParentIndex == INVALID_INDEX_VALUE)
			{
				return;
			}

			IndexType uncleIndex = INVALID_INDEX_VALUE;
			if (mNodes[grandParentIndex].child0 == parentIndex)
			{
				uncleIndex = mNodes[grandParentIndex].child1;
			}
			else
			{
				uncleIndex = mNodes[grandParentIndex].child0;
			}
			BVH_ASSERT(uncleIndex != INVALID_INDEX_VALUE);

			Validate();

			// Measure the cost of rotating or not

			IndexType siblingIndex = INVALID_INDEX_VALUE;
			if (mNodes[parentIndex].child0 == index)
			{
				siblingIndex = mNodes[parentIndex].child1;
			}
			else
			{
				siblingIndex = mNodes[parentIndex].child0;
			}
			BVH_ASSERT(siblingIndex != INVALID_INDEX_VALUE);

			AABB tmpBounds;
			tmpBounds = mNodes[siblingIndex].bounds.ComputeMerged(mNodes[index].bounds);
			float costA = tmpBounds.ComputeCostHeuristic();

			tmpBounds = mNodes[siblingIndex].bounds.ComputeMerged(mNodes[uncleIndex].bounds);
			float costB = tmpBounds.ComputeCostHeuristic();

			if (costA > costB)
			{
				// Swap the uncle with the child
				if (mNodes[grandParentIndex].child0 == uncleIndex)
				{
					mNodes[grandParentIndex].child0 = index;
				}
				else
				{
					mNodes[grandParentIndex].child1 = index;
				}

				mNodes[index].parent = grandParentIndex;

				if (mNodes[parentIndex].child0 == index)
				{
					mNodes[parentIndex].child0 = uncleIndex;
				}
				else
				{
					mNodes[parentIndex].child1 = uncleIndex;
				}

				mNodes[uncleIndex].parent = parentIndex;

				// Refit the original parent
				mNodes[parentIndex].bounds = tmpBounds;
				PrecomputeNodeData(mNodes[parentIndex]);

				Validate();
			}
		}

		IndexType InsertIntoBVH(const AABB& bounds)
		{
			IndexType siblingIndex = INVALID_INDEX_VALUE;
			if (mNodes.size() > 0)
			{
				siblingIndex = FindBestSibling(bounds);
			}

			if (siblingIndex == INVALID_INDEX_VALUE)
			{
				// Must be the first node, add it
				mRootIndex = AllocateNode(bounds);
				if (mRootIndex == INVALID_INDEX_VALUE)
				{
					return INVALID_INDEX_VALUE;
				}

				Validate();
				return mRootIndex;
			}
			else
			{
				AABB newBounds = bounds.ComputeMerged(mNodes[siblingIndex].bounds);

				IndexType newParent = AllocateNode(newBounds);
				if (newParent == INVALID_INDEX_VALUE)
				{
					return INVALID_INDEX_VALUE;
				}

				IndexType leafIndex = AllocateNode(bounds);
				if (leafIndex == INVALID_INDEX_VALUE)
				{
					ReleaseNode(newParent);
					return INVALID_INDEX_VALUE;
				}

				IndexType oldParent = mNodes[siblingIndex].parent;
				mNodes[newParent].parent = oldParent;


				if (oldParent != INVALID_INDEX_VALUE)
				{
					if (mNodes[oldParent].child0 == siblingIndex)
					{
						mNodes[oldParent].child0 = newParent;
					}
					else
					{
						mNodes[oldParent].child1 = newParent;
					}

					mNodes[newParent].child0 = siblingIndex;
					mNodes[newParent].child1 = leafIndex;
					mNodes[siblingIndex].parent = newParent;
					mNodes[leafIndex].parent = newParent;
				}
				else
				{
					mNodes[newParent].child0 = siblingIndex;
					mNodes[newParent].child1 = leafIndex;
					mNodes[siblingIndex].parent = newParent;
					mNodes[leafIndex].parent = newParent;
					mRootIndex = newParent;
				}

				IndexType tmp = mNodes[leafIndex].parent;
				while (tmp != INVALID_INDEX_VALUE)
				{
					IndexType c0 = mNodes[tmp].child0;
					IndexType c1 = mNodes[tmp].child1;

					AABB combinedChildBounds = mNodes[c0].bounds.ComputeMerged(mNodes[c1].bounds);
					mNodes[tmp].bounds = combinedChildBounds;
					PrecomputeNodeData(mNodes[tmp]);

					Rotate(tmp);
					Validate();

					tmp = mNodes[tmp].parent;
				}

				Validate();
				return leafIndex;
			}
		}

		void RemoveFromBVH(IndexType index)
		{
			if (index == mRootIndex)
			{
				mRootIndex = INVALID_INDEX_VALUE;
				ReleaseNode(index);
				Validate();
				return;
			}

			IndexType parent = mNodes[index].parent;
			BVH_ASSERT(parent != INVALID_INDEX_VALUE);

			IndexType sibling = INVALID_INDEX_VALUE;

			if (mNodes[parent].child0 == index)
			{
				sibling = mNodes[parent].child1;
			}
			else
			{
				sibling = mNodes[parent].child0;
			}

			BVH_ASSERT(sibling != INVALID_INDEX_VALUE);

			if (parent == mRootIndex)
			{
				mRootIndex = sibling;
				ReleaseNode(parent);
				ReleaseNode(index);

				// As this is now the root, there is no parent
				mNodes[sibling].parent = INVALID_INDEX_VALUE;

				Validate();
				return;
			}

			// Fuse the grandparent and sibling kicking the current index and parent out of the tree allowing us to free them

			IndexType grandParent = mNodes[parent].parent;
			BVH_ASSERT(grandParent != INVALID_INDEX_VALUE);

			if (mNodes[grandParent].child0 == parent)
			{
				mNodes[grandParent].child0 = sibling;
			}
			else
			{
				mNodes[grandParent].child1 = sibling;
			}

			mNodes[sibling].parent = grandParent;

			ReleaseNode(parent);
			ReleaseNode(index);

			// Refit the tree
			IndexType tmp = grandParent;
			while (tmp != INVALID_INDEX_VALUE)
			{
				IndexType c0 = mNodes[tmp].child0;
				IndexType c1 = mNodes[tmp].child1;

				AABB combinedChildBounds = mNodes[c0].bounds.ComputeMerged(mNodes[c1].bounds);
				mNodes[tmp].bounds = combinedChildBounds;
				PrecomputeNodeData(mNodes[tmp]);

				tmp = mNodes[tmp].parent;
			}

			Validate();
		}

		void DrawTree(const DrawInterface* callbacks, IndexType currentNode, int depth) const
		{
			if (currentNode == INVALID_INDEX_VALUE)
				return;

			bool isLeaf = mNodes[currentNode].isLeaf();

			if (isLeaf)
			{
				callbacks->DrawBounds(mNodes[currentNode].bounds, math::Float4(1, 0, 0, 1));
			}


			DrawTree(callbacks, mNodes[currentNode].child0, depth + 1);
			DrawTree(callbacks, mNodes[currentNode].child1, depth + 1);
		}

		void ConvexHullQueryNode(const ConvexHullView& hull, IndexType currentNode, QueryCallbackFunc callback, void* user, QueryStats* stats, DrawInterface* drawInterface, bool fullyInside) const
		{
			if (mRootIndex == INVALID_INDEX_VALUE)
				return;

			Validate();

			const Node& node = mNodes[currentNode];

			AABBConvexHullStorage permutedAABB = node.bounds.ComputeConvexHullStorage();

			ConvexHullView aabbView = permutedAABB.GetView();

			if (stats)
			{
				stats->totalNodesVisited++;
			}

			bool intersects = false;

			// Perform a very broad approximation of our intersection test by computing a bounding sphere and testing it against our hull
			// This is the simplest form of intersection we can do, it allows us to reject future tests if its fully outside or fully inside
			// And then we only need to perform the high accuracy tests in the event that it intersects with the hull
			if (!intersects)
			{
				if (stats)
				{
					stats->numHSITests++;
				}

				math::Float4 center = node.bounds.ComputeCenter();
				math::Float4 halfExtents = node.bounds.ComputeHalfExtents();
				float radius = sqrt(halfExtents.Dot3(halfExtents));

				Timer timer;
				int r = HullSphereIntersection(hull, center, radius);
				if (stats)
				{
					stats->totalHSITime += timer.getMicroseconds();
				}

				if (r > 0)
				{
					// Fully outside, If the sphere is fully outside, then we know no other test will pass
					return;
				}
				else if (r < 0)
				{
					// Fully inside, If the sphere is fully inside, can bypass all future checks
					intersects = true;
					fullyInside = true;
				}
			}

			if (!intersects)
			{
				if (stats)
				{
					stats->numVIHTests++;
				}

				Timer timer;

				int r = CountVerticesInsideHull(hull, aabbView.permutedVertices);

				if (stats)
				{
					stats->totalVIHTime += timer.getMicroseconds();
				}

				if (r > 0)
				{
					// All vertices are outside, 
					// this could be a false positive so we need to continue and do a SAT test
				}
				else if (r < 0)
				{
					if (stats)
					{
						// By calling true here, it means we skipped 1 SAT test
						stats->separatingAxisTheoremSkipped++;
					}

					// Fully inside, If the sphere is fully inside, can bypass all future checks
					intersects = true;
					fullyInside = true;
				}
				else
				{
					// Some vertices are inside, which means SAT is going to pass so we can skip it
					// However we are not fully inside so we have to continue checking the children
					intersects = true;
				}
			}
			else
			{
				// If we reached here, it means we skipped all intersection tests and we are just bypassing through
				if (stats)
				{
					stats->intersectionsBypassed++;
				}
			}

			if (!intersects)
			{
				Timer timer;

				intersects = IntersectsWith(hull, aabbView, false);

				if (stats)
				{
					stats->totalSATTime += timer.getMicroseconds();
				}

				if (stats)
				{
					if (intersects)
					{
						stats->successfulSATIntersections++;
					}
					else
					{
						stats->failedSATIntersections++;
					}
				}
			}

			if (intersects)
			{
				if (node.isLeaf())
				{
					if (stats)
					{
						stats->leafCount++;
					}

					if (drawInterface)
					{
						drawInterface->DrawBounds(mHandles[node.handleIndex].originalBounds, math::Float4(0, 1, 0, 1));
						drawInterface->DrawBounds(node.bounds, math::Float4(0, 1, 0, 0.5f));
					}

					// Leaf
					callback(mHandles[node.handleIndex].userData, mHandles[node.handleIndex].originalBounds, user);

				}
				else
				{
					if (drawInterface)
					{
						if (fullyInside)
						{
							drawInterface->DrawBounds(node.bounds, math::Float4(0.2f, 0.2f, 0.9f, 0.9f));
						}
						else
						{
							drawInterface->DrawBounds(node.bounds, math::Float4(0.3f, 0.3f, 0.3f, 0.1f));
						}
					}

					BVH_ASSERT(node.child1 != INVALID_INDEX_VALUE);
					ConvexHullQueryNode(hull, node.child0, callback, user, stats, drawInterface, fullyInside);
					ConvexHullQueryNode(hull, node.child1, callback, user, stats, drawInterface, fullyInside);
				}
			}
		}

		struct OptimizeHelper
		{
			IndexType order;
			IndexType index;
		};

		void OptimizeGenerateOrder(std::vector<OptimizeHelper>& data, IndexType nodeIndex, int& counter)
		{
			data[nodeIndex].order = counter;
			counter += 1;

			data[nodeIndex].index = nodeIndex;

			if (mNodes[nodeIndex].child0 != INVALID_INDEX_VALUE)
			{
				OptimizeGenerateOrder(data, mNodes[nodeIndex].child0, counter);
			}

			if (mNodes[nodeIndex].child1 != INVALID_INDEX_VALUE)
			{
				OptimizeGenerateOrder(data, mNodes[nodeIndex].child1, counter);
			}
		}

		void OptimizeInternal()
		{
			std::vector<OptimizeHelper> data;
			data.resize(mNodes.size());

			std::vector<int> oldToNewLUT;
			oldToNewLUT.resize(mNodes.size());

			std::vector<Node> newNodes;
			newNodes.resize(mNodes.size());

			int counter = 0;
			OptimizeGenerateOrder(data, mRootIndex, counter);

			std::sort(data.begin(), data.end(), [](const OptimizeHelper& a, const OptimizeHelper& b)
				{
					return a.order < b.order;
				});

			for (int i = 0; i < data.size(); i++)
			{
				oldToNewLUT[data[i].index] = i;
			}

			mRootIndex = oldToNewLUT[mRootIndex];

			for (int i = 0; i < data.size(); i++)
			{
				IndexType oldIndex = data[i].index;

				newNodes[i] = mNodes[oldIndex];

				if (newNodes[i].parent != INVALID_INDEX_VALUE)
				{
					newNodes[i].parent = oldToNewLUT[newNodes[i].parent];
				}

				if (newNodes[i].child0 != INVALID_INDEX_VALUE)
				{
					newNodes[i].child0 = oldToNewLUT[newNodes[i].child0];
				}

				if (newNodes[i].child1 != INVALID_INDEX_VALUE)
				{
					newNodes[i].child1 = oldToNewLUT[newNodes[i].child1];
				}
			}

			for (int i = 0; i < mHandles.size(); i++)
			{
				mHandles[i].node = oldToNewLUT[i];
			}

			mNodes = std::move(newNodes);
		}

		void Validate() const
		{
#if USE_VALIDATION
			if (!mValidationEnabled)
			{
				return;
			}

			for (size_t i = 0; i < mNodes.size(); i++)
			{
				auto& iter = mNodes[i];

				if (iter.isDead)
					continue;

				if (iter.child0 == INVALID_INDEX_VALUE)
				{
					ASSERTF(iter.child1 == INVALID_INDEX_VALUE, "No node can only have 1 child");
				}

				if (iter.child1 == INVALID_INDEX_VALUE)
				{
					ASSERTF(iter.child0 == INVALID_INDEX_VALUE, "No node can only have 1 child");
				}

				if (iter.child0 != INVALID_INDEX_VALUE && iter.child0 >= NumericCast<IndexType>(mNodes.size()))
				{
					BVH_FAIL("Index out of range");
				}

				if (iter.child1 != INVALID_INDEX_VALUE && iter.child1 >= NumericCast<IndexType>(mNodes.size()))
				{
					BVH_FAIL("Index out of range");
				}

				if (iter.parent != INVALID_INDEX_VALUE && iter.parent >= NumericCast<IndexType>(mNodes.size()))
				{
					BVH_FAIL("Index out of range");
				}

				if (iter.child0 != INVALID_INDEX_VALUE)
				{
					ASSERTF(mNodes[iter.child0].parent == i, "Mismatch child/parent relationship");
					ASSERTF(mNodes[iter.child0].isDead == false, "Pointing to dead node");
				}

				if (iter.child1 != INVALID_INDEX_VALUE)
				{
					ASSERTF(mNodes[iter.child1].parent == i, "Mismatch child/parent relationship");
					ASSERTF(mNodes[iter.child1].isDead == false, "Pointing to dead node");
				}

				if (iter.parent != INVALID_INDEX_VALUE)
				{
					ASSERTF(mNodes[iter.parent].child0 == i || mNodes[iter.parent].child1 == i, "Mismatch child/parent relationship");
					ASSERTF(mNodes[iter.parent].isDead == false, "Pointing to dead node");
				}

				if (mRootIndex != i && iter.parent == INVALID_INDEX_VALUE)
				{
					BVH_FAIL("Can!have a non root node without a parent");
				}

				if (iter.child1 == INVALID_INDEX_VALUE)
				{
					//ASSERTF(iter.userData != 0, "Expecting valid user data on leaf node");
				}
			}
#endif
		}

	private:

		IndexType mRootIndex = INVALID_INDEX_VALUE;

		// Validate the entire tree every time something happens
		bool mValidationEnabled = false;

		// Node structure is very small and light weight allowing us to order it so 
		// parent/child relationships between common nodes fit into a single cacheline and are very quick to traverse
		std::vector<Node> mNodes;

		// Indirect mappings from handle to index and index to handle
		// This indirection allows remapping the internal indices to optimize during runtime
		std::vector<HandleData> mHandles;

		// Free lists to handle quick allocation
		std::vector<IndexType> mFreeNodes;
		std::vector<Handle> mFreeHandles;

	};


}
