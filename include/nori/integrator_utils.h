#pragma once

#include "scene.h"

NORI_NAMESPACE_BEGIN

	class IntegratorUtils {
		public:
			static Color3f processDirectIllumination(Scene const* scene, const Vector3f& ray_d,
			                                            const Intersection& its, const Point2f& sample);

			static float balanceHeuristic(int nf, float f_pdf, int ng, float g_pdf);

			static float powerHeuristic(int nf, float f_pdf, int ng, float g_pdf);
	};

NORI_NAMESPACE_END
