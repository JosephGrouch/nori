#include "nori/integrator_utils.h"
#include "nori/emitter.h"
#include "nori/bsdf.h"

NORI_NAMESPACE_BEGIN

	Color3f IntegratorUtils::processDirectIllumination(Scene const* scene, const Vector3f& ray_d,
	                                                      const Intersection& its, const Point2f& sample) {
		Color3f result = Color3f::Zero();
		// Select a random emitting mesh from those in the scene
		float emitter_pd;
		Mesh* emitter = scene->getRandomEmitterMesh(emitter_pd);

		// Sample a position on the emitter mesh
		Vector3f emitter_p;
		Normal3f emitter_n;
		// Probability here will be 1 / (total surface area of emitter)
		const float surface_pd = emitter->sampleSurface(sample, emitter_p, emitter_n);

		// Sampled emitter point to initial intersection components
		const Vector3f its_to_emitter = emitter_p - its.p;
		const Normal3f wo = its_to_emitter.normalized();
		const Normal3f wi = -wo;
		
		Intersection its_rev;

		// Apply the ray tracing operator r(p,-wc) using shadow ray from initial ray
		// intersection towards sampled emitter
		if (scene->rayIntersect(Ray3f(its.p, wo), its_rev) && its_rev.mesh == emitter) {
			Color3f le = emitter->getEmitter()->eval(emitter_p, emitter_n, wi);

			// Calculate the geometric term using the following formula:
			// G(x<->y) := V(x<->y) (|nx*(x->y)|*|ny*(y->x)|) / (||x-y||^2)
			//				 			(a)			(b)				(c)
			const float a = fabs(its.shFrame.n.dot(wo));
			const float b = fabs(emitter_n.dot(wi));
			const float c = its_to_emitter.dot(its_to_emitter);

			const float geo_term = a * b / c;

			// Evaluate the Fr term
			BSDFQueryRecord bsdf_rec{its.toLocal(wo), its.toLocal(-ray_d), ESolidAngle};
			const Color3f fr = its.mesh->getBSDF()->eval(bsdf_rec);

			// Remember that the value must be scaled by probability density
			result = fr * geo_term * le / emitter_pd / surface_pd;
		}

		return result;
	}

	inline float IntegratorUtils::balanceHeuristic(const int nf, const float f_pdf, const int ng, const float g_pdf) {
		const float f = static_cast<float>(nf) * f_pdf;
		const float g = static_cast<float>(ng) * g_pdf;
		return f / (f + g);
	}

	inline float IntegratorUtils::powerHeuristic(const int nf, const float f_pdf, const int ng, const float g_pdf) {
		const float f = static_cast<float>(nf) * f_pdf;
		const float f_sq = f * f;
		const float g = static_cast<float>(ng) * g_pdf;
		const float g_sq = g * g;
		return f_sq / (f_sq + g_sq);
	}

NORI_NAMESPACE_END
