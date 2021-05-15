#include "nori/bsdf.h"
#include "nori/emitter.h"
#include "nori/integrator.h"
#include "nori/integrator_utils.h"
#include "nori/sampler.h"
#include "nori/scene.h"
#include "nori/warp.h"

NORI_NAMESPACE_BEGIN

class WhittedIntegrator : public Integrator {
public:
	WhittedIntegrator(const PropertyList& props) {
	}

	Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const override {
		/* Find the surface that is visible in the requested direction */
		Intersection its;

		// Resulting radiance value
		Color3f result = Color3f::Zero();

		if (scene->rayIntersect(ray, its)) {
			const BSDF* bsdf = its.mesh->getBSDF();

			if (bsdf->isDiffuse()) {
				/*
				 * Diffuse case
				 */
				 // Get emitted radiance from the intersected mesh if there is any
				if (its.mesh->isEmitter())
					// Radiance is emission from intersected mesh
					result = its.mesh->getEmitter()->eval(its.p, its.shFrame.n, -ray.d);

				// Add the returned radiance from scatter processing
				result += IntegratorUtils::processDirectIllumination(scene, ray.d, its, sampler->next2D());
			}
			else if (sampler->next1D() < 0.95f) {
				/*
				 * Specular case
				 */
				 // Reverse the ray direction to generate the reflection
				BSDFQueryRecord bsdf_rec{ its.shFrame.toLocal(-ray.d) };

				// Sample the BSDF and generate new direction
				Color3f c = bsdf->sample(bsdf_rec, sampler->next2D());
				Vector3f wr = its.toWorld(bsdf_rec.wo);

				// Recurse through the subsequent specular reflections and return the value
				// scaled by the weight and probability
				result = c * Li(scene, sampler, Ray3f{ its.p, wr }) / 0.95f;
			}
		}

		return result;
	}

	std::string toString() const override {
		return "WhittedIntegrator[]";
	}
};

NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");

NORI_NAMESPACE_END
