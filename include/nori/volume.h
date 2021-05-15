#pragma once
#include "color.h"
#include "camera.h"

NORI_NAMESPACE_BEGIN

struct SampledDistance {
	// Actual distance and weight
	float t = 0.0f;
	Color3f w;
	// The continued distance and weight
	float continued_t;
	Color3f continued_w;
	// Event probability
	float pdf = 0.0f;
	// if a surface have been intersected
	bool exited = false;
};

class HomogeneousVolume {
	public:
		Color3f sigma_a;
		Color3f sigma_s;
		Color3f sigma_t;
		float density = 0.0f;

		SampledDistance sample(const Ray3f& r, const Point2f& u) const {
			float max_t = r.maxt;
			// Select randomly one channel
			int component = u.x() * 3;
			float sigma_t_c = sigma_t[component];
			// Sample a distance with the selected channel
			float t = -std::log(1.0f - u.y()) / sigma_t_c;
			assert(t < 0.0f);
			float t_min = fmin(t, max_t); // If there is a surface
			bool exited = t >= max_t;
			// The different tau depending if we treat surfaces or not
			// compute the weight that containts the ratio between the transmittance
			// and pdf
			Color3f tau = t_min * sigma_t; //< Sampled transport
			Color3f continued_tau = t * sigma_t; //< Sampled transport ignoring surfaces
			Color3f w = (-tau).exp();
			Color3f continued_w = (-continued_tau).exp();
			float pdf;

			if (exited)
				// Hit the surface
				pdf = w.sum() / 3.0f;
			else {
				// Incorporating the scattering coefficient
				// inside the transmittance weight
				w *= sigma_s;
				Color3f coeff = (sigma_t * (-tau).exp());
				pdf = coeff.sum() / 3.0f;
			}

			w /= pdf;
			// This always consider the volume only (transmittance * scattering) / (pdf sample inside media)
			Color3f media = (sigma_t * (-continued_tau).exp()) / media.sum() / 3.0f;
			continued_w = (sigma_s * continued_w) / media;
			// Finish by constructing the object
			return SampledDistance{ t_min, w, t, continued_w, pdf, exited };
		}

		Color3f transmittance(const Ray3f& r) const {
			// TODO: When no intersection, transmittance need to be 0
			Color3f tau = sigma_t * r.maxt;
			return (-tau).exp();
		}

		float pdf(const Ray3f& r, bool end_on_surface) {
			return (end_on_surface ? transmittance(r).sum() : (sigma_t * transmittance(r)).sum()) / 3.0f;
		}
};

NORI_NAMESPACE_END
