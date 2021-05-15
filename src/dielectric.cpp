/*
	This file is part of Nori, a simple educational ray tracer

	Copyright (c) 2015 by Wenzel Jakob

	Nori is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License Version 3
	as published by the Free Software Foundation.

	Nori is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/bsdf.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
public:
	Dielectric(const PropertyList& propList) {
		/* Interior IOR (default: BK7 borosilicate optical glass) */
		m_intIOR = propList.getFloat("intIOR", 1.5046f);

		/* Exterior IOR (default: air) */
		m_extIOR = propList.getFloat("extIOR", 1.000277f);
	}

	Color3f eval(const BSDFQueryRecord&) const {
		/* Discrete BRDFs always evaluate to zero in Nori */
		return Color3f(0.0f);
	}

	float pdf(const BSDFQueryRecord&) const {
		/* Discrete BRDFs always evaluate to zero in Nori */
		return 0.0f;
	}

	Color3f sample(BSDFQueryRecord& bRec, const Point2f& sample) const override {
		//throw NoriException("Unimplemented!");
		bRec.measure = EDiscrete;

		// Using the fresnel utility function here calculates values needed later, so inlining
		// may be optimization
		float fr = fresnel(bRec.wi.z(), m_extIOR, m_intIOR);

		// Use sample to determine if this will be reflection or refraction
		if (sample.x() < fr) {
			//
			// Direct case
			// 
			// This results in the same value as the more complicated multiplication
			bRec.wo.x() = -bRec.wi.x();
			bRec.wo.y() = -bRec.wi.y();
			bRec.wo.z() = bRec.wi.z();
			bRec.eta = 1.0f;
		}
		else {
			//
			// Refracting case
			//
			bool is_entering = bRec.wi.z() <= 0.0f;
			float eta_i = m_extIOR;
			float eta_t = m_intIOR;

			if (is_entering) {
				std::swap(eta_i, eta_t);
				bRec.wi.z() = -bRec.wi.z();
			}

			// Compute refraction using Snell's law
			const float eta = eta_i / eta_t;
			bRec.eta = eta;

			float sth_t_sq = eta * eta * Frame::sinTheta2(bRec.wi);
			float cth_t = std::sqrt(1.0f - sth_t_sq);

			bRec.wo.x() = eta * -bRec.wi.x();
			bRec.wo.y() = eta * -bRec.wi.y();
			// Outgoing z-component is just negative cosTheta
			bRec.wo.z() = -cth_t;

			// Flip the reflection z component back if opposite side
			if (is_entering)
				bRec.wo.z() = -bRec.wo.z();
		}

		return Color3f::Ones();
	}

	std::string toString() const {
		return tfm::format(
			"Dielectric[\n"
			"  intIOR = %f,\n"
			"  extIOR = %f\n"
			"]",
			m_intIOR, m_extIOR);
	}

private:
	// Refractive index of the interior
	float m_intIOR;
	// Refractive index of the side that contains the surface normal
	float m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");

NORI_NAMESPACE_END
