#include "nori/emitter.h"
#include "nori/mesh.h"

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter {
public:
	AreaLight(const PropertyList& props) {
		// Obtain radiance from scene properties
		m_radiance = props.getColor("radiance");
	}

	Color3f eval(const Point3f& p, const Normal3f& n, const Vector3f& d) const override {
		// Return radiance if it is positive in direction d, otherwise return zero
		return d.dot(n) > 0.0f
			? m_radiance
			: Color3f::Zero();
	}

	Color3f getEmission() const override {
		return m_radiance;
	}

	std::string toString() const override {
		return "AreaLight[\n"
			"	radiance: " + m_radiance.toString() + "\n"
			"]";
	}

private:
	Color3f m_radiance = Color3f::Zero();
};

NORI_REGISTER_CLASS(AreaLight, "area");

NORI_NAMESPACE_END
