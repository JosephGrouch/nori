#include <nori/volume.h>
#include <nori/common.h>
#include <nori/emitter.h>
#include "nori/block.h"
#include "nori/integrator.h"
#include "nori/mesh.h"
#include "nori/scene.h"
#include "nori/warp.h"

NORI_NAMESPACE_BEGIN

	enum PlaneType {
		UV = 0,
		VT,
		UT,
		UAlphaT,
	};

	enum SinglePlaneStrategy {
		SPS_UV,
		SPS_VT,
		SPS_UT,
		SPS_Average,
		SPS_DiscreteMIS,
		SPS_UAlpha,
		SPS_ContinousMIS,
		SPS_SMISAll,
		SPS_SMISJacobian
	};

	// Converts a mesh to a rectangular emitter to be used as a photon plane
	struct RectangularLightSource {
		Point3f o;
		Vector3f n;
		Vector3f u;
		Vector3f v;
		float u_l;
		float v_l;
		Color3f m_radiance;

		RectangularLightSource(Mesh* mesh) {
			if (mesh->getVertexCount() != 3 && mesh->getIndices().size() != 2)
				throw NoriException("Emitter must be rectangular.");

			MatrixXf verts = mesh->getVertexPositions();

			// TODO CHECK THIS!!!
			o = verts.col(0);
			u = verts.col(1) - o;
			v = verts.col(2) - o;
			u_l = u.norm();
			v_l = v.norm();
			u.normalize();
			v.normalize();
			n = u.cross(v);

			m_radiance = mesh->getEmitter()->getEmission();
		}
	};

	struct PhotonPlaneIts {
		float t_cam;
		float t0;
		float t1;
		float inv_det;
	};

	struct SinglePhotonPlane {
		// Plane origin
		Point3f o;
		// First edge (normalized)
		Vector3f d0;
		// Second edge (normalized)
		Vector3f d1;
		// First edge length
		float length0;
		// Second edge length
		float length1;
		// Random number used to generate the plane (TODO: Unused?)
		Point2f sample;
		// How the plane have been generated
		PlaneType plane_type;
		// This factor will vary between the different light sources
		Color3f weight;
		size_t id_emitter;
		float sample_alpha;

		Point3f light_position(const RectangularLightSource& light, const PhotonPlaneIts& plane_its) const {
			if (plane_type == UV)
				return light.o + light.u * plane_its.t0 + light.v * plane_its.t1;

			return o + d0 * plane_its.t0;
		}

		Color3f contrib(Vector3f d) {
			float jacobian = fabs(d1.cross(d0).dot(d));
			return weight / jacobian;
		}

		static Point2f planeIntersect(const RectangularLightSource& light, const Vector2f& d, const Point2f& o) {
			Vector2f t_0 = {-o.x() / d.x(), -o.y() / d.y()};
			Vector2f t_1 = {light.u_l - t_0.x(), light.v_l - t_0.y()};
			Vector2f t_max_coord = {fmax(t_0.x(), t_1.x()), fmax(t_0.y(), t_1.y())};
			return o + d * fmin(t_max_coord.x(), t_max_coord.y());
		}

		bool rayIntersect(const Ray3f& r, PhotonPlaneIts& plane_its) {
			Vector3f e0 = d0 * length0;
			Vector3f e1 = d1 * length1;

			Vector3f p = r.d.cross(e1);
			float det = e0.dot(p);
			if (fabs(det) < Epsilon) {
				return false;
			}

			Vector3f t = r.o - o;
			float t0 = t.dot(p) / det;
			if (t0 < 0.0f || t0 > 1.0f) {
				return false;
			}

			Vector3f q = t.cross(e0);
			float t1 = r.d.dot(q) / det;
			if (t0 < 0.0f || t0 > 1.0f) {
				return false;
			}

			float t_cam = e1.dot(q) / det;
			if (t_cam <= r.mint || t_cam >= r.maxt) {
				return false;
			}

			// Must scale to the correct distance for correct transmittance sampling
			t1 = t1 * length1;
			t0 = t0 * length0;

			plane_its.t_cam = t_cam;
			plane_its.t0 = t0;
			plane_its.t1 = t1;
			plane_its.inv_det = 1.0f / det;

			return true;
		}

		SinglePhotonPlane(PlaneType t, RectangularLightSource& light, Vector3f& d, Point2f& s2, float s_alpha,
		                  float t_sampled, size_t id, const Color3f& sigma_s) {
			sample_alpha = s_alpha;
			id_emitter = id;
			plane_type = t;
			sample = s2;

			switch (t) {
				case UV:
					o = light.o + d * t_sampled;
					d0 = light.u;
					d1 = light.v;
					length0 = light.u_l;
					length1 = light.v_l;
					weight = M_PI * light.m_radiance / sigma_s;
					break;

				case VT:
					o = light.o + light.u * light.u_l * s2.x();
					d0 = light.v;
					d1 = d;
					length0 = light.u_l;
					length1 = t_sampled;
					weight = M_PI * light.v_l * light.m_radiance;
					break;

				case UT:
					o = light.o + light.v * light.v_l * sample.y();
					d0 = light.u;
					d1 = d;
					length0 = light.u_l;
					length1 = t_sampled;
					weight = M_PI * light.v_l * light.m_radiance;
					break;

				case UAlphaT:
					Point2f o_plane = {sample.x() * light.u_l, sample.y() * light.v_l};
					float alpha = M_PI * sample_alpha;
					Point2f d_plane = {cosf(alpha), sinf(alpha)};
					auto p1_2d = planeIntersect(light, d_plane, o_plane);
					auto p2_2d = planeIntersect(light, -d_plane, o_plane);

					Point3f p1 = light.o + p1_2d.x() * light.u + p1_2d.y() * light.v;
					Point3f p2 = light.o + p2_2d.x() * light.u + p2_2d.y() * light.v;

					Point3f u_plane = p2 - p1;
					float u_plane_length = u_plane.norm();
					u_plane.normalize();

					o = p1;
					d0 = u_plane;
					d1 = d;
					length0 = u_plane_length;
					length1 = t_sampled;
					weight = M_PI * light.m_radiance * (light.u_l * light.v_l);
					break;
			}
		}
	};


	struct IntegratorSinglePlane : Integrator {
		int emitter_count = 0;
		int number_plane_gen = 0;
		int n_samples = 0;
		size_t nb_primitive;
		SinglePlaneStrategy strategy;
		bool stratified;
		std::vector<SinglePhotonPlane> planes;
		std::vector<RectangularLightSource> rect_lights;

		static SinglePhotonPlane generate_plane(PlaneType t, std::vector<RectangularLightSource>& light, size_t id,
		                                        Sampler* sampler, const HomogeneousVolume* m) {
			Vector3f d_out = Vector3f::Zero();

			while (d_out.z() == 0.0f)
				d_out = Warp::squareToCosineHemisphere(sampler->next2D());

			// TODO Check that this actually works or default to alternative implementation
			Frame frame(light[id].n);
			auto d = frame.toWorld(d_out);
			Ray3f ray_med{light[id].o, d};
			auto mrec = m->sample(ray_med, sampler->next2D());

			auto sample = sampler->next2D();
			return SinglePhotonPlane{t, light[id], d, sample, sampler->next1D(), mrec.continued_t, id, m->sigma_s};
		}

		void preprocess(const Scene* scene) override {
			Sampler* sampler = static_cast<Sampler*>(NoriObjectFactory::createInstance("independent", PropertyList()));

			for (Mesh* mesh : scene->getMeshes()) {
				if (mesh->isEmitter()) {
					++emitter_count;
					rect_lights.emplace_back(mesh);
				}
			}

			const HomogeneousVolume* m = scene->getVolume();

			while (planes.size() < nb_primitive) {
				size_t id_emitter = sampler->next1D() * rect_lights.size();

				switch (strategy) {
					case SPS_UT:
						planes.emplace_back(generate_plane(UT, rect_lights, id_emitter, sampler, m));
						break;
					case SPS_VT:
						planes.emplace_back(generate_plane(VT, rect_lights, id_emitter, sampler, m));
						break;
					case SPS_UV:
						planes.emplace_back(generate_plane(UV, rect_lights, id_emitter, sampler, m));
						break;
					case SPS_DiscreteMIS:
						planes.emplace_back(generate_plane(UT, rect_lights, id_emitter, sampler, m));
						planes.emplace_back(generate_plane(VT, rect_lights, id_emitter, sampler, m));
						planes.emplace_back(generate_plane(UV, rect_lights, id_emitter, sampler, m));
						break;
					case SPS_UAlpha:
					case SPS_ContinousMIS:
					case SPS_SMISAll:
					case SPS_SMISJacobian:
						planes.emplace_back(generate_plane(UAlphaT, rect_lights, id_emitter, sampler, m));
				}

				++number_plane_gen;
			}

			// Build the BVH to speedup the computation...
			//BVHAccel bvh_plane = BHVAccel::create(planes);
		}

		Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const override {
			Sampler* sampler_ecmis = static_cast<Sampler*>(NoriObjectFactory::createInstance("independent",
				PropertyList()));
			const HomogeneousVolume* m = scene->getVolume();

			Color3f c = Color3f::Zero();
			//for (plane_its, b_id) in bvh_plane.gather(ray) {
			//for (auto test in bvh_plane.gather(ray)) {
			//    PhotonPlaneIts plane_its;
			//    SinglePhotonPlane plane = bvh_plane.elements[b_id];
			// This code is if we do not use BVH
			for (auto plane : planes) {
				PhotonPlaneIts plane_its;
				if (!plane.rayIntersect(ray, plane_its))
					continue;

				Point3f p_hit = ray.o + ray.d * plane_its.t_cam;
				Point3f p_light = plane.light_position(rect_lights[plane.id_emitter], plane_its);
				RectangularLightSource rect_light = rect_lights[plane.id_emitter];
				//if (accel.visible(p_hit, p_light)) {
				Ray3f ray_tr{ray.o, ray.d};
				ray_tr.maxt = plane_its.t_cam;
				Color3f transmittance = m->transmittance(ray_tr);
				Color3f rho = Color3f(0.25f / M_PI);
				float w = 0.0f;
				switch (strategy) {
					case SPS_UT:
					case SPS_UV:
					case SPS_VT:
					case SPS_UAlpha:
					case SPS_ContinousMIS:
					case SPS_SMISAll:
					case SPS_SMISJacobian:
						w = 1.0f;
						break;
					case SPS_Average:
						w = 1.0f / 3.0f;
						break;
					case SPS_DiscreteMIS:
						// Need to compute all possible shapes
						Vector3f d = p_hit - p_light;
						// TODO: Not used
						float t_sampled = d.norm();
						d.normalize();
						auto t_planes = std::vector<SinglePhotonPlane>();
						t_planes.emplace_back(UV, rect_light, d, plane.sample, 0.0f, t_sampled, plane.id_emitter,
						                      m->sigma_s);
						t_planes.emplace_back(UT, rect_light, d, plane.sample, 0.0f, t_sampled, plane.id_emitter,
						                      m->sigma_s);
						t_planes.emplace_back(VT, rect_light, d, plane.sample, 0.0f, t_sampled, plane.id_emitter,
						                      m->sigma_s);
						// FIXME: Normally this code is unecessary
						// 	As we can reuse the plane retrived.
						// 	However, it seems to have a miss match between photon planes
						//	contribution calculation.
						int debug_id = plane.plane_type;
						// NOT IMPLEMENTED
						assert(debug_id != UAlphaT);

						// TODO THIS IS A HACKY TRANSLATION, so must check for bugs
						float avg = t_planes[debug_id].contrib(ray.d).sum() / 3.0f;
						//float inv = 1.0f / avg;
						float sum = 0.0f;
						for (auto p : t_planes) {
							float pc = p.contrib(ray.d).sum() / 3.0f;
							if (pc != 0.0f && isfinite(pc))
								sum += 1.0f / pc;
						}
						//w = inv / sum
						if (avg != 0.0f && sum != 0.0f)
							w = 1.0f / avg / sum;
						break;
				}

				// Compute the plane weighted contribution
				Color3f contrib;
				switch (strategy) {
						// Deng et al. CMIS
					case SPS_ContinousMIS:
						// Here we use their integration from
						// Normally, all the jacobian simplifies
						// So it is why we need to have a special estimator
						contrib = plane.weight / 2.0f / INV_PI / sqrtf(powf(rect_light.u.cross(plane.d1).dot(ray.d),
						                                                    2.0f) + powf(rect_light.v.cross(plane.d1).
							                                               dot(ray.d), 2.0f));

						// SMIS
						//SPS_SMISAll(n_samples)
						//SPS_SMISJacobian(n_samples) = > {
					case SPS_SMISAll:
					case SPS_SMISJacobian: assert(n_samples > 0);

						// Compute wrap random number for generating fake planes that generate same
						// path configuration. Indeed, we need to be sure that the new planes alpha plane
						// cross the same point on the light source
						Point3f p_l = p_light - rect_light.o;
						// Mitigate floating point precision issue
						Point2f sample_wrap = {
							clamp(p_l.dot(rect_light.u) / rect_light.u_l, 0.0f, 1.0f),
							clamp(p_l.dot(rect_light.v) / rect_light.v_l, 0.0f, 1.0f)
						};
						// --------------------
						// Create N-1 fake planes with the same path configuration
						// using stratified sampling.

						// We will compute the inverse norm of the SMIS weight
						// We initialize the variable with the actual plane that
						// the ray intersected
						float inv_norm;
						switch (strategy) {
							case SPS_SMISAll:
								// J(..) * p(e)
								// we do not include A (light source area)
								// as it cancel out
								inv_norm = fabs(plane.d1.cross(plane.d0).dot(ray.d)) * plane.length0;
								break;
							case SPS_SMISJacobian:
								// J(..)
								inv_norm = fabs(plane.d1.cross(plane.d0).dot(ray.d));
								break;
							default:
								throw NoriException("Invalid sampling strategy.");
						}

						// Generate the other planes
						float offset = plane.sample_alpha;
						for (int i = 0; i < n_samples - 1; ++i) {
							// The new alpha
							float new_alpha;
							if (stratified) {
								float rat = (static_cast<float>(i) + 1.0f) / static_cast<float>(n_samples);
								new_alpha = offset + (rat - truncf(rat));
							}
							else {
								new_alpha = sampler_ecmis->next1D();
							}
							assert(new_alpha >= 0.0f && new_alpha <= 1.0f);

							// This construct a new plane (call previous snipped)
							// to generate the plane
							SinglePhotonPlane new_plane(UAlphaT, rect_light, plane.d1, sample_wrap, new_alpha, 0.0f,
							                            plane.id_emitter, m->sigma_s);

							// Accumulate the weight
							switch (strategy) {
								case SPS_SMISAll:
									inv_norm += fabs(new_plane.d1.cross(new_plane.d0).dot(ray.d)) * new_plane.length0;
									break;
								case SPS_SMISJacobian:
									inv_norm += fabs(new_plane.d1.cross(new_plane.d0).dot(ray.d));
									break;
								default:
									throw NoriException("Invalid sampling strategy.");
							}
						}

						// The SMIS weight and the plane contribution
						// First each weight have J(..) * p(e) or J(..) at the numerator
						// however this factor cancel out with the plane evaluation: f(..) / (J(..) * p(e))
						// With that, all planes contribute to the same :)
						w = 1.0f / inv_norm;
						switch (strategy) {
							case SPS_SMISAll:
								// Note we need to also cancel out lenght0
								// as it will cancel out inside the SMIS_weight
								contrib = plane.weight * plane.length0 * static_cast<float>(n_samples);
								break;
							case SPS_SMISJacobian:
								contrib = plane.weight * static_cast<float>(n_samples);
								break;
							default:
								throw NoriException("Invalid sampling strategy.");

								contrib *= w;
						}

						// Default: evaluate and weight the contrib
						// plane.contrib(..) =  plane.weight / jacobian
						// w: other MIS compute above
						contrib = w * plane.contrib(ray.d); // Do nothing and just evaluate the plane
				}

				// Compute the rest of the term
				// and accumulate them
				c += rho * transmittance * m->sigma_s * contrib * static_cast<float>(emitter_count) * (1.0f /
					static_cast<float>(number_plane_gen));
				//}

				return c;
			}
		}
	};

NORI_NAMESPACE_END
