/*2
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

#include <nori/common.h>
#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f& sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f& sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToTent(const Point2f& sample) {
    throw NoriException("Warp::squareToTent() is not yet implemented!");
}

float Warp::squareToTentPdf(const Point2f& p) {
    throw NoriException("Warp::squareToTentPdf() is not yet implemented!");
}

Point2f Warp::squareToUniformDisk(const Point2f& sample) {
    //throw NoriException("Warp::squareToUniformDisk() is not yet implemented!");
    float r = sqrtf(sample.x());
    float phi = 2.0f * M_PI * sample.y();

    return Point2f{ r * cos(phi), r * sin(phi) };
}

float Warp::squareToUniformDiskPdf(const Point2f& p) {
    //throw NoriException("Warp::squareToUniformDiskPdf() is not yet implemented!");
    float x_sq = p.x() * p.x();
    float y_sq = p.y() * p.y();
    float len_sq = x_sq + y_sq;

    return len_sq > 1.0f ? 0.0f : INV_PI;
}

Vector3f Warp::squareToUniformSphere(const Point2f& sample) {
    throw NoriException("Warp::squareToUniformSphere() is not yet implemented!");
}

float Warp::squareToUniformSpherePdf(const Vector3f& v) {
    throw NoriException("Warp::squareToUniformSpherePdf() is not yet implemented!");
}

Vector3f Warp::squareToUniformHemisphere(const Point2f& sample) {
    throw NoriException("Warp::squareToUniformHemisphere() is not yet implemented!");
}

float Warp::squareToUniformHemispherePdf(const Vector3f& v) {
    throw NoriException("Warp::squareToUniformHemispherePdf() is not yet implemented!");
}

// Sample a cosine-weighted vector on the unit hemisphere with respect to solid angles.
Vector3f Warp::squareToCosineHemisphere(const Point2f& sample) {
    //throw NoriException("Warp::squareToCosineHemisphere() is not yet implemented!");
    float phi = 2.0f * M_PI * sample.x();
    float r = sqrtf(sample.y());

    float x = r * cosf(phi);
    float y = r * sinf(phi);
    float z = sqrtf(1.0f - sample.y());

    return Vector3f{ x, y, z };
}

// Density of squareToCosineHemisphere() with respect to solid angles. 
float Warp::squareToCosineHemispherePdf(const Vector3f& v) {
    //throw NoriException("Warp::squareToCosineHemispherePdf() is not yet implemented!");
    return fmax(0.0f, v.z() * INV_PI);
}

Vector3f Warp::squareToBeckmann(const Point2f& sample, float alpha) {
    //throw NoriException("Warp::squareToBeckmann() is not yet implemented!");
    float const theta = 2.0f * M_PI * sample.x();

    float tph_sq = fmax(0.0f, -alpha * alpha * logf(1.0f - sample.y()));
    float cph = 1.0f / sqrtf(1.0f + tph_sq);
    float sph = sqrtf(1.0f - cph * cph);
    float cth = cosf(theta);
    float sth = sinf(theta);

    return Vector3f{ sph * cth, sph * sth, cph };
}

float Warp::squareToBeckmannPdf(const Vector3f& m, float alpha) {
    //throw NoriException("Warp::squareToBeckmannPdf() is not yet implemented!");
    float result = 0.0f;

    float cth = m.z();

    if (cth > 0.0f) {
        float cth_sq = cth * cth;
        float tan_theta_sq = (1.0f - cth_sq) / cth_sq;

        float alpha_sq = alpha * alpha;

        // Note that scalar of 2 has been canceled out of the operands
        float const longitudinal = expf(-tan_theta_sq / alpha_sq) / (alpha_sq * cth_sq * cth);
        float const azimuthal = INV_PI;

        result = longitudinal * azimuthal;
    }

    return result;
}


NORI_NAMESPACE_END
