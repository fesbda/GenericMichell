#ifndef MESH_H
#define MESH_H

#include "tools.h"

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/quaternion.hpp>

#include <limits>

struct Transformation
{

private:
    glm::dvec3 translation {};
    glm::dvec3 rot_origin {};
    glm::dvec3 euler_angles {};
    glm::dmat3x3 rot_mat {};

public:
    Transformation() { set_rotation({0, 0, 0}); }

    void set_translation(glm::dvec3 pos) { translation = pos; }

    void set_centre(glm::dvec3 coord) { rot_origin = coord; }

    void set_rotation(glm::dvec3 angles)
    {
        euler_angles = angles;
        rot_mat = glm::toMat3(glm::quat(angles));
    }

    glm::dvec3 get_translation() const { return translation; }

    glm::dvec3 get_rotation() const { return euler_angles; }

    glm::dvec3 get(glm::dvec3 local_coord) const { return translation + rot_origin + rot_mat * (local_coord - rot_origin); }
};


class Mesh
{

private:
    std::vector<glm::dvec3> vertices;
    std::vector<glm::ivec3> face_vertices;
    std::vector<glm::dvec3> face_normals;
    std::vector<double> face_areas;
    std::vector<glm::dvec3> vertex_normals;

    glm::dvec3 bmin, bmax;

public:
    Mesh() { }

    Mesh(std::vector<glm::dvec3> const& verts, std::vector<glm::ivec3> const& face_vers)
        : vertices(verts)
        , face_vertices(face_vers)
    {
        update();
    }

    ~Mesh() {}

    void set_vertices(std::vector<glm::dvec3> const& verts) { this->vertices = verts; }

    void set_faces(std::vector<glm::ivec3> const& face_verts) { this->face_vertices = face_verts; }

    int add_vertex(double x, double y, double z)
    {
        vertices.push_back({x, y, z});
        return int(vertices.size()) - 1;
    }

    int add_face(int v1, int v2, int v3)
    {
        face_vertices.push_back({v1, v2, v3});
        return int(face_vertices.size()) - 1;
    }

    void update()
    {
        face_normals.resize(face_vertices.size());
        face_areas.resize(face_vertices.size());

        // calc faces info: normals, areas
        #pragma omp parallel for
        for (std::size_t fi = 0; fi < face_vertices.size(); fi++) {

            auto nodes = face_vertices[fi];
            auto v1 = vertices[nodes.x];
            auto v2 = vertices[nodes.y];
            auto v3 = vertices[nodes.z];
            auto crs = glm::cross(v2 - v1, v3 - v1);
            double const crs_len = glm::length(crs);
            face_areas[fi] = 0.5 * crs_len;
            // Guard degenerate (zero-area) faces: glm::normalize(0) is NaN, which would poison
            // every downstream area-weighted force. A zero-area face carries zero weight anyway.
            // (glm::normalize kept for the non-degenerate branch so the result is bit-identical.)
            face_normals[fi] = crs_len > 1e-300 ? glm::normalize(crs) : glm::dvec3(0.0);
        }


        // bounding box
        constexpr double dmax = std::numeric_limits<double>::max();
        constexpr double dlow = std::numeric_limits<double>::lowest();
        bmin = {dmax, dmax, dmax};
        bmax = {dlow, dlow, dlow};

        for (std::size_t vi = 0; vi < vertices.size(); vi++) {

            auto v = vertices[vi];

            bmin.x = std::min(bmin.x, v.x);
            bmin.y = std::min(bmin.y, v.y);
            bmin.z = std::min(bmin.z, v.z);

            bmax.x = std::max(bmax.x, v.x);
            bmax.y = std::max(bmax.y, v.y);
            bmax.z = std::max(bmax.z, v.z);
        }
    }

    void reverse_normals()
    {
        for (auto& n : face_normals) { n = -n; }
    }

    const std::vector<glm::ivec3>& get_faces() const { return face_vertices; }

    const std::vector<glm::dvec3>& get_vertices() const { return vertices; }

    const std::vector<glm::dvec3>& get_face_normals() const { return face_normals; }

    const std::vector<double>& get_face_areas() const { return face_areas; }

    glm::dvec3 get_min() const { return bmin; }

    glm::dvec3 get_max() const { return bmax; }

    glm::dvec3 get_size() const { return bmax - bmin; }

    glm::dvec3 get_vertex(std::size_t id, Transformation const& t) const { return t.get(vertices[id]); }

};


#endif // MESH_H
