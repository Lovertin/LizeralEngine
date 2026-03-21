#pragma once
#include <vector>
#include "runtime/core/math/vector2.h"
#include "runtime/core/math/vector3.h"
// #include "runtime/function/res_type/mesh/mesh.h"
#include "runtime/function/res_type/Model/Mesh.h"

using vec2 = Lizeral::Vector2;
using vec3 = Lizeral::Vector3;

namespace Lizeral{

    class ImportModel{
    private:
        int numVertices;
        std::vector<vec3> vertices;
        std::vector<vec2> texCoords;
        std::vector<vec3> normalVecs;
    public:
        ImportModel();
        ImportModel(const char *filePath);
        int getNumVertices();
        std::vector<vec3> getVertices();
        std::vector<vec2> getTextureCoords();
        std::vector<vec3> getNormals();

        std::vector<Vertex> GetMeshVertices();
    };

    class ModelImporter
    {
    private:
        std::vector<float> vertVals;
        std::vector<float> triangleVerts;
        std::vector<float> textureCoords;
        std::vector<float> stVals;
        std::vector<float> normals;
        std::vector<float> normVals;
    public:
        ModelImporter();
        void parseOBJ(const char *filePath);
        int getNumVertices();
        std::vector<float> getVertices();
        std::vector<float> getTextureCoordinates();
        std::vector<float> getNormals();
    };

}