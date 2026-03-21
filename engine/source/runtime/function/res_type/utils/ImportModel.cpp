#include <fstream>
#include <sstream>
#include <iostream>
#include "ImportModel.h"

namespace Lizeral{

    ImportModel::ImportModel(){}

    ImportModel::ImportModel(const char *filePath) {
        ModelImporter modelImporter = ModelImporter();
        modelImporter.parseOBJ(filePath);
        numVertices = modelImporter.getNumVertices();
        std::vector<float> verts = modelImporter.getVertices();
        std::vector<float> tcs = modelImporter.getTextureCoordinates();
        std::vector<float> normals = modelImporter.getNormals();

        for (int i = 0; i < numVertices; i++) {
            vertices.push_back(vec3(verts[i*3], verts[i*3+1], verts[i*3+2]));
            texCoords.push_back(vec2(tcs[i*2], tcs[i*2+1]));
            normalVecs.push_back(vec3(normals[i*3], normals[i*3+1], normals[i*3+2]));
        }
    }

    int ImportModel::getNumVertices() { return numVertices; }
    std::vector<vec3> ImportModel::getVertices() { return vertices; }
    std::vector<vec2> ImportModel::getTextureCoords() { return texCoords; }
    std::vector<vec3> ImportModel::getNormals() { return normalVecs; }

    std::vector<Vertex> ImportModel::GetMeshVertices() {
        std::vector<Vertex> meshVertices;
        meshVertices.reserve(numVertices); 

        for (int i = 0; i < numVertices; i++) {
            Vertex v;
            v.Position = vertices[i];
            
            if (i < texCoords.size()) {
                v.TexCoords = texCoords[i];
            } else {
                v.TexCoords = vec2(0.0f, 0.0f);
            }

            if (i < normalVecs.size()) {
                v.Normal = normalVecs[i];
            } else {
                v.Normal = vec3(0.0f, 1.0f, 0.0f); 
            }

            meshVertices.push_back(v);
        }

        return meshVertices;
    }


    ModelImporter::ModelImporter() {}

    void ModelImporter::parseOBJ(const char *filePath) {
        float x, y, z;
        std::string content;
        std::ifstream fileStream(filePath, std::ios::in);
        if (!fileStream.is_open()) {
            std::cerr << "ERROR: Failed to open model file: " << filePath << std::endl;
            return; 
        }
        std::string line = "";
        while (std::getline(fileStream, line)) {
            // getline(fileStream, line);
            if (line.compare(0, 2, "v ") == 0) {
                std::stringstream ss(line.erase(0, 1));
                ss >> x; ss >> y; ss >> z;
                vertVals.push_back(x);
                vertVals.push_back(y);
                vertVals.push_back(z);
            }
            if (line.compare(0, 2, "vt") == 0) {
                std::stringstream ss(line.erase(0, 2));
                ss >> x; ss >> y;
                stVals.push_back(x);
                stVals.push_back(y);
            }
            if (line.compare(0, 2, "vn") == 0) {
                std::stringstream ss(line.erase(0, 2));
                ss >> x; ss >> y; ss >> z;
                normVals.push_back(x);
                normVals.push_back(y);
                normVals.push_back(z);
            }
            if (line.compare(0, 2, "f ") == 0) {
                std::string oneCorner, v, t, n;
                std::stringstream ss(line.erase(0, 2));
                for (int i = 0; i < 3; i++) {
                    getline(ss, oneCorner, ' ');
                    std::stringstream oneCornerSS(oneCorner);
                    getline(oneCornerSS, v, '/');
                    getline(oneCornerSS, t, '/');
                    getline(oneCornerSS, n, '/');

                    int vertRef = (stoi(v) - 1) * 3;
                    int tcRef = (stoi(t) - 1) * 2;
                    int normRef = (stoi(n) - 1) * 3;

                    triangleVerts.push_back(vertVals[vertRef]);
                    triangleVerts.push_back(vertVals[vertRef + 1]);
                    triangleVerts.push_back(vertVals[vertRef + 2]);

                    textureCoords.push_back(stVals[tcRef]);
                    textureCoords.push_back(stVals[tcRef + 1]);

                    normals.push_back(normVals[normRef]);
                    normals.push_back(normVals[normRef + 1]);
                    normals.push_back(normVals[normRef + 2]);
                }
            }
        }
    }

    int ModelImporter::getNumVertices() { return (triangleVerts.size()/3); }
    std::vector<float> ModelImporter::getVertices() { return triangleVerts; }
    std::vector<float> ModelImporter::getTextureCoordinates() { return textureCoords; }
    std::vector<float> ModelImporter::getNormals() { return normals; }
}


