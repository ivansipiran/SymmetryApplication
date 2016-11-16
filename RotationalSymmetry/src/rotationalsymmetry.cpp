#include "RotationalSymmetry/rotationalsymmetry.h"

namespace Engine{


RotationalSymmetry::RotationalSymmetry(Harris3D::Mesh* in, Util::PropertySet* p)
{
    //We perform a deep copy of the input mesh into a Harris3D::mesh
    input = in;
    prop = p;
}

Harris3D::Mesh* RotationalSymmetry::applyOpenVDB(){

        std::vector<openvdb::Vec3s> points;
        std::vector<openvdb::Vec3I> triangles;

        float diag = input->get_diagonal();
        float factor = atof(prop->getProperty("factor-volumetric").c_str());

        openvdb::math::Transform::Ptr xform = openvdb::math::Transform::createLinearTransform(factor * diag);

        for (int i = 0; i < input->get_number_vertices(); i++){
            openvdb::Vec3s xyz(0.0, 0.0, 0.0);
            xyz = xform->worldToIndex(openvdb::Vec3s(input->get_vertices()[i].x(),input->get_vertices()[i].y(),input->get_vertices()[i].z()));
            points.push_back(xyz);
        }
        for (int i = 0; i < input->get_number_triangles(); i++){
            openvdb::Vec3I triangle(input->get_triangles()[i].get_vertex_at(0),input->get_triangles()[i].get_vertex_at(1),input->get_triangles()[i].get_vertex_at(2));
            triangles.push_back(triangle);
        }

        std::vector<openvdb::Vec3s> newPoints;
        std::vector<openvdb::Vec3I> newTriangles;
        std::vector<openvdb::Vec4I> newQuads;
        openvdb::FloatGrid::Ptr levelset = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(*xform, points, triangles);
        openvdb::tools::volumeToMesh(*levelset,newPoints,newTriangles,newQuads);

        Harris3D::Mesh* newMesh = new Harris3D::Mesh();
        //input->cleanMesh();
        newMesh->set_number_vertices(newPoints.size());
        newMesh->set_number_triangles(newQuads.size() * 2);

        for (unsigned i = 0; i < newPoints.size(); i++){
            newMesh->add_vertex(i, factor*diag*newPoints[i].x(), factor*diag*newPoints[i].y(), factor*diag*newPoints[i].z());

        }
        for (unsigned i = 0; i < newQuads.size(); i++){
            newMesh->add_triangle(2*i,newQuads[i].z(),newQuads[i].y(),newQuads[i].x());
            newMesh->add_triangle(2*i+1,newQuads[i].w(),newQuads[i].z(),newQuads[i].x());
        }

        newMesh->compute_normals();
        newMesh->get_diagonal();
        SimpleMesh::IO::OFFWriter w("ja.off");
        w.write_mesh(*newMesh);
        return newMesh;
}

Harris3D::Mesh* RotationalSymmetry::decimate(Harris3D::Mesh* in){
    int numFaces = atoi(prop->getProperty("num-triangles-decimation").c_str());
    Harris3D::Mesh* out = new Harris3D::Mesh();

    if(in->get_number_triangles()>numFaces){
            for(int i=0;i<in->get_number_vertices();i++){
                Simplify::Vertex v;
                v.p.x=in->get_vertices()[i].x();
                v.p.y=in->get_vertices()[i].y();
                v.p.z=in->get_vertices()[i].z();
                Simplify::vertices.push_back(v);
            }

            for(int i=0;i<in->get_number_triangles();i++){
                Simplify::Triangle t;
                t.v[0] = in->get_triangles()[i].get_vertex_at(0);
                t.v[1] = in->get_triangles()[i].get_vertex_at(1);
                t.v[2] = in->get_triangles()[i].get_vertex_at(2);
                Simplify::triangles.push_back(t);
            }

            //cout << "Input: " << Simplify::triangles.size() << " triangles " << Simplify::vertices.size() << " vertices\n" << endl;

            Simplify::simplify_mesh(numFaces);

            //cout << "Output: " << Simplify::triangles.size() << " triangles " << Simplify::vertices.size() << " vertices\n" << endl;

            //Reconstruct mesh
            //in->cleanMesh();
            out->set_number_vertices(Simplify::vertices.size());

            int numFaces = 0;
            for(unsigned int i = 0; i < Simplify::triangles.size(); i++){
                if(!Simplify::triangles[i].deleted)
                    numFaces++;
            }
            out->set_number_triangles(numFaces);

            for (unsigned i = 0; i < Simplify::vertices.size(); i++){
                out->add_vertex(i, Simplify::vertices[i].p.x, Simplify::vertices[i].p.y, Simplify::vertices[i].p.z);
            }

            int cont = 0;
            for (unsigned i = 0; i < Simplify::triangles.size(); i++){
                if(!Simplify::triangles[i].deleted){
                    out->add_triangle(cont,Simplify::triangles[i].v[0],Simplify::triangles[i].v[1],Simplify::triangles[i].v[2]);
                    cont++;
                }
            }
            Simplify::vertices.clear();
            Simplify::triangles.clear();
            out->compute_normals();
            out->get_diagonal();
        }else{
            out->deepCopy(in);
        }

    return out;
}

Harris3D::Mesh* RotationalSymmetry::transferClassification(Harris3D::Mesh* in, Harris3D::Mesh *original){
    //Crear un Kd-tree para original
        SimpleMesh::Vertex* vertices = original->get_vertices();
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
        cloud->width = original->get_number_vertices();
        cloud->height = 1;
        cloud->points.resize(cloud->width * cloud->height);

        for(size_t i = 0; i < cloud->points.size(); i++){
            cloud->points[i].x = vertices[i].x();
            cloud->points[i].y = vertices[i].y();
            cloud->points[i].z = vertices[i].z();
        }

        pcl::KdTreeFLANN<pcl::PointXYZ> tree;
        tree.setInputCloud(cloud);

        for(int i = 0; i < in->get_number_vertices(); i++){
            pcl::PointXYZ sp;

            sp.x = in->get_vertices()[i].x();
            sp.y = in->get_vertices()[i].y();
            sp.z = in->get_vertices()[i].z();

            int K = 1;
            vector<int> idx(K);
            vector<float> dist(K);

            if(tree.nearestKSearch(sp, K, idx, dist)>0){
                in->get_vertices()[i].set_color(vertices[idx[0]].red(),vertices[idx[0]].green(), vertices[idx[0]].blue());
            }
        }

        //Eliminar los vértices (y triangulos asociados) que estan etiquetados con azul
        vector<bool> markVertices;
        markVertices.resize(in->get_number_vertices(), false);

        vector<bool> markTriangles;
        markTriangles.resize(in->get_number_triangles(), false);

        vector<bool> markBoundary;
        markBoundary.resize(in->get_number_vertices(), false);

        int numVertices = 0;
        for(int i = 0; i < in->get_number_vertices(); i++){
            if(in->get_vertices()[i].blue()>0.95){
                markVertices[i] = true;
                numVertices++;
            }
        }

        int numTriangles = 0;
        for(int i = 0; i < in->get_number_triangles(); i++){
            vector<unsigned int> p = in->get_triangles()[i].get_vertices();
            if(markVertices[p[0]] || markVertices[p[1]] || markVertices[p[2]]){
                markTriangles[i] = true;
                markBoundary[p[0]] = true;
                markBoundary[p[1]] = true;
                markBoundary[p[2]] = true;
                numTriangles++;
            }
        }

        numVertices = in->get_number_vertices() - numVertices;
        numTriangles = in->get_number_triangles() - numTriangles;

        Harris3D::Mesh* newMesh = new Harris3D::Mesh();
        newMesh->set_number_vertices(numVertices);
        newMesh->set_number_triangles(numTriangles);

        vector<unsigned int> mapping;

        unsigned int numRemoved = 0;
        for(unsigned int i = 0; i < in->get_number_vertices(); i++){
            if(!markVertices[i]){
                newMesh->add_vertex(i - numRemoved, in->get_vertices()[i].x(), in->get_vertices()[i].y(), in->get_vertices()[i].z());
                newMesh->add_normal(i - numRemoved, in->get_vertices()[i].nx(), in->get_vertices()[i].ny(), in->get_vertices()[i].nz());
                if(markBoundary[i])
                    newMesh->add_color(i - numRemoved, 1.0, 0.0, 0.0);
                else
                    newMesh->add_color(i - numRemoved, 0.0, 0.0, 0.0);

            }else{
                numRemoved++;
            }
            mapping.push_back(i - numRemoved);
        }

        int contTriangles = 0;
        for(unsigned int i = 0; i < in->get_number_triangles(); i++){
            if(!markTriangles[i]){
                vector<unsigned int> p = in->get_triangles()[i].get_vertices();
                newMesh->add_triangle(contTriangles++, mapping[p[0]], mapping[p[1]], mapping[p[2]]);
            }
        }

        newMesh->get_diagonal();

       return newMesh;
}

void RotationalSymmetry::applyRANSAC(Harris3D::Mesh* in, vector<float>& center, vector<float>& directionVector){
        float thresholdRANSAC = atof(prop->getProperty("threshold-ransac").c_str());
        vector<unsigned int> interestPoints;
        Harris3D::HarrisDetector hd(in, prop);
        hd.detectInterestPoints(interestPoints);

        std::cout << "Interest points:" << interestPoints.size() << std::endl;
        vector<unsigned int> filteredPoints;
        if(in->has_color()){
            for(unsigned int i = 0; i < interestPoints.size(); i++){
                if(in->get_vertices()[interestPoints[i]].red() < 0.95)
                    filteredPoints.push_back(interestPoints[i]);
            }
        }else{
            filteredPoints.assign(interestPoints.begin(), interestPoints.end());
        }

        //Point cloud analysis
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);

        cloud->width    = filteredPoints.size();
        cloud->height   = 1;
        cloud->is_dense = false;
        cloud->points.resize(cloud->width * cloud->height);

        for(unsigned i = 0; i < cloud->points.size(); ++i){
            cloud->points[i].x  = in->get_vertices()[filteredPoints[i]].x();
            cloud->points[i].y  = in->get_vertices()[filteredPoints[i]].y();
            cloud->points[i].z  = in->get_vertices()[filteredPoints[i]].z();

        }

        std::vector<int> inliers;
        pcl::SampleConsensusModelCircle3D<pcl::PointXYZ>::Ptr model_circle(new pcl::SampleConsensusModelCircle3D<pcl::PointXYZ> (cloud));
        pcl::RandomSampleConsensus<pcl::PointXYZ> ransac(model_circle);
        ransac.setDistanceThreshold(thresholdRANSAC * in->get_diagonal());
        ransac.computeModel();
        ransac.getInliers(inliers);

        std::cout << "Num. inliers:" << inliers.size() << std::endl;

        //Obtain coefficients
        Eigen::VectorXf model_coefficients;
        ransac.getModelCoefficients(model_coefficients);

        center.push_back(model_coefficients[0]);
        center.push_back(model_coefficients[1]);
        center.push_back(model_coefficients[2]);

        directionVector.push_back(model_coefficients[4]);
        directionVector.push_back(model_coefficients[5]);
        directionVector.push_back(model_coefficients[6]);
}

Harris3D::Mesh* RotationalSymmetry::alignMeshes(Harris3D::Mesh* in, vector<float> center, vector<float> directionVector){
        float alpha = atof(prop->getProperty("alpha").c_str());
        int rotationFactor = atoi(prop->getProperty("rotation-factor").c_str());
        int iter_non_rigid = atoi(prop->getProperty("iter-non-rigid").c_str());
        float factor_spacing = atof(prop->getProperty("factor-spacing").c_str());
        float compatible_angle = atof(prop->getProperty("compatible-angle").c_str());
        int enable_icp = atoi(prop->getProperty("enable-icp").c_str());

        //Computar un kdtree para la malla de entrada. Se computa solo una vez y se reutiliza siempre
        SimpleMesh::Vertex* vertices = in->get_vertices();
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
        cloud->width = in->get_number_vertices();
        cloud->height = 1;
        cloud->points.resize(cloud->width * cloud->height);

        for(size_t i = 0; i < cloud->points.size(); i++){
            cloud->points[i].x = vertices[i].x();
            cloud->points[i].y = vertices[i].y();
            cloud->points[i].z = vertices[i].z();
        }

        pcl::KdTreeFLANN<pcl::PointXYZ> treeTarget;
        treeTarget.setInputCloud(cloud);


        Harris3D::Mesh * bigMesh = new Harris3D::Mesh();
        vector<SimpleMesh::Vertex> auxVertex;
        vector<SimpleMesh::Triangle> auxTriangles;

        int nLoop = 0;
        for (int s=360/rotationFactor; s<360;s+=360/rotationFactor){
            float angle = (float) s;

            //Find rotated mesh
            Harris3D::Mesh * auxMesh = rotateMesh(in, angle, center, directionVector);
            float spacing = computeSpacing(auxMesh);
            //cout << "Spacing:" <<  spacing << endl;

            if(enable_icp==1){
                //auxMesh contiene el resultado
                cout << "Empieza Sparse ICP" << endl;
                //Aplicamos SparseICP
                Eigen::Matrix3Xd source(3, in->get_number_vertices());
                Eigen::Matrix3Xd target(3, in->get_number_vertices());
                Eigen::Matrix3Xd normals(3, in->get_number_vertices());

                for(unsigned int j = 0; j < in->get_number_vertices(); j++){
                    source(0,j) = auxMesh->get_vertices()[j].x();
                    source(1,j) = auxMesh->get_vertices()[j].y();
                    source(2,j) = auxMesh->get_vertices()[j].z();

                    target(0,j) = in->get_vertices()[j].x();
                    target(1,j) = in->get_vertices()[j].y();
                    target(2,j) = in->get_vertices()[j].z();

                    normals(0,j) = in->get_vertices()[j].nx();
                    normals(1,j) = in->get_vertices()[j].ny();
                    normals(2,j) = in->get_vertices()[j].nz();
                }

                ICP::Parameters par;
                par.p = 0.5;
                par.f = ICP::PNORM;
                par.max_icp = 50;
                par.max_outer = 50;

                double a = 0;
                ICP::point_to_plane(source, target, normals, a, par);

                for(int j = 0; j < in->get_number_vertices(); j++){
                    auxMesh->get_vertices()[j].setX(source(0,j));
                    auxMesh->get_vertices()[j].setY(source(1,j));
                    auxMesh->get_vertices()[j].setZ(source(2,j));
                }
            }

            //stringstream ss;
            //ss << "output" << nLoop << ".off";
            //outputMesh(auxMesh, ss.str().c_str());

            //Calculate translations
            for (int w = 0; w < iter_non_rigid; w++){
                //Calculamos un kdtree para la malla a procesar
                vertices = auxMesh->get_vertices();
                pcl::PointCloud<pcl::PointXYZ>::Ptr cloudAux (new pcl::PointCloud<pcl::PointXYZ>);
                cloudAux->width = auxMesh->get_number_vertices();
                cloudAux->height = 1;
                cloudAux->points.resize(cloudAux->width * cloudAux->height);

                for(size_t i = 0; i < cloudAux->points.size(); i++){
                    cloudAux->points[i].x = vertices[i].x();
                    cloudAux->points[i].y = vertices[i].y();
                    cloudAux->points[i].z = vertices[i].z();
                }

                pcl::KdTreeFLANN<pcl::PointXYZ> treeSource;
                treeSource.setInputCloud(cloudAux);

                int n = 3 * auxMesh->get_number_vertices();
                Eigen::VectorXf x(n), b(n);
                Eigen::SparseMatrix<float> a(n,n);
                x.setZero();
                b.setZero();
                a.setZero();
                std::vector<Eigen::Triplet<float>> tripletList;
                //cout << "Llenando matriz esparsa" << endl;

                for (int j = 0; j<auxMesh->get_number_vertices(); j++){
                    SimpleMesh::Vertex v = auxMesh->get_vertices()[j];
                    int outIndex;

                    float aux = 0.0;
                    std::set<unsigned int>::iterator it;
                    for (it = v.get_vertices().begin(); it != v.get_vertices().end(); ++it){
                        float voronoiTriangle = computeVoronoiTriangle(auxMesh, j, *it);
                        float edge = computeDistance(&auxMesh->get_vertices()[j], &auxMesh->get_vertices()[*it]);

                        float coef = (2 * alpha * voronoiTriangle) / edge;
                        aux += coef;

                        tripletList.push_back(Eigen::Triplet<float>(3*j,3*(*it),-coef));
                        tripletList.push_back(Eigen::Triplet<float>(3*j+1,3*(*it)+1,-coef));
                        tripletList.push_back(Eigen::Triplet<float>(3*j+2,3*(*it)+2,-coef));
                    }

                    if(getClosestCompatiblePoint(auxMesh, j, in, treeTarget, spacing, outIndex, factor_spacing, compatible_angle)){
                        SimpleMesh::Vertex q = input->get_vertices()[outIndex];
                        float bCoef = 2*(1 - alpha)*computeOmega(treeSource, &q, &v, vertices, spacing)*computeVoronoiCell(auxMesh, j);
                        //if(v.red()>0.95)
                        //    bCoef /= 2;

                        aux += bCoef;

                        //B
                        b(3*j) = bCoef * (q.x() - v.x());
                        b(3*j+1) = bCoef * (q.y() - v.y());
                        b(3*j+2) = bCoef * (q.z() - v.z());
                    }

                    //A
                    tripletList.push_back(Eigen::Triplet<float>(3*j,3*j,aux));
                    tripletList.push_back(Eigen::Triplet<float>(3*j+1,3*j+1,aux));
                    tripletList.push_back(Eigen::Triplet<float>(3*j+2,3*j+2,aux));
                }

                //cout << "Optimizacion" << endl;
                //X
                Eigen::BiCGSTAB<Eigen::SparseMatrix<float> > solver;
                a.setFromTriplets(tripletList.begin(), tripletList.end());
                solver.compute(a);
                x = solver.solve(b);

                //Translate mesh
                for(int i = 0; i<auxMesh->get_number_vertices(); i++){
                    auxMesh->get_vertices()[i].setX(auxMesh->get_vertices()[i].x() + x[3*i]);
                    auxMesh->get_vertices()[i].setY(auxMesh->get_vertices()[i].y() + x[3*i+1]);
                    auxMesh->get_vertices()[i].setZ(auxMesh->get_vertices()[i].z() + x[3*i+2]);
                }
                auxMesh->compute_normals();
            }

            //Add vertices to newVertices
            int offset = nLoop * in->get_number_vertices();
            for (int i=0; i< in->get_number_vertices();i++){
                auxVertex.push_back(SimpleMesh::Vertex(auxMesh->get_vertices()[i].x(), auxMesh->get_vertices()[i].y(), auxMesh->get_vertices()[i].z()));
            }

            //Obtain faces
            int faceOffset = nLoop * in->get_number_triangles();
            for (int i=0; i< in->get_number_triangles();i++){
                SimpleMesh::Triangle auxFace;
                auxFace.add_vertex(auxMesh->get_triangles()[i].get_vertex_at(0)+offset);
                auxFace.add_vertex(auxMesh->get_triangles()[i].get_vertex_at(1)+offset);
                auxFace.add_vertex(auxMesh->get_triangles()[i].get_vertex_at(2)+offset);
                auxTriangles.push_back(auxFace);
            }
            //cout << "Loop: " << nLoop << endl;
            nLoop++;
            //break;
        }

        bigMesh->set_number_vertices(auxVertex.size());
        bigMesh->set_number_triangles(auxTriangles.size());

        for(int i = 0; i < auxVertex.size(); i++){
            bigMesh->add_vertex(i, auxVertex[i].x(), auxVertex[i].y(), auxVertex[i].z());
        }

        for(int i = 0; i < auxTriangles.size(); i++){
            bigMesh->add_triangle(i, auxTriangles[i].get_vertex_at(0), auxTriangles[i].get_vertex_at(1), auxTriangles[i].get_vertex_at(2));
        }

        return bigMesh;
}

Harris3D::Mesh* RotationalSymmetry::rotateMesh(Harris3D::Mesh* mesh, float angle, vector<float>& center, vector<float>& directionVector){
    //Create the mesh
    Harris3D::Mesh * auxMesh = new Harris3D::Mesh();
    auxMesh->set_number_vertices(mesh->get_number_vertices());
    auxMesh->set_number_triangles(mesh->get_number_triangles());

    //Obtain vertices
    for (int i=0; i< mesh->get_number_vertices();i++){
        float * rotatedVertex = Rotation::rotateVertex(mesh->get_vertices()[i].x()-center[0],mesh->get_vertices()[i].y()-center[1],mesh->get_vertices()[i].z()-center[2],angle,directionVector[0],directionVector[1],directionVector[2]);
        auxMesh->add_vertex(i, rotatedVertex[0]+center[0], rotatedVertex[1]+center[1], rotatedVertex[2]+center[2]);
        auxMesh->add_color(i, mesh->get_vertices()[i].red(), mesh->get_vertices()[i].green(), mesh->get_vertices()[i].blue());
    }

    //Obtain faces
    for (int i=0; i< mesh->get_number_triangles();i++){
         auxMesh->add_triangle(i,mesh->get_triangles()[i].get_vertex_at(0),mesh->get_triangles()[i].get_vertex_at(1),mesh->get_triangles()[i].get_vertex_at(2));
    }

    auxMesh->compute_normals();
    auxMesh->get_diagonal();
    return auxMesh;
}

//Calcular el average spacing
float RotationalSymmetry::computeSpacing(Harris3D::Mesh* mesh){
    int numEdges = 0;
    float sumEdges = 0.0;

    for(int i = 0; i < mesh->get_number_triangles(); i++){
        sumEdges += computeDistance(&mesh->get_vertices()[mesh->get_triangles()[i].get_vertex_at(0)], &mesh->get_vertices()[mesh->get_triangles()[i].get_vertex_at(1)]);
        sumEdges += computeDistance(&mesh->get_vertices()[mesh->get_triangles()[i].get_vertex_at(1)], &mesh->get_vertices()[mesh->get_triangles()[i].get_vertex_at(2)]);
        sumEdges += computeDistance(&mesh->get_vertices()[mesh->get_triangles()[i].get_vertex_at(0)], &mesh->get_vertices()[mesh->get_triangles()[i].get_vertex_at(2)]);
        numEdges += 3;
    }

    return sumEdges/numEdges;
}

//Calcula la distancia entre dos puntos
float RotationalSymmetry::computeDistance(SimpleMesh::Vertex *i, SimpleMesh::Vertex *j){

    float sum = 0.0;
    float d = i->x()-j->x();

    sum += d*d;
    d = i->y()-j->y();  sum += d*d;
    d = i->z()-j->z();  sum += d*d;


    return std::sqrt(sum);
}

//Calcula triangulo Voronoi de un edge
float RotationalSymmetry::computeVoronoiTriangle(Harris3D::Mesh* mesh, int i, int j){
    int neighbor1 = -1;
    int neighbor2 = -1;
    SimpleMesh::Vertex vertex1 = mesh->get_vertices()[i];
    SimpleMesh::Vertex vertex2 = mesh->get_vertices()[j];
    std::set<unsigned int>::iterator it1;
    std::set<unsigned int>::iterator it2;

    //find if they are neighbors
    for (it1 = vertex1.get_vertices().begin(); it1 != vertex1.get_vertices().end(); ++it1){
        for (it2 = vertex2.get_vertices().begin(); it2 != vertex2.get_vertices().end(); ++it2){
            if(*it1==*it2){
                if (neighbor1 == -1){
                    neighbor1 = *it1;
                }
                else{
                    neighbor2 = *it1;
                    break;
                }
            }
        }
    }

    if (neighbor1!=-1 && neighbor2!=-1){
        SimpleMesh::Vertex barycenter1 = computeBarycenter(mesh, vertex1.get_index(), vertex2.get_index(), neighbor1);
        SimpleMesh::Vertex barycenter2 = computeBarycenter(mesh, vertex1.get_index(), vertex2.get_index(), neighbor2);

        return computeAreaTriangle(&vertex1, &barycenter1, &barycenter2);
    }
    else{
        return 0.0;
    }
}

//Calcula el baricentro de un triángulo dado por 3 vértices
Harris3D::Vertex RotationalSymmetry::computeBarycenter(Harris3D::Mesh* mesh, int i, int j, int k){
    //find barycenter
    Harris3D::Vertex barycenter;
    barycenter.setX((mesh->get_vertices()[i].x() + mesh->get_vertices()[j].x() + mesh->get_vertices()[k].x())/3);
    barycenter.setY((mesh->get_vertices()[i].y() + mesh->get_vertices()[j].y() + mesh->get_vertices()[k].y())/3);
    barycenter.setZ((mesh->get_vertices()[i].z() + mesh->get_vertices()[j].z() + mesh->get_vertices()[k].z())/3);

    return barycenter;
}

//Calcula el área de un triángulo
float RotationalSymmetry::computeAreaTriangle(SimpleMesh::Vertex * i, SimpleMesh::Vertex * j, SimpleMesh::Vertex * k){
    //Heron formula
    float a = computeDistance(i,j);
    float b = computeDistance(i,k);
    float c = computeDistance(j,k);
    float s = (a+b+c)/2;

    return std::sqrt(s*(s-a)*(s-b)*(s-c));
}

//Buscar en el kdtree el punto mas cercano al vertice[index] en inputMesh
bool RotationalSymmetry::getClosestCompatiblePoint(Harris3D::Mesh *inputMesh, int index, Harris3D::Mesh* targetMesh, pcl::KdTreeFLANN<pcl::PointXYZ> kdtree, float spacing, int& outIndex, float factor_spacing, float compatible_angle){
    outIndex = -1;
    pcl::PointXYZ sp;
    SimpleMesh::Vertex vertex = inputMesh->get_vertices()[index];

    sp.x = vertex.x();
    sp.y = vertex.y();
    sp.z = vertex.z();

    float radius = factor_spacing * spacing;
    vector<int> idx;
    vector<float> dist;

    //Harris3D::Vertex* r = new Harris3D::Vertex;
    SimpleMesh::Vertex candidate;

    if(kdtree.radiusSearch(sp, radius, idx, dist)>0){
        for(size_t i = 0; i < idx.size(); i++){
            candidate = targetMesh->get_vertices()[idx[i]];
            if(computeNormalDotProduct(&vertex, &candidate)>compatible_angle){
                outIndex = idx[i];
                break;
            }
        }
    }

    if(outIndex!=-1)
        return true;

    return false;

}

//Calcula el producto interno entre dos normales
//Bug: Normales no están definidas en mallas creadas dinámicamente
float RotationalSymmetry::computeNormalDotProduct(SimpleMesh::Vertex * i, SimpleMesh::Vertex * j){

    return (i->nx()*j->nx()+i->ny()*j->ny()+i->nz()*j->nz());
}

float RotationalSymmetry::computeOmega(pcl::KdTreeFLANN<pcl::PointXYZ> kdtree, SimpleMesh::Vertex* q, SimpleMesh::Vertex* v, SimpleMesh::Vertex* verts, float h){
    pcl::PointXYZ sp;

    sp.x = q->x();
    sp.y = q->y();
    sp.z = q->z();

    int K = 1;
    vector<int> idx(K);
    vector<float> dist(K);

    SimpleMesh::Vertex* r = new SimpleMesh::Vertex;

    if(kdtree.nearestKSearch(sp, K, idx, dist)>0){
        r->setX(verts[idx[0]].x());
        r->setY(verts[idx[0]].y());
        r->setZ(verts[idx[0]].z());
    }

    float distF = computeDistance(v,r);
    return std::exp(-(distF*distF)/(h * h));
}

//Calcula celula Voronoi de un vértice
float RotationalSymmetry::computeVoronoiCell(Harris3D::Mesh* mesh, int i){
    float area = 0.0;
    std::set<unsigned int>::iterator it;
    SimpleMesh::Vertex vertex = mesh->get_vertices()[i];

    for (it = vertex.get_vertices().begin(); it != vertex.get_vertices().end(); ++it){
        float aux = computeVoronoiTriangle(mesh, vertex.get_index(), *it);
        if (aux > 0.0){
            area += aux;
        }
    }

    return area;
}

SimpleMesh::Mesh* RotationalSymmetry::computeCompletion(){

    Harris3D::Mesh* resultOpenVDB = nullptr;
    //Preservar la malla original - Realizar copia profunda
    int performOpenVDB = atoi(prop->getProperty("performOpenVDB").c_str());
    if(performOpenVDB){
        std::cout << "OpenVDB" << std::endl;
        resultOpenVDB = applyOpenVDB();
    }

    cout << resultOpenVDB->get_number_vertices() << endl;
    Harris3D::Mesh* resultDecimation = nullptr;

    int performSimplification = atoi(prop->getProperty("performSimplification").c_str());
    if(performSimplification){
        std::cout << "Decimate" << std::endl;
        if(performOpenVDB)
            resultDecimation = decimate(resultOpenVDB);
        else
            resultDecimation = decimate(input);
    }

    cout << resultDecimation->get_number_vertices() << endl;
    SimpleMesh::IO::OFFWriter w("output.off");
    w.write_mesh(*resultDecimation);

    Harris3D::Mesh* resultPrepro = input;
    if(performOpenVDB)
        resultPrepro = resultOpenVDB;
    if(performSimplification)
        resultPrepro = resultDecimation;



    Harris3D::Mesh* resultClass = resultPrepro;
    if(performOpenVDB || performSimplification){
        std::cout << "Transfer Classification" << std::endl;
        if(input->has_color())
            resultClass = transferClassification(resultPrepro, input);
    }

        std::cout << "RANSAC" << std::endl;
        vector<float> center;
        vector<float> directionVector;

        applyRANSAC(resultClass, center, directionVector);

        std::cout << "Alignment" << std::endl;
        SimpleMesh::Mesh* output = alignMeshes(resultPrepro, center, directionVector);

        if(performOpenVDB)
            delete resultOpenVDB;
        if(performSimplification)
            delete resultDecimation;

        if(input->has_color())
            delete resultClass;
        return output;
}

}
