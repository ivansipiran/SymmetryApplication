#include <cassert>
#include <cstdlib>
#include <cfloat>
#include <algorithm>
#include <vector>

#include <Eigen/Dense>
#include "CGAL/Simple_cartesian.h"
#include <CGAL/point_generators_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_3.h>

#include "Harris3D/Mesh.h"
#include "Harris3D/HarrisDetector.h"
//#include "Util/Clock.h"

typedef CGAL::Simple_cartesian<float> K;
typedef K::Point_3 Point_3;
typedef CGAL::Search_traits_3<K> Traits;
typedef CGAL::Random_points_in_cube_3<Point_3>       Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Kd_tree<Traits> Tree;
typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;

using namespace Eigen;

namespace Harris3D{

struct HarrisResult{
    SimpleMesh::Vertex* vertex;
    float response;
};

HarrisDetector :: HarrisDetector(){
	typeNeighborhood = ADAPTIVE;
	fractionDiagonal = 0.01;
	numberRingNeighbor = 0;

	k = 0.04;

	numberRingsDetection = 1;

	typeSelection = FRACTION;
	paramSelection = 0.01;

	filteringSteps = 0;
	object = NULL;
	prop = NULL;
}

HarrisDetector :: HarrisDetector(Mesh* obj, Util::PropertySet* pr){
	object = obj;
	prop = pr;
    object->initializeVertexInformation();
	processOptions();
}

void HarrisDetector :: processOptions(){
	/*We evaluate the properties*/
	//typeNeighborhood
	std::string val = prop->getProperty("type-neighborhood");
	if(val.empty()){
		std::cout << "Option type-neighborhood was not set. ADAPTIVE is assumed." << std::endl;
		typeNeighborhood = ADAPTIVE;
	}
	else{
		if(!val.compare("spatial"))
			typeNeighborhood = SPATIAL;
		else if(!val.compare("rings"))
			typeNeighborhood = RINGS;
		else if(!val.compare("adaptive"))
			typeNeighborhood = ADAPTIVE;
		else{
			std::cout << "Option type-neighborhood is not recognized. ADAPTIVE is assumed." << std::endl;
			typeNeighborhood = ADAPTIVE;
		}
	}

	//typeSelection
	val = prop->getProperty("interest-points-selection");
	if(val.empty()){
		std::cout << "Option interest-points-detection was not set. FRACTION is assumed." << std::endl;
		typeSelection = FRACTION;
	}else{
		if(!val.compare("clustering"))
			typeSelection = CLUSTERING;
		else if(!val.compare("fraction"))
			typeSelection = FRACTION;
		else if(!val.compare("number"))
			typeSelection = NUMBER;
		else{
			std::cout << "Option interest-points-detection is not recognized. FRACTION is assumed." << std::endl;
		}
	}

	//parameter-neighborhood
	val = prop->getProperty("parameter-neighborhood");
	if(typeNeighborhood == SPATIAL || typeNeighborhood == ADAPTIVE){
		if(val.empty()){
			std::cout << "Using SPATIAL or ADAPTIVE: Option parameter-neighborhood was not set. Using 0.01 by default" << std::endl;
			fractionDiagonal = 0.01;
		}
		else
			fractionDiagonal = atof(val.c_str());
	}else if(typeNeighborhood == RINGS){
		if(val.empty()){
			std::cout << "Using RINGS: Option parameter-neighborhood was not set. Using 2 by default" << std::endl;
			numberRingNeighbor = 2;
		}
		else
			numberRingNeighbor = atoi(val.c_str());
	}

	//k
	val = prop->getProperty("K");
	if(val.empty()){
		std::cout << "Parameter K was not set. Using 0.04 by default" << std::endl;
		k = 0.04;
	}
	else{
		k = atof(val.c_str());
	}

	//numberRingsDetection
	val = prop->getProperty("ring-maxima-detection");
	if(val.empty()){
		std::cout << "Parameter ring-maxima-detection was not set. Using 1 by default" << std::endl;
		numberRingsDetection = 1;
	}else{
		numberRingsDetection = atoi(val.c_str());
	}

	//parameter-selection
	val = prop->getProperty("parameter-selection");
	if(typeSelection == FRACTION || typeSelection == CLUSTERING){
		if(val.empty()){
			std::cout << "Parameter parameter-selection was not set. Using 0.01 by default." << std::endl;
			paramSelection = 0.01;
		}else{
			paramSelection = atof(val.c_str());
		}
	}else if(typeSelection == NUMBER){
		if(val.empty()){
			std::cout << "Parameter parameter-selection was not set. Using 100 by default." << std::endl;
			numberKeypoints = 100;
		}else{
			numberKeypoints = atoi(val.c_str());
		}
	}

	val = prop->getProperty("filtering-steps");
	if(val.empty()){
		std::cout << "Parameter filtering-steps was not set. Using 0 by default" << std::endl;
		filteringSteps = 0;
	}else{
		filteringSteps = atoi(val.c_str());
	}
}

void HarrisDetector :: showOptions(){
	std::cout << "type-neighborhood = ";
	switch(typeNeighborhood){
		case SPATIAL: std::cout << "SPATIAL" << std::endl;
					  break;
		case ADAPTIVE: std::cout << "ADAPTIVE" << std::endl;
						break;
		case RINGS:	std::cout << "RINGS" << std::endl;
					break;
	}

	std::cout << "parameter-neighborhood = ";
	if(typeNeighborhood == SPATIAL || typeNeighborhood == ADAPTIVE)
		std::cout << fractionDiagonal << std::endl;
	else
		std::cout << numberRingNeighbor << std::endl;

	std::cout <<"K = "<< k << std::endl;
	std::cout <<"ring-maxima-detection = " << numberRingsDetection << std::endl;

	std::cout << "interest-points-selection = ";
	switch(typeSelection){
		case FRACTION: std::cout << "FRACTION" << std::endl;
						break;
		case CLUSTERING: std::cout << "CLUSTERING" << std::endl;
						break;
	}

	std::cout << "parameter-selection = " << paramSelection << std::endl;
}

bool comp(HarrisResult i, HarrisResult j){
    return i.response > j.response;
}

void HarrisDetector :: detectInterestPoints(std::vector<unsigned int>& interestPoints){
	assert(object);

	float max = 0.0;
	int rad;
	int cont = 0;
    float diag = object->get_diagonal();

    SimpleMesh::Vertex* vertices = object->get_vertices();
    SimpleMesh::Triangle* faces = object->get_triangles();

	std::cout << "Starting Harris" << std::endl;
	Tree tree;

	if(typeNeighborhood == SPATIAL){
		std::vector<Point_3> points;
        for(int i = 0; i < object->get_number_vertices(); i++){
			points.push_back(Point_3(vertices[i].x(), vertices[i].y(), vertices[i].z()));
		}
		tree.insert(points.begin(), points.end());
	}

    float* responses = new float[object->get_number_vertices()];
    memset(responses, 0, sizeof(float)*object->get_number_vertices());
	//Process each vertex

	/*This loop computes the Harris response for each vertex*/
	//#pragma omp parallel for
    for(int i = 0; i< object->get_number_vertices(); i++){
		//cout << "Vertex " << i << "->";
        std::vector<SimpleMesh::Vertex*> set;

		//Calcular radio adaptativo
		if(typeNeighborhood == ADAPTIVE){
            //rad = vertices[i].getRadius(vertices, diag * fractionDiagonal, set);
            rad = object->getRadiusAtVertex(i, diag * fractionDiagonal, set);
		}
		else if(typeNeighborhood == SPATIAL){
			Point_3 query(vertices[i].x(), vertices[i].y(), vertices[i].z());
			Fuzzy_sphere fs(query, diag * fractionDiagonal, 0.0);

			std::vector<Point_3> puntos;
			tree.search(back_inserter(puntos), fs);

			for(int i = 0; i < puntos.size(); i++){
                SimpleMesh::Vertex* p = new SimpleMesh::Vertex(puntos[i].x(), puntos[i].y(), puntos[i].z());
				set.push_back(p);
			}
		}else if(typeNeighborhood == RINGS){
            //vertices[i].getNeighborhood(numberRingNeighbor, set, vertices);
            object->getNeighborhoodAtVertex(i, numberRingNeighbor, set);
		}

		if(set.size() < 6){
			cont++;
			//vertices[i].setResponse(DBL_MIN);
			responses[i] = FLT_MIN;
			continue;
		}
		//Process "set", first, calculate the centroid
		float xc = 0, yc = 0, zc = 0;
		for(int j = 0; j< set.size();j++){
			xc += set[j]->x();
			yc += set[j]->y();
			zc += set[j]->z();
		}

		xc /= set.size();
		yc /= set.size();
		zc /= set.size();

		//Translate the vertex, in order the centroid is in [0 0 0]
		for(int j = 0; j< set.size();j++){
			set[j]->setX(set[j]->x() - xc);
			set[j]->setY(set[j]->y() - yc);
			set[j]->setZ(set[j]->z() - zc);
		}

		//Aplicar PCA para encontrar una pose de la nube de puntos de manera que la mayor dispersión de los puntos esté en el plano XY
		//La media de las 3 coordenadas ya es (0, 0, 0), así que en realidad no es necesario calcularla, directamente calculamos la
		//matriz de covarianza
		float A[9];
		memset(A, 0, sizeof(float)*9);

		for(int j = 0; j < set.size(); j++){
			float x = set[j]->x();
			float y = set[j]->y();
			float z = set[j]->z();

			A[0] += x*x;	A[1] += x*y;	A[2] += x*z;
							A[4] += y*y;	A[5] += y*z;
											A[8] += z*z;
		}
		A[3] = A[1];	A[6] = A[2];	A[7] = A[5];
		for(int j = 0; j < 9; j++)
			A[j] /= (set.size()-1);

		Matrix3d m;
		m(0,0) = A[0];	m(0,1) = A[1];	m(0,2) = A[2];
		m(1,0) = A[3];	m(1,1) = A[4];	m(1,2) = A[5];
		m(2,0) = A[6];	m(2,1) = A[7];	m(2,2) = A[8];

		SelfAdjointEigenSolver<MatrixXd>	eigensolver(m);
		Matrix3d evec = eigensolver.eigenvectors();

		//cout << "The eigenvalues are " << eigensolver.eigenvalues() << endl;
		float x_1 = evec(0,2);	float x_2 = evec(1,2);	float x_3 = evec(2,2);
		float y_1 = evec(0,1);	float y_2 = evec(1,1);	float y_3 = evec(2,1);
		float z_1 = evec(0,0);	float z_2 = evec(1,0);	float z_3 = evec(2,0);

		float x2 = set[0]->x() - xc;
		float y2 = set[0]->y() - yc;
		float z2 = set[0]->z() - zc;

		if((z_1*x2 + z_2*y2 + z_3*z2) < 0){
			z_1 = -z_1;
			z_2 = -z_2;
			z_3 = -z_3;

			float aux_x1 = x_1;
			float aux_x2 = x_2;
			float aux_x3 = x_3;

			x_1 = y_1;
			x_2 = y_2;
			x_3 = y_3;
			y_1 = aux_x1;
			y_2 = aux_x2;
			y_3 = aux_x3;

		}
		//Realizamos la rotacion, con el nuevo sistema de coordenadas
		for(int j = 0; j < set.size(); j++){
			float x = set[j]->x();
			float y = set[j]->y();
			float z = set[j]->z();

			set[j]->setX(x_1*x + x_2*y + x_3*z);
			set[j]->setY(y_1*x + y_2*y + y_3*z);
			set[j]->setZ(z_1*x + z_2*y + z_3*z);
		}

		//Movemos todos los puntos para que el punto de analisis se encuentre en el origen del plano XY. Solo movemos en las coordenadas X e Y
		float x = set[0]->x();
		float y = set[0]->y();

		for(int j = 0; j < set.size(); j++){
			set[j]->setX(set[j]->x() - x);
			set[j]->setY(set[j]->y() - y);
		}

		//Fit a quadratic surface
		float C[36];
		float D[6];

		memset(C, 0, sizeof(float)*36);
		memset(D, 0, sizeof(float)*6);

		for(int j = 0; j < set.size(); j++){
			float x = set[j]->x();
			float y = set[j]->y();
			float z = set[j]->z();

			float x2 = x*x;
			float y2 = y*y;
			float x3 = x2*x;
			float y3 = y2*y;

			C[0] += x*x3;	C[1] += x3*y;	C[2] += x2*y2;	C[3] += x3;		C[4] += x2*y;	C[5] += x2;
							C[7] += x2*y2;	C[8] += x*y3;	C[9] += x2*y;	C[10] += x*y2;	C[11] += x*y;
											C[14] += y*y3;	C[15] += x*y2;	C[16] += y3;	C[17] += y2;
															C[21] += x2;	C[22] += x*y;	C[23] += x;
																			C[28] += y2;	C[29] += y;
			D[0] += z*x2;	D[1] += z*x*y;	D[2] += z*y2;	D[3] += z*x;	D[4] += z*y;	D[5] += z;
		}

		C[6] = C[1];
		C[12] = C[2];	C[13] = C[8];
		C[18] = C[3];	C[19] = C[9];	C[20] = C[15];
		C[24] = C[4];	C[25] = C[10];	C[26] = C[16];	C[27] = C[22];
		C[30] = C[5];	C[31] = C[11];	C[32] = C[17];	C[33] = C[23];	C[34] = C[29];

		C[35] = (float)set.size();

		MatrixXd m1(6,6);
		VectorXd  b1(6);

		for(int j = 0; j < 6; j++){
			for(int k = 0; k < 6; k++)
				m1(j,k) = C[j*6 + k];
			b1(j) = D[j];
		}

		ColPivHouseholderQR<MatrixXd> dec(m1);
		VectorXd x1 = dec.solve(b1);

		float c20_2 = x1(0);	//p1
		float c11 = x1(1);		//p2
		float c02_2 = x1(2);	//p3
		float c10 = x1(3);		//p4
		float c01 = x1(4);		//p5
		float c0 = x1(5);		//p6

		//double c20 = c20_2*2;
		//double c02 = c02_2*2;
		float c20 = c20_2;
		float c02 = c02_2;

		float fx2 = c10*c10 + 2*c20*c20 + 2*c11*c11; //A
		float fy2 = c01*c01 + 2*c11*c11 + 2*c02*c02; //B
		float fxfy = c10*c01 + 2*c20*c11 + 2*c11*c02; //C

		//double k = 0.04;
		float resp = fx2*fy2 - fxfy*fxfy - k*(fx2 + fy2)*(fx2 + fy2);

		//vertices[i].setResponse(resp);
        //cout << resp << endl;
		responses[i] = resp;

		if(resp > max)
			max = resp;

		for(int j = 0; j < set.size(); j++)
			delete set[j];
		set.clear();
	}

    float* auxResp = new float[object->get_number_vertices()];

	std::cout << "Executing filtering in " << filteringSteps << " steps" << std::endl;
	//Perform the filtering before the use of the responses
	for(int i = 0; i < filteringSteps; i++){
        memcpy(auxResp, responses, sizeof(float)*object->get_number_vertices());

        for(int j = 0; j < object->get_number_vertices(); j++){
			std::set<unsigned int> neighbors = vertices[j].get_vertices();
			float sum=0.0;
			for(std::set<unsigned int>::iterator it = neighbors.begin(); it!=neighbors.end(); it++){
				sum += auxResp[*it];
			}
			responses[j] = sum/neighbors.size();
		}
	}

	delete[] auxResp;
    //for(int i = 0; i < object->get_number_vertices(); i++)
    //	vertices[i].setResponse(responses[i]);
    object->setResponses(responses);

	delete[] responses;
    std::vector<HarrisResult> candidatePoints;

	//Search for local maximum
    for(unsigned int i = 0; i< object->get_number_vertices();i++){
        //vertices[i].processMaximum(vertices, numberRingsDetection);
        object->processMaximumAtVertex(i, numberRingsDetection);
        if(object->getInterestAtVertex(i)){
            HarrisResult aux;
            aux.vertex = &vertices[i];
            aux.response = object->getResponseAt(i);
            candidatePoints.push_back(aux);
		}
	}

	std::cout << "Candidates:" << candidatePoints.size() << std::endl;
	std::sort(candidatePoints.begin(), candidatePoints.end(), comp);

	if(typeSelection == FRACTION){
		//Seleccionar los puntos de mayor respuesta
        int numPoints = paramSelection * object->get_number_vertices() > candidatePoints.size()? candidatePoints.size():paramSelection*object->get_number_vertices();
        for(int i = 0; i < numPoints; i++){
            SimpleMesh::Vertex* aux = candidatePoints[i].vertex;
            interestPoints.push_back(aux->get_index());
        }
	}else if(typeSelection == CLUSTERING){
		//Aplicar proceso de Clustering
		for(int i = 0; i < candidatePoints.size(); i++){
			bool found = false;
			int j = 0;
			while(j < interestPoints.size() && !found){
                SimpleMesh::Vertex* aux = candidatePoints[i].vertex;
                float distX = vertices[interestPoints[j]].x() - aux->x();
                float distY = vertices[interestPoints[j]].y() - aux->y();
                float distZ = vertices[interestPoints[j]].z() - aux->z();

				if(sqrt(distX*distX + distY*distY + distZ*distZ) < (paramSelection * diag))
					found = true;
				j++;
			}
            if(!found){
                SimpleMesh::Vertex* aux = candidatePoints[i].vertex;
                interestPoints.push_back(aux->get_index());
            }
		}
	}else if(typeSelection == NUMBER){
		if(candidatePoints.size() < numberKeypoints)
			numberKeypoints = candidatePoints.size();
        for(int i = 0; i < numberKeypoints; i++){
            SimpleMesh::Vertex* aux = candidatePoints[i].vertex;
            interestPoints.push_back(aux->get_index());
        }
	}

	//for(int i = 0; i < interestPoints.size(); i++)
	//	cout << interestPoints[i].getIndex() << endl;
	//cout << "Interest points:" <<interestPoints.size() << endl;
}

}
