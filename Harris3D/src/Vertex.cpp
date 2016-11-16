#include "Harris3D/Vertex.h"
#include "Util/util.h"
#include <iostream>
#include <cmath>
#include <map>

namespace Harris3D{

void Vertex::getNeighborhood(int rad, std::vector<Vertex*>& V, Vertex* vertices){
	std::queue<Vertex*> Q;
	std::vector<Vertex*> marked;

	Q.push(this);
	this->setMark(true);
	this->setDepth(0);
	marked.push_back(this);

	while(!Q.empty()){
		Vertex* v0 = Q.front();
		Q.pop();
		V.push_back(new Vertex(v0->x(), v0->y(), v0->z())); //Indeed, copy vertex information rather than return the same vertex

		int dep = v0->getDepth();
		if(dep <= rad){
			std::set<unsigned int> listVertices = v0->get_vertices();
			std::set<unsigned int> :: iterator it;
			for(it = listVertices.begin(); it!=listVertices.end(); it++){
				Vertex* v1 = &vertices[*it];
					if(!v1->isMarked()){
						Q.push(v1);
						v1->setMark(true);
						v1->setDepth(dep + 1);
						marked.push_back(v1);
					}
				}
			}
	}

	std::vector<Vertex*>::iterator ini = marked.begin();

	while(ini<marked.end()){
		(*ini)->setMark(false);
		(*ini)->setDepth(0);
		ini++;
	}
}

int Vertex :: getRadius(Vertex*  vertices, float radius, std::vector<Vertex*>& V){
	std::vector<Vertex*> marked; //Store the marked vertices
	std::map<unsigned int, float> distances; //Store the distances relatives to the current vertex
	std::map<unsigned int, Vertex*> markedRing; //Elements in a new ring
	std::queue<Vertex*> Q;
	float maxDistance = 0.0;
	int rad = -1;

	Q.push(this);
	this->setMark(true);
	markedRing.insert(std::pair<unsigned int, Vertex*>(this->index, this));

	distances[this->index] = 0.0;

	while(!Q.empty()){
		Vertex* v0 = Q.front();
		Q.pop();

		int dep = v0->getDepth();
		if(dep != rad){ //First vertex in the new ring
			std::map<unsigned int, Vertex*>::iterator it;
			float max = 0.0;

			//Mark the previous ring
			for(it = markedRing.begin(); it!=markedRing.end(); it++){
				Vertex* mar = (*it).second;
				mar->setMark(true);
				marked.push_back(mar);
				V.push_back(new Vertex(mar->x(), mar->y(), mar->z()));
				if(distances[(*it).first] > max)
					max = distances[(*it).first];
			}

			rad++;
			markedRing.clear();
			maxDistance = max;
			if(maxDistance > radius)
				break;
		}

		std::set<unsigned int> listVertices = v0->get_vertices();
		std::set<unsigned int> :: iterator it;

		for(it = listVertices.begin(); it!=listVertices.end(); it++){
				Vertex* v1 = &vertices[*it];
				if(!v1->isMarked()){
					if(distances[v1->get_index()] == 0.0){ //Distance is not set
						Q.push(v1);
						v1->setDepth(dep + 1);
					}
					markedRing.insert(std::pair<unsigned int, Vertex*>(v1->get_index(), v1));
					float dist = v0->distanceL2(v1);
					float newDistance = distances[v0->get_index()] + dist;
					if(distances[v1->get_index()] == 0.0){ //First time on this vertex
						distances[v1->get_index()] = newDistance;
					}else if(newDistance  < distances[v1->get_index()]){
						distances[v1->get_index()] = newDistance;
					}
				}
			}
		//}
	}

	if(!markedRing.empty()){
			std::map<unsigned int, Vertex*>::iterator it;
			float max = 0.0;


			for(it = markedRing.begin(); it!=markedRing.end(); it++){
				Vertex* mar = (*it).second;
				mar->setMark(true);
				marked.push_back(mar);
				V.push_back(new Vertex(mar->x(), mar->y(), mar->z()));
				if(distances[(*it).first] > max)
					max = distances[(*it).first];
			}

			rad++;
			markedRing.clear();
			maxDistance = max;

	}

	//Unmark all vertices
	std::vector<Vertex*>::iterator ini = marked.begin();
	while(ini  < marked.end()){
		(*ini)->setMark(false);
		(*ini)->setDepth(0);
		ini++;
	}
	return rad;
}

void Vertex::processMaximum(Vertex* vertices, int numRings){
		std::set<unsigned int> :: iterator it;
		for(it = adjacentVertices.begin(); it!=adjacentVertices.end(); it++){
			Vertex* v1 = &vertices[*it];
			if(v1!=this){
				if(response < v1->getResponse())
					return;
			}
		}
	isInterest = true;
}

void Vertex :: getPatch(Vertex* vertices, std::vector<unsigned int> indices, std::set<unsigned int>& returned, std::set<unsigned int>& faceR, float radius, Vertex center){
	std::set<unsigned int> waiting;
	std::queue<unsigned int> visited;
	visited.push(this->index); //Este vertice va a la cola

	waiting.insert(indices.begin(), indices.end());
	waiting.erase(this->index); // Eliminamos este vertice del conjunto de faltantes

	while(!waiting.empty() && !visited.empty()){
		int ind = visited.front();
		visited.pop();

		if(!vertices[ind].isMarked()){
			returned.insert(ind);
			vertices[ind].setMark(true);
			waiting.erase(ind);

			std::set<unsigned int> listVertices = vertices[ind].get_vertices();
			std::set<unsigned int> :: iterator it;

			for(it = listVertices.begin(); it!=listVertices.end(); it++){
				unsigned int ind = *it;
				if(!vertices[ind].isMarked()){
					float distX = vertices[ind].x() - center.x();
					float distY = vertices[ind].y() - center.y();
					float distZ = vertices[ind].z() - center.z();
					float dist = sqrt(distX*distX + distY*distY + distZ*distZ);
					if(dist < radius)
						visited.push(*it);
				}
			}

			std::vector<unsigned int> fac = vertices[ind].get_triangles();
			std::vector<unsigned int>::iterator it1;
			for(it1 = fac.begin(); it1!=fac.end(); it1++)
				faceR.insert(*it1);
		}

	}

	std::set<unsigned int>::iterator it2;
	for(it2 = returned.begin(); it2!=returned.end(); it2++)
		vertices[*it2].setMark(false);

}

std::ostream& operator<<(std::ostream& out, Vertex& point){
	out << point.x() <<" "<<point.y()<<" "<<point.z()<<std::endl;

	return out;
}

float Vertex::distanceL2(Vertex* v1){
    float val = (x()-v1->x())*(x()-v1->x()) + (y()-v1->y())*(y()-v1->y()) + (z()-v1->z())*(z()-v1->z());
    return sqrt(val);
}
}
