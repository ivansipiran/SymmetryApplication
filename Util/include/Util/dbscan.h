#include <vector>
#include <list>
#include <cmath>
#include <cstring>

using namespace std;

//Implemented with column-based indexing, so if you need to use it within C/C++, the 
//set of points should have the x-coordinates in the first row and the y-coordinates in
//the second row.

double euclideanDistance(double* array, int index1, int index2, int rows, int cols){
	double sum = 0.0;
	for(int i = 0; i < cols; i++){
		double diff = array[index1*cols + i] - array[index2*cols + i];
		sum += diff*diff;
	}
	return sqrt(sum);
}

void getNeighbors(double* array, int rows, int cols, int index, double eps, list<int>& indices){
	
	for(int i = 0; i < rows; i++){
		double dist = euclideanDistance(array, index, i, rows, cols);
		if(dist <= eps)
			indices.push_back(i);
	}
	
	
}

double getEps(double* array, int rows, int cols, int k){
	
	double prod = 1.0;
	for(int i = 0; i < cols; i++){
		double max = array[i];
		double min = array[i];
		for(int j = 1; j < rows; j++){
			if(array[j*cols + i] > max)
				max = array[j*cols + i];
			if(array[j*cols + i] < min)
				min = array[j*cols + i];
		}
		prod *= (max - min);
	}
	
	return pow(((prod*k)/(rows*pow(M_PI, ((double)cols)/2))), 1.0/cols);
}

void dbscan(double* array, int rows, int cols, int k, double eps, vector<int>& cla, vector<int>& type){
	int no = 1;
	int* visited = new int[rows];
	
	memset(visited, 0, sizeof(int)*rows);
	
	for(int i = 0; i < rows; i++){
		if(!visited[i]){
			list<int> indices;
			getNeighbors(array, rows, cols, i, eps, indices);
			
			int numNeighbors = indices.size();
			if(numNeighbors > 1 && numNeighbors < (k+1)){
				type[i] = 0;
				cla[i] = 0;
			}
			if(numNeighbors == 1){
				type[i] = -1;
				cla[i] = -1;
				visited[i] = 1;
			}
			
			if(numNeighbors > k){
				type[i] = 1;
				list<int>::iterator it = indices.begin();
				while(it!=indices.end()){
					cla[*it] = no;
					it++;
				}
				
				while(!indices.empty()){
					list<int> indicesAux;
					
					int indAux = indices.front();
					getNeighbors(array, rows, cols, indAux, eps, indicesAux);
					indices.pop_front();
					
					int numN = indicesAux.size();
					if(numN > 1){
						list<int>::iterator it1 = indicesAux.begin();
						while(it1!=indicesAux.end()){
							if(!visited[*it1]){
								visited[*it1] = 1;
								indices.push_back(*it1);
								cla[*it1] = no;
							}
							it1++;
						}
						
						if(numN > k)
							type[indAux] = 1;
						else
							type[indAux] = 0;	
						
					}
				}
				no++;
			}
		}
	}
	
	for(int i = 0; i < rows; i++){
		if(cla[i]==0){
			cla[i] = -1;
			type[i] = -1;
		}
	}
	
	delete[] visited;
}
