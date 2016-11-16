#ifndef __HARRISDETECTOR_H_
#define __HARRISDETECTOR_H_

#include "Harris3D/harris_global.h"
#include "Harris3D/Mesh.h"
//#include "Vertex.h"
#include "Util/PropertySet.h"

namespace Harris3D{

class HARRIS_API HarrisDetector{
	private:
		Mesh* object;
        Util::PropertySet* prop;

		int typeNeighborhood;
		float fractionDiagonal;
		int numberRingNeighbor;

		float k;

		int numberRingsDetection;

		int typeSelection;
		float paramSelection;
		int numberKeypoints;

		int filteringSteps;

		void processOptions();

	public:
		enum _type_neighborhood{SPATIAL, ADAPTIVE, RINGS};
		enum _type_selection{FRACTION, CLUSTERING, NUMBER};

		HarrisDetector();
        HarrisDetector(Mesh* obj, Util::PropertySet* pr);
		virtual ~HarrisDetector() {}

		inline void setMesh(Mesh* obj) { object = obj;}
        inline void setProperties(Util::PropertySet* pr) { prop = pr; processOptions();}

		void showOptions();

		void detectInterestPoints(std::vector<unsigned int>& interestPoints);

};
}
#endif
