#ifndef __PROPERTYSET_H_
#define __PROPERTYSET_H_

#include <map>
#include <string>
#include "util_global.h"

namespace Util{

class UTIL_API PropertySet{
	private:
		std::map<std::string, std::string> propertyList;
		std::string trim(std::string& input);

	public:
        PropertySet() {}
        virtual ~PropertySet() {}

		void load(std::string filename);
		std::string getProperty(std::string property);
		void show();
        void addProperty(std::string name, std::string value);
};
}
#endif

