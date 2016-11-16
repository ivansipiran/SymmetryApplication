#include <iostream>
#include <string>
#include <fstream>
#include <locale>
#include <cstdlib>
#include <cassert>
#include "Util/PropertySet.h"

using namespace std;

namespace Util{

void PropertySet :: load (string filename) {
	ifstream in(filename.c_str());
	assert(in);
	int numline = 0;

	while(in){
		string line;
		getline(in, line);
		numline++;
		line = trim(line);
		if(line.size() > 0){
			if(line[0] != '#'){
				size_t pos = line.find('=');
				if(pos != string::npos){
					string a = line.substr(0,pos);
					string b = line.substr(pos + 1);
					propertyList[trim(a)] = trim(b);
				}else{
					cout << "Bad formatting of property in line " << numline << endl;
					exit(EXIT_FAILURE);
				}
			}
		}

	}
	in.close();
}

void PropertySet :: addProperty (std::string name, std::string value){
    propertyList[trim(name)] = trim(value);
}

string PropertySet :: getProperty (string property) {
	string value = propertyList[property];
	if(value.empty())
		propertyList.erase(property);
	return value;
}

void PropertySet :: show () {
	map<string, string>::iterator it;

	for(it = propertyList.begin(); it != propertyList.end(); it++){
		cout << (*it).first << " => " << (*it).second << endl;
	}
}

string PropertySet :: trim (string& input){
	string ret = input;

	string::iterator p = ret.begin();
	while(p != ret.end() && *p == ' '){
		p = ret.erase(p);
	}

	p = ret.end();
	while(p != ret.begin() && !isalnum(*(p-1))){
		p = ret.erase(p - 1);
	}

	return ret;
}

}
