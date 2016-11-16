#include <iostream>
#include <cmath>
#include "Util/util.h"

using namespace std;

namespace Util{

void myPause(){
	cout<<"Press any key to continue...";
	//getchar();
}

void myPause(const char* message){
	cout<<message<<endl;
	myPause();
}

int sgn(double num){
	if(num >= 0.0)
		return 1;
	else
		return -1;
}

void skipline(istream &in){
	char c;
	while(in>>noskipws>>c && c!='\n');
	in>>skipws;
}

bool isNumber(string str){
	int length = str.size();
	bool isChar = true;

	int i = 0;

	while(isChar && i < length){
		if(!isdigit(str[i]))
			isChar = false;
		i++;
	}

	return isChar;
}
}
