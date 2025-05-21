#include <regex>
#include <string>
#include "FSPath.h"


bool FSPath::isAbsolute(const char* path)
{
#ifdef WIN32
	return std::regex_search(path, std::regex("^([A-Z]|[a-z])\\:[\\\\\\/]"));
#else
	return std::regex_search(path, std::regex("^\\/"));
#endif
}

bool FSPath::isPath(const char* path)
{
	return std::regex_search(path, std::regex("\\/|\\\\"));
}

void FSPath::filePath(char* filename, char* path)
{
	std::string name(filename);

	int index = name.find_last_of("/\\");

    if(index == std::string::npos)
    {
        strcpy(path,"");
        return;
    }

#ifdef WIN32
    char sep = '\\';
#else
    char sep = '/';
#endif

    sprintf(path,"%s%c", name.substr(0,index).c_str(), sep);
}


















