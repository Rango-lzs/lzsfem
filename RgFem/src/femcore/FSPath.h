
#pragma once

class FSPath
{
public:
	static bool isAbsolute(const char* path);
	static bool isPath(const char* path);

	static void filePath(char* filename, char* path);
};

