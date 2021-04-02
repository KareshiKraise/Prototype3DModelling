#ifndef _LOADER_
#define _LOADER_
#include <vector>
#include <string>


	class IOManager
	{
	public:
		static bool readFileToBuffer(std::string filePath, std::vector<unsigned char>& buffer);
		static bool readFileToBuffer(std::string filePath, std::string& buffer);

	};



#endif