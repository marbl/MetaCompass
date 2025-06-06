#ifndef FAI_PARSER_H
#define FAI_PARSER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>

class FAI_parser{
	std::string filename;
	std::ifstream fai_file;
	int buf_len;
	bool has_next_flag;

	void init();
	std::vector<std::string> get_tokenized_line();

public:
	FAI_parser(std::string filename, int buf_len = 1000);
	std::map<std::string, int> parse();
	void close();
};

#endif