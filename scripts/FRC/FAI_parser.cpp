#include "FAI_parser.h"

using namespace std;

void FAI_parser::init(){
	fai_file.open(filename, ifstream::in);
	//TODO have this actually throw something useful and output it
	if(!fai_file){
		cout << "Could not open FAI file: " << filename << "\n";
		throw 1; 
	}
}

FAI_parser::FAI_parser(string filename, int buf_len) : 
	filename(filename), 
	buf_len(buf_len),
	has_next_flag(true) 
{ init(); }

vector<string> FAI_parser::get_tokenized_line(){
	string buf;
	if(getline(fai_file, buf)){
		vector<string> tokens;
		istringstream ss(buf);
		string temp;
		while(getline(ss, temp, '\t')){
			tokens.push_back(temp);
		}
		return tokens;
	}
	has_next_flag = false;
	return vector<string>();
}

map<string, int> FAI_parser::parse(){
	map<string, int> clens;

	while(has_next_flag){
		vector<string> data = get_tokenized_line();
		if(!has_next_flag) return clens;

		string cname(data[0]);
		int clen = stoi(data[1]);
		clens[cname] = clen;	
	}
	return clens;
}

void FAI_parser::close(){
	fai_file.close();
}