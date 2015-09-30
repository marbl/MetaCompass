#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <fstream>
using std::ifstream;

#include <string>
using std::string;

#include <cstdlib>

void helpmsg();

int main(int argc, char *argv[]) {

	if (argc != 4) {
		helpmsg();
		exit(1);
	}
	
	ifstream ifs(argv[1]);
	if (!ifs) {
		cerr << "Could not open file: " << argv[1] << endl;
		exit(1);
	}

	unsigned int lencut = atoi(argv[2]), width = atoi(argv[3]);
	string eachline, seq, id;
	while (getline(ifs, eachline)) {
		if (eachline[0] == '>') { 
			if (seq.length() >= lencut) { 
				cout << id << endl;
				for (size_t i = 0; i < seq.length(); i += width) {
					cout << seq.substr(i, seq.length() - i < width ? seq.length() - i : width) << endl;
				}
			}
			id = eachline;
			seq.clear();
		}
		else seq += eachline;
	}

	if (seq.length() >= lencut) {
		cout << id << endl;
		for (size_t i = 0; i < seq.length(); i += width) {
			cout << seq.substr(i, seq.length() - i < width ? seq.length() - i : width) << endl;
		}
	}
}


void helpmsg() {
	cerr << endl;

	cerr << "Usage: " << endl;
	cerr << "         ./formatFASTA <fasta file> <min seq length> <# bases per line>" << endl;
	cerr << endl;

	cerr << "Contact: " << endl;
	cerr << "         Have problems? Contact Victoria Cepeda - vcepeda@umiacs.umd.edu" << endl;

	cerr << endl;
}
