#include <iostream>
#include "stdc++.h"
using namespace std; 
using std::endl;
using std::cout;
using std::cerr;
using std::ios;
using std::ifstream;
using std::ofstream;
using std::istringstream;
#include <fstream>
#include <sstream>
#include <array>
#include <ctime>
#include <climits>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <iterator>
#include <string>
#include <set>
//#include<bits/stdc++

typedef multimap<int,pair<double,vector <string> > > MI2D2S2S;
typedef map<double,vector <string> > D2VS;
typedef map<string,pair <int,double> > S2S2I;
typedef multimap<pair <int,double>, string > MI2S2I;
typedef map <string, vector <string> > MS2S;
typedef map<string,vector <double> > S2VD;


struct Cmdopts {
  string taxpairfile;
  string mashfile;
  string taxfile;
  string outprefix;
};

void helpmsg();
void getcmdopts(const int argc, char *argv[], Cmdopts &cmdopts);
void extractmash(const string mashfile, const string outprefix,set <string> &ids,S2VD &mashinfo);
void extracttax(const string taxpairfile, const string outprefix,set <string> & ids,S2VD &mashinfo);
void extractgcftax(const string taxfile, const string outprefix,set <string> & ids,S2VD &mashinfo);

int main(int argc, char *argv[]) {

  // read command line options
  Cmdopts cmdopts;
  getcmdopts(argc, argv, cmdopts);
  
  // extract and output genome sequences
  set <string> ids;
  S2VD mashinfo;
  extractmash(cmdopts.mashfile, cmdopts.outprefix,ids,mashinfo);
  extracttax(cmdopts.taxpairfile, cmdopts.outprefix,ids,mashinfo);
  extractgcftax(cmdopts.taxfile, cmdopts.outprefix,ids,mashinfo);
 
  return 0;
}

void extractmash(const string mashfile, const string outprefix,set <string> & ids,S2VD &mashinfo) {
  ifstream ifs(mashfile.c_str());
  if (!ifs) {
    cerr << "Could not open file " << mashfile << endl;
    exit(1);
  }
  //identity=float(items[0])
  //sharedhashes=items[1]
  //hashes=int(sharedhashes.split("/")[0])
  //medianmultiplicity=int(items[2])
  //pvalue=float(items[3])
  //queryid=items[4]
	  
  string line,acc,ass1,ass2,hashes;
  int pos1,pos2,pos3,totalhashes;
  double ani,sharedhashes;
  istringstream iss;
  ofstream ctgfile;
  while (getline(ifs, line)) {
	  
	  // Find position of '\t' using find() 
	  pos1 = line.find("\t"); 
	  ani = stod(line.substr(0,pos1));///////
	  cout << "ani:"<< ani<< "\t";
	  
	  pos2 = line.find("\t", pos1+1);
	  hashes =line.substr(pos1+1, pos2-pos1-1);////
	  cout << "hashes:"<< hashes<< "\t";
	  
	  pos3= line.find("/", pos1+1);
	  sharedhashes=stod(line.substr(pos1+1, pos3-pos1-1));////
	  cout << "sharedhashes:"<< sharedhashes<< "\t";
	  pos1 = line.find("\t", pos2+1);	
	  totalhashes=stod(line.substr(pos3+1,pos1-pos3-1));////
	  cout << "totalhashes:"<< totalhashes<< "\t";
	  	
	  pos1 = line.find("\t", pos1+1);			  	  

	  pos2 = line.find("\t", pos1+1);	
	  acc=line.substr(pos1+1, pos2-pos1-1);////
	  cout << "acc:"<< acc <<endl;
	  ids.insert(acc);
	  
	  vector <double> mash;
	  mash.push_back(ani);
	  mash.push_back(sharedhashes);
	  mash.push_back(totalhashes);
  
	  mashinfo.insert(pair<string,vector<double> >(acc, mash));
  }
}


// extract and output genome sequences
void extracttax(const string taxpairfile, const string outprefix,set <string> & ids,S2VD &mashinfo) {
  ifstream ifs(taxpairfile.c_str());
  if (!ifs) {
    cerr << "Could not open file " << taxpairfile << endl;
    exit(1);
  }
   //acc1=items[0]
  //acc2=items[1]
  //ass1=items[2]
  //ass2=items[3]
  //tax=int(items[4])
  //ani=float(items[5])		
  string line,acc1,acc2,ass1,ass2;
  int pos1,pos2,tax;
  double ani;
  istringstream iss;
  ofstream ctgfile;
  MI2D2S2S taxinfo;
  set<string>::iterator it; 
  set<string>::iterator it2; 
  
  while (getline(ifs, line)) {
	  // Find position of '\t' using find() 
	  pos1 = line.find("\t"); 
	  acc1 = line.substr(0,pos1);///////
	  cout << "acc1:"<< acc1<< "\t";
	  
	  pos2 = line.find("\t", pos1+1);
	  acc2 = line.substr(pos1+1, pos2-pos1-1);////
	  cout << "acc2:"<< acc2<< "\t";
	  
	  it = ids.find(acc1);
	  it2 = ids.find(acc2);
	  
	  if (it!=ids.end() && it2!=ids.end()){
		  cout << "found accs:"<< acc1<< "\t" <<acc2<<endl;
		  
		  pos1 = line.find("\t", pos2+1);		
		  ass1=line.substr(pos2+1,pos1-pos2 -1);////
		  cout << "ass1:"<< ass1 <<"\t";
	  
		  pos2 = line.find("\t", pos1+1);		
		  ass2=line.substr(pos1+1, pos2-pos1-1);////
		  cout << "ass2:"<< ass2 <<"\t";
	  
		  pos1 = line.find("\t", pos2+1);	
		  tax=stoi(line.substr(pos2+1,pos1-pos2 -1));////
		  cout << "tax:"<< tax <<"\t";
	  
		  pos2 = line.find("\t", pos1+1);		
		  ani=stod(line.substr(pos1+1, pos2-pos1-1));////
		  cout << "ani:"<< ani <<endl;
	  
		  vector <string> ass_acc;
		  ass_acc.push_back(acc1);
		  ass_acc.push_back(ass1);
		  ass_acc.push_back(acc2);
		  ass_acc.push_back(ass2);

		  pair <double, vector <string> > refani=make_pair(ani,ass_acc);
	  
		  taxinfo.insert(pair<int,pair<double,vector<string> > >(tax, refani));
		  ass_acc.clear();
		  
		  //check in pairs()
		  //pick highest ani pair if closely related
		  
	  }
	  else cout << "not found accs:"<< acc1<< "\t" <<acc2<<endl;
		  
  }
  cout << "FIN"<<endl;
  
  //typedef map<pair<double,vector <string> > D2VS;
  D2VS pair1;
  int prev_tax=0;
  set <string> max;
  set <string> min;
  string tempmax;
  string tempmin;
  double maxmash=0;
  double minmash=1;
  
  for (MI2D2S2S::iterator it= taxinfo.begin(); it != taxinfo.end(); ++it){
	  //vector of 4
	  //taxinfo.insert(pair<double,vector<string> > >refani));
	  if(it->first!=prev_tax && prev_tax!=0){
		  cout << "tax::"<< prev_tax<< endl;
		  for (D2VS::iterator it2= pair1.begin(); it2 != pair1.end(); ++it2){//start
			  cout << "for::it2->first:"<< it2->first<< endl;
			  S2VD::const_iterator it3 = mashinfo.find(it2->second[0]); 
			  if(it3!=mashinfo.end()){
				  cout << "mashscreen acc1:"<< it3->second[0]<< "\t"<< it3->second[1]<<"/"<< it3->second[2]<< "\t"<< it3->first<<endl;
				  maxmash=it3->second[0];
				  //minmash=0;
				  tempmax=it2->second[0];
				  tempmin=it2->second[2];
			  }
			  S2VD::const_iterator it4 = mashinfo.find(it2->second[2]); 
			  if(it4!=mashinfo.end()){
				  cout << "mashscreen acc2:"<< it4->second[0]<< "\t"<< it4->second[1]<< "/"<< it4->second[2]<< "\t"<< it4->first<< endl;
				  if (it4->second[0]>maxmash){
					  minmash=maxmash;
					  maxmash=it4->second[0];
				      tempmax=it2->second[2];
					  tempmin=it2->second[0];
				  }
				  else minmash=it4->second[0];
			  }
			  cout << "it2->second[0]:"<< it2->second[0]<< endl;
			  cout << "it2->second[1]:"<< it2->second[1]<< endl;
			  cout << "it2->second[2]:"<< it2->second[2]<< endl;
			  cout << "it2->second[3]:"<< it2->second[3]<< endl;
			  cout << "MAX:"<< maxmash<< "\t"<< tempmax<< endl;
			  cout << "MIN:"<< minmash<< "\t"<< tempmin<< endl;
			  
			  set<string >::iterator itmax; 
			  set<string >::iterator itmin; 
			  itmax = max.find(tempmin);
			  itmin = min.find(tempmax);
				  
			  if (itmin==min.end()){//max was NOT on min set, insert max
				  max.insert(tempmax);		  
			  }
			  if (itmax!=max.end()){//min was on max set, erese min
				  max.erase(tempmin);
				  //ids.erase(tempmin);
			  }
			  ids.erase(tempmin);			  
			  min.insert(tempmin);		  
		  }//fin
		  pair1.clear();
	  }
	  //else cout << "tax::else"<< it->first<< endl;

	  pair1.insert(pair<double,vector<string> >(it->second));
	  prev_tax=it->first;
  	
  }
  //////////////////////////////////////hererepeated from inside if(it->first!=prev_tax && prev_tax!=0){ 
  cout << "tax::"<< prev_tax<< endl;
  for (D2VS::iterator it2= pair1.begin(); it2 != pair1.end(); ++it2){//start
	  cout << "for::it2->first:"<< it2->first<< endl;
	  S2VD::const_iterator it3 = mashinfo.find(it2->second[0]); 
	  if(it3!=mashinfo.end()){
		  cout << "mashscreen acc1:"<< it3->second[0]<< "\t"<< it3->second[1]<<"/"<< it3->second[2]<< "\t"<< it3->first<<endl;
		  maxmash=it3->second[0];
		  //minmash=0;
		  tempmax=it2->second[0];
		  tempmin=it2->second[2];
	  }
	  S2VD::const_iterator it4 = mashinfo.find(it2->second[2]); 
	  if(it4!=mashinfo.end()){
		  cout << "mashscreen acc2:"<< it4->second[0]<< "\t"<< it4->second[1]<< "/"<< it4->second[2]<< "\t"<< it4->first<< endl;
		  if (it4->second[0]>maxmash){
			  minmash=maxmash;
			  maxmash=it4->second[0];
		      tempmax=it2->second[2];
			  tempmin=it2->second[0];
		  }
		  else minmash=it4->second[0];
	  }
	  cout << "it2->second[0]:"<< it2->second[0]<< endl;
	  cout << "it2->second[1]:"<< it2->second[1]<< endl;
	  cout << "it2->second[2]:"<< it2->second[2]<< endl;
	  cout << "it2->second[3]:"<< it2->second[3]<< endl;
	  cout << "OUT:MAX:"<< maxmash<< "\t"<< tempmax<< endl;
	  cout << "OUT:MIN:"<< minmash<< "\t"<< tempmin<< endl;
	  
	  set<string >::iterator itmax;
	  set<string >::iterator itmin; 
	  itmax = max.find(tempmin);
	  itmin = min.find(tempmax);
		  
	  if (itmin==min.end()){//max was NOT on min set, insert maxpair
		  max.insert(tempmax);		  
	  }
	  if (itmax!=max.end()){//min was on max set, erese minpair
		  max.erase(tempmin);
	  }
	  ids.erase(tempmin);
	  min.insert(tempmin);		  
  }//fin
  pair1.clear();

  
}
MI2S2I invertMap(S2S2I &map)
{
	MI2S2I multimap;

	for(S2S2I ::iterator it = map.begin(); it != map.end(); ++it){
		multimap.insert(pair<pair<int, double >,string >(it->second, it->first));
	}

	return multimap;
}
// extract and output genome sequences
void extractgcftax(const string taxfile, const string outprefix,set <string> & ids,S2VD &mashinfo) {
  ifstream ifs(taxfile.c_str());
  if (!ifs) {
    cerr << "Could not open file " << taxfile << endl;
    exit(1);
  }
  string line,acc,ass;
  int pos1,pos2,tax;
  double ani;
  istringstream iss;
  ofstream ctgfile;
  S2S2I gcftax;
  set<string>::iterator it1; 
  //multimap <string, string > MS2S;
  MS2S ass2acc;
  
  while (getline(ifs, line)) {
	  // Find position of '\t' using find() 
	  pos1 = line.find("\t"); 
	  ass = line.substr(0,pos1);
	  cout << "mc.tax ass:"<< ass<< ".\t";
	  pos2 = line.find("\t", pos1+1);
	  acc=line.substr(pos1+1,pos2-pos1 -1);////
	  cout << "mc.tax acc:"<< acc <<".\t";
	  it1 = ids.find(acc);
	  if (it1!=ids.end()){//start1	

		  S2VD::iterator it2= mashinfo.find(acc);
		  if (it2!=mashinfo.end()){
			  
			  pos1 = line.find("\t", pos2+1);	
			  tax=stoi(line.substr(pos2+1,pos1-pos2 -1));////
			  cout << " mc.tax tax:"<< tax << ".\t";
	  
			  MS2S::iterator it0= ass2acc.find(ass);
			  if (it0!=ass2acc.end()){
				  it0->second.push_back(acc);	
			  }
			  else{
				  vector <string> vs;
				  vs.push_back(acc);	  
				  ass2acc.insert(make_pair(ass, vs));	  
			  }
			  
			  //////////
			  ani=it2->second[0];
			  cout << " mc.tax ani:"<< ani << endl;
		  
			  S2S2I::iterator git = gcftax.find(ass);
			  if(git!=gcftax.end()){
				  if ( ani > (git->second).second )
					  (git->second).second=ani;
			  }
			  else{
				  pair <int, double > acc_tax= make_pair(tax,ani);
				  gcftax.insert(pair<string,pair<int,double> >(ass, acc_tax));	  
			  }
		  }// fi (it2!=mashinfo.end())
		  
	  }//end1
	  else cout << "not found acc:" << acc << endl;
		  
  }
  
  
  MI2S2I taxgcf = invertMap(gcftax);




/* 
  if (file.is_open()){
  	// Creating a iterator pointing to start of set
  	set<string >::iterator itset = ids.begin();
 
  	// Iterate till the end of set
  	while (itset != ids.end())
  	{
  		// Print the element
  		file << (*itset) << endl;
  		//Increment the iterator
  		itset++;
  	}
	
  	file.close();
  }	

*/
  
  MI2S2I::reverse_iterator it = taxgcf.rbegin();
  ofstream file (string(outprefix).c_str());
  if (file.is_open()){
    // Create a map iterator and point to the end of map
    int prev_tax=0,i=0;
    // Iterate over the map using Iterator till beginning.
    while (it != taxgcf.rend()) {
  	  if (prev_tax != it->first.first && prev_tax!=0){
  		  i=0;
  	  }		
  	  prev_tax=it->first.first;
  	  i++;
  	  if (i<=10){
  		  MS2S::iterator ita=ass2acc.find(it->second);
  		  if (ita!=ass2acc.end()){
  			  cout << "LAST:"<< it->first.first << "\t"<< it->first.second << "\t"<< it->second << "\t" << i <<endl;
			  
  			  for (vector<string>::const_iterator j = (ita->second).begin(); j != (ita->second).end(); ++j){
  			      cout << "LAST2:"<< *j << "\t";
		          file << *j << endl;
		      }
  			  cout << endl;

  		  }
  	  }
  	  it++;
    }
  	file.close();
  }	
}

// parse command line options
void getcmdopts(const int argc, char *argv[], Cmdopts &cmdopts) {

  if (argc != 5) {
    helpmsg();
    exit(1);
  }

  cmdopts.taxpairfile = argv[1];
  cmdopts.mashfile = argv[2];
  cmdopts.taxfile = argv[3];
  cmdopts.outprefix = argv[4];
}


// print out usage help message
void helpmsg() {
  cerr << endl;

  cerr << "Usage:" << endl;
  cerr << "        ./processmash <tax pair file> <mash screen file> <tax file> <output ids file>" << endl;
  cerr << endl;
  
  cerr << "Contact:" << endl;
  cerr << "        Have problems? Contact Victoria Cepeda - vcepeda@umiacs.umd.edu" << endl;
  cerr << endl;
}
