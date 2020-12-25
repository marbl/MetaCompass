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
//#include<bits/stdc++.h> 
//TODO: check last element in all loops and couts

// An interval has start time and end time 
struct Interval 
{ 
	int start, end; 
};
string pathbin="";
string outprefix="";

typedef map<string,vector <string> > S2VS;//marker2cluster;
typedef multimap<string,Interval> S2Interval;//marker2interval;
typedef map<string,int> S2I;//marker2bases;marker2length;marker2numcog;depth2genome
typedef map<string,float> S2F;//marker2cov
typedef map<string,S2F> S2S2F;//marker2cov
typedef map<string, list<pair<string,float> > > S2LS2F;//genomecogcoverage;
typedef map<string,pair <float, float> > S2F2F;
typedef multimap<string,pair<string,vector <double> > > MS2VF;//genome2stats/pid-mismatches-evalue-score;
typedef map<string,vector <double> > S2VF;//genome2stats/pid-mismatches-evalue-score;
typedef map<string,double > S2D;//mean pid-mismatches-evalue-score;
typedef map<string,string > S2S;//genome2assembly
typedef multimap<string,string > MS2S;//genome2assembly
typedef multimap<pair<string,vector <double> > ,string > MVF2S;//genome2stats/pid-mismatches-evalue-score;

MVF2S invertMap(MS2VF &map)
{
	MVF2S multimap;

	for(MS2VF ::iterator it = map.begin(); it != map.end(); ++it){
		multimap.insert(pair<pair<string, vector <double> >,string >(it->second, it->first));
		
	}

	return multimap;
}

void helpmsg() {
  cerr << endl;

  cerr << "Usage:" << endl;
  cerr << "        ./breath <covthreshold> <outprefix> <pathbin>" << endl;
  cerr << endl;
  //cerr << "Parameters:" << endl;
       //<< "\t-h,--help\t\tShow this help message\n"
       //<< "\t-c,--coverage COVERAGE\tSpecify the minimum breadth of coverage"
       //<< "\t-o,--output-prefix OUTPUT\tSpecify the output directory"
   //          << std::endl;
  cerr << "Contact:" << endl;
  cerr << "        Have problems? Contact Victoria Cepeda - vcepeda@umiacs.umd.edu" << endl;
  cerr << endl;

}

void reverseStr(string& str) 
{ 
    int n = str.length(); 
  
    // Swap character starting from two corners 
    for (int i = 0; i < n / 2; i++) 
        std::swap(str[i], str[n - i - 1]); 
} 

double mean(vector <double> v){
	double mean = accumulate( v.begin(), v.end(), 0.0)/v.size(); 
	return mean;
}
// Function used in sort 
bool Intervalcomp(Interval a, Interval b) 
	{   return a.start > b.end; } 

/* Compares two intervals according to their staring time. 
This is needed for sorting the intervals using library 
function sort(). See http://goo.gl/iGspV */
bool compareInterval(Interval i1, Interval i2)
{ 
	return (i1.start < i2.start); 
} 

void read_clusters(S2VS &marker2cluster)
{
	string inputfile=outprefix+ "/contigs_clusters";
	ifstream file(inputfile);
	if (!file.is_open()) {
		cerr << "Could not open reference sequences file " << inputfile << endl;
		exit(1);
	}
	string line,representative,prev_representative;
	vector <string> clustervector;
	while (getline(file, line))
	{	// Process str
		size_t found=line.find('>');
		if (found!=string::npos){
			representative = line.erase(0,1);
			if ( !clustervector.empty()){
				marker2cluster.insert(pair<string,vector <string> >(prev_representative, clustervector));
				//cout<< "REPRESENTATIVE:\t"<<prev_representative << "\t"<<endl;
				//for (vector<string>::const_iterator it=clustervector.begin(); it!=clustervector.end(); ++it)
				//{cout <<"clustervector:"<<*it <<endl;}
				clustervector.clear();
			}
			prev_representative=representative;
		}
		else{
			clustervector.push_back(line);//add marker to cluster vector
		}				
	}
	//add last element
	if ( ! clustervector.empty() ){
		marker2cluster.insert(pair<string,vector <string> >(prev_representative, clustervector));
		/*cout<< "REPRESENTATIVE:\t"<<prev_representative << "\t"<<endl;
		for (vector<string>::const_iterator it=clustervector.begin(); it!=clustervector.end(); ++it) 
		{cout <<*it <<endl;}*/
		clustervector.clear();
	}
	for (S2VS::const_iterator it=marker2cluster.begin(); it!=marker2cluster.end(); ++it) {
		clustervector=it->second;
		cout <<"marker2cluster:"<<endl;
		cout <<">"<<it->first <<"\t"  <<clustervector.size() <<endl;
		for (vector<string>::const_iterator it2=clustervector.begin(); it2!=clustervector.end(); ++it2) 
		{cout <<*it2 <<endl;}
	   // for (unsigned i=0; i<clustervector.size(); ++i)
	     // {cout << clustervector[i] <<endl;}
	}
}

void read_blast_file(S2Interval &marker2interval, S2VS &marker2cluster, MS2VF &genome2stats,set <string> &keys)
{
	string inputfile=outprefix+ "/mc.blastn.all";
    //cout << "BLAST FILE:" << inputfile <<endl;
	ifstream file (inputfile);
	if (!file.is_open()) {
		cerr << "Could not open reference sequences file " << inputfile << endl;
		exit(1);
	}
	size_t  pos1 = 0, pos2 = 0, alglength, breadth_start, breadth_end, aux;
	string marker,line,genome,assembly;
	double pident,mismatch,evalue,score;
	vector <double> statsvector;
	//unordered 
	//set <string> keys;//set to store the unique keys
	while (getline(file, line))
	{
//////////////////
		// Process blast line
		pos1 = line.find("\t");           // 1 column read id, query
		pos2 = line.find("\t", pos1+1);   // 2 column marker id, subj
		marker = line.substr(pos1+1, pos2-pos1-1);
	    //cout << "marker:"<< marker<< "\t";
		pos1 = line.find("\t", pos2+1);   // 3 pident
		pident=stod(line.substr(pos2+1, pos1-pos2-1));
	    //cout << "pident:"<< pident<< "\t";
		//VECTOR 1
		statsvector.push_back(pident);
		pos2 = line.find("\t", pos1+1);   // 4 length
		alglength= stoi(line.substr(pos1+1, pos2-pos1-1));
	    //cout << "alglength:"<< alglength<< "\t";
		pos1 = line.find("\t", pos2+1);   // 5 mismatch
		mismatch=stod(line.substr(pos2+1, pos1-pos2-1));
	    //cout << "mismatch:"<< mismatch<< "\t";
		//VECTOR 2
		statsvector.push_back(mismatch);
		pos2 = line.find("\t", pos1+1);   // 6 gapopen
	    //cout << "gapoen:"<< line.substr(pos1+1, pos2-pos1-1)<< "\t";
		pos1 = line.find("\t", pos2+1);   // 7 qstart
	    //cout << "qstart:"<< line.substr(pos2+1, pos1-pos2-1)<< "\t";
		pos1 = line.find("\t", pos1+1);   // 8 qend
		pos2 = line.find("\t", pos1+1);   // 9 sstart start subj align
		breadth_start = stoi(line.substr(pos1+1, pos2-pos1-1));
	    //cout << "breadth_start:"<< breadth_start<< "\t";
		pos1 = line.find("\t", pos2+1);   // 10 send end subj align
		breadth_end = stoi(line.substr(pos2+1, pos1-pos2-1));
	    //cout << "breadth_end:"<< breadth_end<< "\t";
		
		if (breadth_end < breadth_start){
			aux=breadth_start;
			breadth_start=breadth_end;
			breadth_end=aux;
		}
		pos2 = line.find("\t", pos1+1);   // 11 evalue
		evalue=stod(line.substr(pos1+1, pos2-pos1-1));
	    //cout << "evalue:"<< evalue<< "\t";
		//VECTOR 3
		statsvector.push_back(evalue);
		pos1 = line.find("\t", pos2+1);   // 12 score
		score=stod(line.substr(pos2+1, pos1-pos2-1));
	    //cout << "score:"<< score<< endl;
		//VECTOR 4
		statsvector.push_back(score);
//////////////////
		//create interval
		Interval coordinates;
		coordinates.start=breadth_start;
		coordinates.end=breadth_end;
		marker2interval.insert(pair<string,Interval>(marker, coordinates));
		keys.insert(marker);
		//create genome stats
		pos1 = marker.find("_");//NC
		pos2 = marker.find("_", pos1+1);//[_0-9].[0-9]
		genome = marker.substr(0, pos2);
		pos1 = marker.find("_",pos2+1);//GCF
		pos1 = marker.find("_", pos1+1);//[_0-9].[0-9]
		assembly=marker.substr(pos2+1, pos1-pos2-1);
		pair <string, vector <double> > assemblystats;
		assemblystats = make_pair(assembly,statsvector);
		genome2stats.insert(pair<string,pair <string, vector <double> > >(genome, assemblystats));
	    cout << "genome2stats:"<< marker << "\t"<<genome<< "\t"<<assemblystats.first<< "\t";
		for (vector<double>::const_iterator it4=assemblystats.second.begin(); it4!=assemblystats.second.end(); ++it4){
			cout <<*it4 <<"\t";			
		}
		cout << endl;
		//if representative  of a cluster found in blast, add rest of members of cluster
		S2VS::const_iterator it2 = marker2cluster.find(marker); 
		if(it2!=marker2cluster.end()){
			vector<string> names = marker2cluster.at(marker);
		    cout << "representative genome2stats:"<< marker << "\t"<<genome<< "\t"<<assemblystats.first<<endl;
			for (vector<string>::const_iterator it3=names.begin(); it3!=names.end(); ++it3){
				//cout << *it3;
				string name=it3->c_str();
			    cout << "name:"<< name << "\t";
				marker2interval.insert(pair<string,Interval>(name, coordinates));
				keys.insert(name);
				//create genome stats
				pos1 = name.find("_");//NC
				pos2 = name.find("_", pos1+1);//[_0-9].[0-9]
				genome = name.substr(0, pos2);
				pos1 = name.find("_",pos2+1);//GCF
				pos1 = name.find("_", pos1+1);//[_0-9].[0-9]
				assembly=name.substr(pos2+1, pos1-pos2-1);
				pair <string, vector <double> > assemblystats;
				assemblystats = make_pair(assembly,statsvector);
				genome2stats.insert(pair<string,pair <string, vector <double> > >(genome, assemblystats));
			    cout << "genome2stats2:"<< name << "\t"<< genome<< "\t"<<assemblystats.first<< "\t";
				for (vector<double>::const_iterator it5=assemblystats.second.begin(); it5!=assemblystats.second.end(); ++it5){
					cout <<*it5 <<"\t";
				}
				cout << endl;
			}
		    cout << "END"<<endl;
		}
		statsvector.clear();	
	}
   // cout << "END BLAST FILE:" << inputfile <<endl;
}

void process_marker_ln(S2I & marker2length){
	string inputfile=pathbin+"/refseq/markers/markers.length";
	ifstream file(inputfile);
	if (!file.is_open()) {
		cerr << "Could not open reference sequences file " << inputfile << endl;
		exit(1);
	}
	size_t  pos1 = 0, length;
	string line,marker;	
	while (getline(file, line)){
		pos1 = line.find("\t");
		marker = line.substr(0, pos1);
		length = stoi(line.substr(pos1+1, line.size()-pos1-1));
		marker2length.insert(pair<string,int>(marker, length));			
	    //cout << "process_marker_ln:\t" << marker <<"\t"<< length << endl;
	}
	
}
void process_genome_marker_ln(S2I & genome_marker2length){
	string inputfile=pathbin+"/refseq/markers/genome2markers.length";
	ifstream file(inputfile);
	if (!file.is_open()) {
		cerr << "Could not open reference sequences file " << inputfile << endl;
		exit(1);
	}
	size_t  pos1 = 0, length;
	string line,genome;	
	while (getline(file, line)){
		pos1 = line.find("\t");
		genome = line.substr(0, pos1);
		length = stoi(line.substr(pos1+1, line.size()-pos1-1));
		genome_marker2length.insert(pair<string,int>(genome, length));			
	    //cout << "process_marker_ln:\t" << marker <<"\t"<< length << endl;
	}
	
}
void process_assembly_marker_ln(S2I & assembly_marker2length){
	string inputfile=pathbin+"/refseq/markers/assembly2markers.length";
	ifstream file(inputfile);
	if (!file.is_open()) {
		cerr << "Could not open reference sequences file " << inputfile << endl;
		exit(1);
	}
	size_t  pos1 = 0, length;
	string line,assembly;	
	while (getline(file, line)){
		pos1 = line.find("\t");
		assembly = line.substr(0, pos1);
		length = stoi(line.substr(pos1+1, line.size()-pos1-1));
		assembly_marker2length.insert(pair<string,int>(assembly, length));			
	    //cout << "process_marker_ln:\t" << marker <<"\t"<< length << endl;
	}
	
}

void process_marker_cog(S2I &marker2numcog){
	string inputfile=pathbin+"/refseq/markers/markers.count";	
	ifstream file(inputfile);
	if (!file.is_open()) {
		cerr << "Could not open reference sequences file " << inputfile << endl;
		exit(1);
	}
	size_t  pos1 = 0, cogs;
	string line,genome;	
	while (getline(file, line)){
		pos1 = line.find("\t"); 
		genome = line.substr(0, pos1);// 
		cogs = stoi(line.substr(pos1+1, line.size()-pos1-1));
		marker2numcog.insert(pair<string,int>(genome, cogs));
	    //cout << "process_marker_cog:\t" << genome <<"\t"<< cogs << endl;
	}
}

//https://www.geeksforgeeks.org/merging-intervals/
void mergeIntervals(Interval arr[], size_t n, int & index) 
{ 
	// Sort Intervals in decreasing order of 
	// start time 
	sort(arr, arr+n, Intervalcomp); 
	//cout << "\n Sorted Intervals are: "; 
	//for (int i = 0; i < n; i++) 
	//	cout << "[" << arr[i].start << ", " << arr[i].end << "] "; 
	//int index = 0; // Stores index of last element 
	// in output array (modified arr[])
	// Traverse all input Intervals 
	for (size_t i=0; i<n; i++) 
	{ 
		// If this is not first Interval and overlaps 
		// with the previous one 
		if (index != 0 && arr[index-1].start <= arr[i].end) 
		{ 
			while (index != 0 && arr[index-1].start <= arr[i].end) 
			{ 
				// Merge previous and current Intervals 
				arr[index-1].end = max(arr[index-1].end, arr[i].end); 
				arr[index-1].start = min(arr[index-1].start, arr[i].start); 
				index--; 
			} 
		} 
		else // Doesn't overlap with previous, add to 
			// solution 
			arr[index] = arr[i]; 

		index++; 
	} 

	// Now arr[0..index-1] stores the merged Intervals 
}

//per marker store alignlength, intervale.//after storing post-processing interval in stack, it's faster
//also calculate here marker_totalcount and total marker algn			
//outfile << marker + "\t" << alglength << "\t" << breadth_start << "\t" << breadth_end<< endl;
//
//create_interval_arrays(marker2interval,marker2intervals);
//taking too long
//key or marker
void create_interval_array(string key,S2Interval & marker2interval, size_t size,S2I & marker2bases){
	Interval *array=new Interval[size];//size of map
	//time_t nowtime;
	size_t i=0;
	int index1=0;
    cout << "key:"<<key << "\t"<< endl;
	S2Interval::iterator itr1=marker2interval.lower_bound(key);
	S2Interval::iterator itr2=marker2interval.upper_bound(key);
	while(itr1!=itr2){
		if(itr1->first==key){
			cout <<itr1->first <<"\t"<<itr1->second.start<<"\t"<<itr1->second.end<<endl;
			array[i].start=itr1->second.start;
			array[i].end=itr1->second.end;
			i++;
		}
		itr1++;
	}
    //nowtime = time(NULL);
    //cout << endl << "# mergeIntervals ... ";
	mergeIntervals(array, size, index1);
    //cout << "success!" << " (" << time(NULL) - nowtime << " seconds)" << endl;
	int total_bases=0;
	for (int i = 0; i < index1; i++){
		total_bases+=array[i].end-array[i].start+1;
	}
	marker2bases.insert(pair<string,int>(key, total_bases));
    cout << "bases! " << total_bases << endl;
	
}
void create_interval_arrays(S2Interval & marker2interval,S2I & marker2bases,set <string> &keys){
	//string prev="none";
	size_t size=0;
	time_t nowtime;
	//bounds:
	//replace by set keys
	for(set <string> ::iterator it = keys.begin(); it != keys.end(); ++it){
		string key=*it;
		size=marker2interval.count(key);
	    cout << "key:"<<key << "\tsize:"<<size<< endl;
		//new dec 2019
	    nowtime = time(NULL);
	    cout << "# create_interval_array ... "<< key<<endl;
		cout <<endl;
		
		create_interval_array(key,marker2interval,size,marker2bases);				
	    cout << "success!" << " (" << time(NULL) - nowtime << " seconds)" << endl;
		cout <<endl;
	}

}	
void marker_breadth_cov(S2F & marker2cov,S2F & genome2cov, S2F2F & assembly2cov, S2S & genome2assembly, S2I & marker2bases, S2I & marker2lenght, S2I & genome_marker2length, S2I & assembly_marker2length, float covthreshold,string outprefix){
	//calculate breadth of coverage per marker gene
	pair <string,float> cogcoverage;
	string genome,assembly, cog;
	float coverage,coverage2_num,coverage2_den;//breadth of coverage
	//typedef map<string,pair <float, float> > S2F2F	
	size_t pos1,pos2;
	for(S2I ::iterator it = marker2bases.begin(); it != marker2bases.end(); ++it){		
		S2I::iterator it2 = marker2lenght.find(it->first);//only one per marker, not multimap
		if(it2!=marker2lenght.end()){
		//process each marker separately
			pos1 = it->first.find("_");//NC
			pos2 = it->first.find("_", pos1+1);//[_0-9].[0-9]
			genome = it->first.substr(0, pos2);
			pos1 = it->first.find("_",pos2+1);//GCF
			pos1 = it->first.find("_", pos1+1);//[_0-9].[0-9]
			//pos2 = it->first.find("_", pos2+1);//[_0-9].[0-9]
			assembly=it->first.substr(pos2+1, pos1-pos2-1);
			pos2 = it->first.find("_", pos1+1);//[_0-9].[0-9]
			cog = it->first.substr(pos1+1, pos2-pos1-1);
			coverage=(float)it->second/it2->second;//bases/length
			coverage2_num=(float)it->second;//bases/length			
			coverage2_den=(float)it2->second;//bases/length	
		//end process each marker separately
			//insert marker2cov
			marker2cov.insert(pair<string,float>(it->first, coverage));
			//DOOOOOOOOO only if marker coverage is more than this!!!1/3 default
			if(coverage>=covthreshold){
			//insert genome maps:genome2assembly 
			genome2assembly.insert(pair<string,string>(genome, assembly));//good <3
			//no sirve de naaa 
			//cogcoverage = make_pair(cog,coverage);
			
			//create maps:genome2cov
			S2F::iterator it4 = genome2cov.find(genome+"\t"+assembly);
			if(it4!=genome2cov.end()){
				it4->second=it4->second + coverage2_num;
			}
			else{
				genome2cov.insert(pair<string,float>(genome+"\t"+assembly, coverage2_num));
			}			
			//end create maps:genome2cov

			//create maps:assembly2cov
			S2F2F::iterator it5 = assembly2cov.find(assembly);
			if(it5!=assembly2cov.end()){
				float num=(it5->second).first;
				(it5->second).first=num + coverage2_num;
			}
			else{
				pair <float,float> bases2length;
								
				S2I::iterator itlen = assembly_marker2length.find(assembly);
				if(itlen!=assembly_marker2length.end()){
					float total_ln=itlen->second;
					bases2length = make_pair(coverage2_num,total_ln);//coverage2_den);
					//if (coverage2_num/total_ln>=0.0){//potentially useful
					assembly2cov.insert(pair<string,pair<float,float> >(assembly, bases2length));
					//}
				}//endif
			}//endelse			
			//end create maps:assembly2cov
		}//DOOOOOOOOO
		}//end if
	}//end for
}//end void

int main(int argc, char *argv[])
{ 

	if (argc < 4) {
      helpmsg();
      exit(1);
    }

    float covthreshold = atof(argv[1]);
    outprefix = argv[2];
    pathbin = argv[3];
	
	cout << "PARAMETERS:\ncovthreshold: " << covthreshold << "\noutprefix: " << outprefix << "\npathbin: "<< pathbin <<endl;
	
	S2VS marker2cluster;
    time_t nowtime = time(NULL);
    cerr << endl << "# Reading and process clusters file ... ";
	read_clusters(marker2cluster);	
    cerr << "success!" << " (" << time(NULL) - nowtime << " seconds)" << endl;
	
	//PROCESS BLAST FILE
	S2Interval marker2interval;
	MS2VF genome2stats;
	set <string> keys;//set to store the unique keys
	
    nowtime = time(NULL);
    cerr << "# Read and process blast file ... ";
	read_blast_file(marker2interval,marker2cluster,genome2stats,keys);
    cerr << "success!" << " (" << time(NULL) - nowtime << " seconds)" << endl;
	
	S2I marker2bases;
    nowtime = time(NULL);
    cerr << "# create_interval_arrays ... ";
	create_interval_arrays(marker2interval,marker2bases,keys);
    cerr << "success!" << " (" << time(NULL) - nowtime << " seconds)" << endl;
	
	
	S2I marker2length;
	S2I marker2numcog;
	
	
    nowtime = time(NULL);
    cerr << "# process_marker_ln ... ";
	process_marker_ln(marker2length);
    cerr << "success!" << " (" << time(NULL) - nowtime << " seconds)" << endl;
	
    nowtime = time(NULL);
    cerr << "# process_marker_cog ... ";
	process_marker_cog(marker2numcog);
    cerr << "success!" << " (" << time(NULL) - nowtime << " seconds)" << endl;
	////////////////////////////////
	
	
	S2I genome_marker2length;
	S2I assembly_marker2length;
	process_genome_marker_ln(genome_marker2length);
	process_assembly_marker_ln(assembly_marker2length);
	
	S2F marker2cov;
	S2F genome2cov;
	S2F2F assembly2cov;
	S2LS2F genomecogcoverage;
	S2S genome2assembly;
    nowtime = time(NULL);
    cerr << "# marker_breadth_cov ... ";
	marker_breadth_cov(marker2cov,genome2cov,assembly2cov,genome2assembly,marker2bases,marker2length,genome_marker2length,assembly_marker2length,covthreshold,outprefix);
    cerr << "success!" << " (" << time(NULL) - nowtime << " seconds)" << endl;

//////////////////////create intervalfile marker_intervals. File contains all marker coordinates found in blastn results
	//write files
	S2I depth,depth2assembly;
	ofstream intervalfile (string(outprefix + "/marker_intervals").c_str());
	//typedef multimap<string,Interval> S2Interval;//marker2interval;
	string genome,assembly,prev_genome="none",prev_assembly="none";
	size_t pos1,pos2,depth_count = 0,assembly_count=0;
	if (intervalfile.is_open()){
		for(S2Interval ::iterator it = marker2interval.begin(); it != marker2interval.end(); ++it){
			intervalfile << it->first <<"\t"<< (it->second).start  <<"\t"<< (it->second).end <<endl;
			pos1 = (it->first).find("_");//NC
			pos2 = (it->first).find("_", pos1+1);//[_0-9].[0-9]
			genome = (it->first).substr(0, pos2);
			pos1 = it->first.find("_",pos2+1);//GCF
			pos1 = it->first.find("_", pos1+1);//[_0-9].[0-9]
			assembly=it->first.substr(pos2+1, pos1-pos2-1);
			
			if(genome.compare(prev_genome)!=0 && prev_genome.compare("none")!=0)
			{
				depth.insert(pair<string,int >(prev_genome, depth_count));
				depth_count=0;
			}
			if(assembly.compare(prev_assembly)!=0 && prev_assembly.compare("none")!=0)
			{
				depth2assembly.insert(pair<string,int >(prev_assembly, assembly_count));
				assembly_count=0;
			}
			prev_genome=genome;
			prev_assembly=assembly;
			assembly_count++;
			depth_count++;
		}
		//add last element
		depth.insert(pair<string,int >(prev_genome, depth_count));
		depth2assembly.insert(pair<string,int >(prev_assembly, assembly_count));
		
		intervalfile.close();
	}
	else cout << "Unable to open file" << outprefix << "/marker_intervals";

//////////////////////create breadth_marker_bases

	string file=outprefix + "/marker_bases";
	//ofstream myfile (string(outprefix + "/breadth_marker_bases").c_str());
	ofstream myfile (file);
	
	if (myfile.is_open()){
		for(S2I ::iterator it = marker2bases.begin(); it != marker2bases.end(); ++it){
			myfile << it->first <<"\t"<< it->second <<endl;
		}
		myfile.close();
	}
	else cout << "Unable to open file " << file <<endl;

//////////////////////create genomecogcoverage
	
	ofstream myfile2 (string(outprefix + "/genome2numcogs").c_str());
	S2I assembly2cog;
	if (myfile2.is_open()){//marker2cov sorted by accession, not assembly!!!!!!!!probably need to iterate again for assemblies.
		string prev_genome="none",prev_assembly="none",assembly_id,genome;
		size_t  pos1 = 0, pos2 = 0,i=0;
		for(S2F ::iterator it = marker2cov.begin(); it != marker2cov.end(); ++it){
			pos1 = it->first.find("_");//NC_
			pos2 = it->first.find("_", pos1+1);//NC_[_0-9].[0-9]_
			genome = it->first.substr(0, pos2);
			pos1 = it->first.find("_",pos2+1);//NC_[_0-9].[0-9]_GCF_
			pos1 = it->first.find("_", pos1+1);//NC_[_0-9].[0-9]_GCF_[_0-9].[0-9]_
			assembly_id=it->first.substr(pos2+1, pos1-pos2-1);
			if (prev_genome.compare(genome)!=0 && prev_genome.compare("none")!=0){
				myfile2 << prev_genome <<"\t"<<prev_assembly <<"\t"<<i<< endl;
				S2I::iterator it2 = assembly2cog.find(prev_assembly);
				if(it2!=assembly2cog.end()){
					it2->second+=i;
				}
				else{
					assembly2cog.insert(pair<string,int>(prev_assembly, i));
				}
				i=0;
			}
			
			prev_genome=genome;
			prev_assembly=assembly_id;
			i++;

		}
		myfile2 << prev_genome <<"\t"<<prev_assembly <<"\t"<<i<< endl;
		S2I::iterator it2 = assembly2cog.find(prev_assembly);
		if(it2!=assembly2cog.end()){
			it2->second+=i;
		}
		else{
			assembly2cog.insert(pair<string,int>(prev_assembly, i));
		}	
		myfile2.close();
	}
	else cout << "Unable to open file "<<outprefix <<"/genome2numcogs";
	
	ofstream myfile3 (string(outprefix + "/assembly2numcogs").c_str());
	if (myfile3.is_open()){//marker2cov sorted by accession, not assembly!!!!!!!!probably need to iterate again for assemblies.
		for(S2I ::iterator it = assembly2cog.begin(); it != assembly2cog.end(); ++it){
			myfile3 <<it->first<<"\t"<<it->second<< endl;
		}
		myfile3.close();
	}
	else cout << "Unable to open file "<<outprefix <<"/assembly2numcogs";
		
	ofstream statsfile (string(outprefix + "/genome_stats").c_str());
	//typedef multimap<string,Interval> S2Interval;//marker2interval;

	S2D pident;
	S2D mismatch;
	S2D evalue;
	S2D score;	
	vector <double> p,m,e,s;
	prev_genome="";prev_assembly="";
	if (statsfile.is_open()){
		//multimap<string,vector <float> > MS2VF
		for(MS2VF ::iterator it = genome2stats.begin(); it != genome2stats.end(); ++it){
			//statsfile << it->first <<"\t"<< it->second.first <<"\t"<< (it->second.second)[0]<<"\t"<< (it->second.second)[1]<<"\t"<< (it->second.second)[2]<<"\t"<< (it->second.second)[3] <<endl;
			if (prev_genome.compare(it->first)!=0 && prev_genome.compare("")!=0){
				statsfile << prev_genome <<"\t"<< prev_assembly <<"\t"<< mean(p) <<"\t"<< mean(m)<<"\t"<< mean(e)<<"\t"<< mean(s)<<endl;
				p.clear();m.clear();e.clear();s.clear();
			}
			p.push_back((it->second.second)[0]);
			m.push_back((it->second.second)[1]);
			e.push_back((it->second.second)[2]);
			s.push_back((it->second.second)[3]);
			prev_genome=it->first;
			prev_assembly=it->second.first;
		}
		statsfile << prev_genome <<"\t"<< prev_assembly <<"\t"<< mean(p) <<"\t"<< mean(m)<<"\t"<< mean(e)<<"\t"<< mean(s)<<endl;
		
		p.clear();m.clear();e.clear();s.clear();
		statsfile.close();
	}
	else cout << "Unable to open file genome_stats";
	
// invert the map
	MVF2S genome2stats2 = invertMap(genome2stats);
	
	ofstream statsfile2 (string(outprefix + "/assembly_stats").c_str());
	if (statsfile2.is_open()){
		prev_assembly="";
		p.clear();m.clear();e.clear();s.clear();

		for(MVF2S ::iterator it = genome2stats2.begin(); it != genome2stats2.end(); ++it){
			if (prev_assembly.compare(it->first.first)!=0 && prev_assembly.compare("")!=0){
				statsfile2 << prev_assembly <<"\t"<< mean(p) <<"\t"<< mean(m)<<"\t"<< mean(e)<<"\t"<< mean(s)<<endl;
				p.clear();m.clear();e.clear();s.clear();
			}		
			p.push_back((it->first.second)[0]);
			m.push_back((it->first.second)[1]);
			e.push_back((it->first.second)[2]);
			s.push_back((it->first.second)[3]);
			//prev=it->second;
			prev_assembly=it->first.first;	
		}
		statsfile2 << prev_assembly <<"\t"<< mean(p) <<"\t"<< mean(m)<<"\t"<< mean(e)<<"\t"<< mean(s)<<endl;
		p.clear();m.clear();e.clear();s.clear();
		statsfile2.close();
	}
	else cout << "Unable to open file assembly_stats";
//////////////////////create depth_genome_coverage

	ofstream depthfile (string(outprefix + "/depth_genome_coverage").c_str());
	if (depthfile.is_open()){
		for(S2I ::iterator it = depth.begin(); it != depth.end(); ++it){
			depthfile << it->first <<"\t"<< it->second <<endl;
		}
		depthfile.close();
	}
	else cout << "Unable to open file" <<outprefix << "/depth_genome_coverage";
	
//////////////////////create depth_assembly_coverage	
	ofstream depthfile2 (string(outprefix + "/depth_assembly_coverage").c_str());
	if (depthfile2.is_open()){
		for(S2I ::iterator it = depth2assembly.begin(); it != depth2assembly.end(); ++it){
			depthfile2 << it->first <<"\t"<< it->second <<endl;
		}
		depthfile2.close();
	}
	else cout << "Unable to open file" <<outprefix << "/depth_assembly_coverage";
	
//////////////////////create breadth_marker_coverage	
	
	ofstream myfile4 (string(outprefix + "/breadth_marker_coverage").c_str());
	
	if (myfile4.is_open()){
		for(S2F ::iterator it = marker2cov.begin(); it != marker2cov.end(); ++it){
			myfile4 << it->first <<"\t"<< it->second <<endl;
		}
		myfile4.close();
	}
	else cout << "Unable to open file breadth_marker_coverage";
//////////////////////create breadth_marker_coverage	

	ofstream myfile6 (string(outprefix + "/breadth_genome_bases_cov").c_str());
	if (myfile6.is_open()){
		for(S2F ::iterator it =  genome2cov.begin(); it != genome2cov.end(); ++it){
			myfile6 << it->first<<"\t"<< it->second <<endl;
		}
		myfile6.close();
	}
	else cout << "Unable to open file breadth_genome_bases_cov";
	
	

	/*for(S2F ::iterator it =  genome2cov.begin(); it != genome2cov.end(); ++it){
			float coverage2_num,coverage2_den;//breadth of coverage
			size_t pos1;
			pos1 = it->first.find("\t");
			pos2 = it->first.find("\t", pos1+1);
			assembly = it->first.substr(pos1+1, pos2-pos1-1);
			coverage2_num=it->second;
			
			S2F2F::iterator it2 = assembly2cov.find(assembly);
			if(it2!=assembly2cov.end()){
				float num=(it2->second).first;
				(it2->second).first=num+coverage2_num;
			}
			else{
				pair <float,float> bases2length;
								
				S2I::iterator itlen = assembly_marker2length.find(assembly);
				if(itlen!=assembly_marker2length.end()){
					int total_ln=itlen->second;
					bases2length = make_pair(coverage2_num,total_ln);//coverage2_den);
					//if (coverage2_num/total_ln>=0.0){//potentially useful
					assembly2cov.insert(pair<string,pair<float,float> >(assembly, bases2length));
					cout <<  assembly << endl;
					//}
				}//endif
			}//endelse			
	}//endfor
	
	*/
	
	
	
		
	ofstream file6 (string(outprefix + "/breadth_assembly_bases_cov").c_str());
	if (file6.is_open()){
		for(S2F2F ::iterator it =  assembly2cov.begin(); it != assembly2cov.end(); ++it){
			file6 << it->first<<"\t"<< (it->second).first <<"\t"<< (it->second).second<<"\t"<< (it->second).first/(it->second).second <<endl;
		}
		file6.close();
	}	
	else cout << "Unable to open file breadth_assembly_bases_cov";
	
	ofstream myfile5 (string(outprefix + "/mc.assembly.cov_cogs").c_str());
	if (myfile5.is_open()){
		for(S2F2F ::iterator it =  assembly2cov.begin(); it != assembly2cov.end(); ++it){
			if((it->second).first/(it->second).second >=0.3){
				S2I::iterator it2 = assembly2cog.find(it->first);
				if(it2!=assembly2cog.end()){
					if (it2->second >=10){
						myfile5 <<  it->first <<"\t"<< (it->second).first/(it->second).second<<"\t" <<it2->second<< endl;}
				}
			}
		}
		myfile5.close();
	}
	else cout << "Unable to open file "<<outprefix <<"/mc.assembly.ids";
	return 0;  	
}