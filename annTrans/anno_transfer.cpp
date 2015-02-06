#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <regex>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include <limits>
#include <cstdlib>

using namespace std;

struct featureLine{
	string seqname;
	string source;
	string type;
	int start;
	int end;
	string score;
	char strand;
	string phase;
	string attr;
};

struct _sortByStart {
        bool operator() (const featureLine& i, const featureLine& j){
                return i.start < j.start;
        }
} sortByStart;




void readFile(string filename, vector<featureLine>& lines); 

void overlapByRate(vector<vector<int> >& ret, vector<featureLine>& query, vector<featureLine>& subject, double queryoverlap, double subjectoverlap);
void anyOverlap(vector<vector<int> >& ret, vector<featureLine>& query, vector<featureLine>& subject) {
  stable_sort(subject.begin(), subject.end(), sortByStart);
  for (vector<featureLine>::iterator i = query.begin();i != query.end();i++) {
	 ret.push_back(vector<int>());
	for (vector<featureLine>::iterator j = subject.begin();j != subject.end();j++) {
		 if (j->start > i->end) break;
		 int indexj = j - subject.begin();
		 if (j->end >= i->start) {
			ret[ret.size() - 1].push_back(indexj);
		}
	 }
  }
}

void readFasta(string filename, map<string, string>& output);
void readTrinotate(string filename, map<string, map<string, string> >& output);

void listTrinotate(ostream& file, string id, int start, int end, map<string, string>& fastamap, string seqname);

void writeOutput(string filename, vector<featureLine>& output);

inline string string_match(string str, string pattern) {
	string start = pattern + "=";
	string end = ";";
	size_t spos = str.find(start) + start.size();
	size_t epos = str.find(end,spos);
	return str.substr(spos, epos - spos);
}

vector<string> goodseq;

int main(int argc, char* argv[]) {
  //if (argc != 4) {
  //	exit(1);
  //}
  
  clock_t timer = clock();
  ios_base::sync_with_stdio(false);
  
  bool force_unsupported_features = false;
  
  string targetFile = "";
  string sourceFile = "";
  string seqFile = "";
  string outFile = "";
  double targetOverlap = -1.0;
  double sourceOverlap = -1.0;
  bool in_test_mode = false;
  bool targetPref = false;
  while (1) {
	 int option_index = 0;
	 struct option options[] = {
		{"target-pref", no_argument, NULL, 0},
		{"target", required_argument, NULL, 1},
		{"source", required_argument, NULL, 2},
		{"genome-sequence", required_argument, NULL, 3},
		{"output", required_argument, NULL, 4},
		{"force-unsupported-features", no_argument, NULL, 5},
		{"target-overlap", required_argument, NULL, 6},
		{"source-overlap", required_argument, NULL, 7},
		{"test", no_argument, NULL, 8},
		{0,0,0,0}
	 };
	 int c = getopt_long(argc, argv, "", options, &option_index);
	 if (c == -1) break;
	 switch (c) {
		case 0:
		  targetPref = true;
		  break;
		case 1:
		  targetFile = optarg;
		  break;
		case 2:
		  sourceFile = optarg;
		  break;
		case 3:
		  seqFile = optarg;
		  break;
		case 4:
		  outFile = optarg;
		  break;
		case 5:
		  force_unsupported_features = true;
		  break;
		case 6:
		  targetOverlap = atof(optarg);
		  break;
		case 7:
		  sourceOverlap = atof(optarg);
		  break;  
		case 8:
		  in_test_mode = true;
		  break;
		default:
		  break;	
	 }
  }
  //cout << targetFile << "\n" << sourceFile << "\n" << seqFile << "\n" << outFile << "\n" << targetPref << "\n" << force_unsupported_features << "\n" << targetOverlap << "\n" << sourceOverlap << endl;
  //string targetFile = "";
  //  string sourceFile = "";
  //    string seqFile = "";
  //      string outFile = "";
  //
  if (targetFile == "" || sourceFile == "" || seqFile == "" || outFile == "" || sourceOverlap == -1.0 || targetOverlap == -1.0) {
	cerr << "Argument incomplete: please see usage" << endl;
	return 1;
  }
  vector<string> features;
  features.push_back("gene"); 
  vector<featureLine> targetLines, sourceLines;
  // If not in test mode, the duplicated ids are fixed for trinotate
  if (!in_test_mode) {
	system(("./fix_id.sh " + targetFile).c_str());
	targetFile = "fixed.gff3";
  }
  readFile(targetFile, targetLines);
  readFile(sourceFile, sourceLines);
  
  {	
	 ifstream seqlist("R64PlusDBL.goodcontigs");
	 string str;
	 while (seqlist.good()) {
		str.clear();
		getline(seqlist, str, '\t');
		goodseq.push_back(str);
		seqlist.ignore(256, '\n');
	 }
	 seqlist.close();
  }
  
  vector<featureLine> output;
  map<string, int> namemap;
  map<string, string> fastamap;
  readFasta(seqFile, fastamap);
  ofstream tempSeq("tempSeq.fasta");
  
  //cout << targetLines[3].attr << "\n";
  //cout << string_match(targetLines[3].attr, "ID") << endl;
  //readFile(argv[1], targetLines);
  //readFile(argv[2], sourceLines);
  vector<vector<int> > pairs;
  vector<int> unsupported;
  if (targetPref) sourceOverlap = -1.0;
  if (force_unsupported_features) {
	 anyOverlap(pairs, targetLines, sourceLines);
	for (int i = 0;i < sourceLines.size();i++) unsupported.push_back(i);
  }
  else {
	 overlapByRate(pairs, targetLines, sourceLines, targetOverlap, sourceOverlap);
  }
  
  vector<int> listLines;
  int lastncRNA = -1;
  bool transferred = false;
  vector<string> dupLog;
  vector<string> transLog;
  cout << "Time spent for preparation: " << (clock() - timer) / CLOCKS_PER_SEC << "\n";
  for (vector<vector<int> >::iterator i = pairs.begin();i != pairs.end();i++) {
	 bool annotated = false;
	 int indexi = i - pairs.begin();
	 transferred = false;
	 cout << "Processing line " << indexi << ": " << targetLines[indexi].seqname << ", " << targetLines[indexi].start << "-" << targetLines[indexi].end << ", " << targetLines[indexi].attr << "\n";
	 if (!i->empty()/* && find(goodseq.begin(), goodseq.end(), targetLines[indexi].seqname) != goodseq.end()*/) {
		for (vector<int>::iterator j = i->begin();j != i->end();j++) {
		  int indexj = *j;
		  //cout << "i:" << targetLines[indexi].start << "-" << targetLines[indexi].end << ", j:" << sourceLines[indexj].start << "-" << sourceLines[indexj].end << "\n";
		  /*if (targetLines[indexi].start == 11202 && sourceLines[indexj].start == 11202) {
			c *out << targetLines[indexi].seqname << "\n" << sourceLines[indexj].seqname << "\n";
			cout << features[0] == sourceLines[indexj].type;
		}*/
		  if (targetLines[indexi].seqname == sourceLines[indexj].seqname && find(features.begin(), features.end(), sourceLines[indexj].type) != features.end()) {
			 annotated = true;
			 // In force unsupported feature mode, remove the supported from the list
			 if (force_unsupported_features) {
				vector<int>::iterator supported = find(unsupported.begin(), unsupported.end(), indexj);
				if (supported != unsupported.end()) unsupported.erase(supported, supported + 1);
			}
			 if (transferred && targetPref) {
				string id = string_match(sourceLines[indexj].attr, "ID");
				dupLog.push_back(id);
				continue;
			 }
			 if (targetPref) {
				string id = string_match(sourceLines[indexj].attr, "ID");
				
				transLog.push_back(id);
			 }
			 if (targetLines[indexi].strand == sourceLines[indexj].strand) {
				if (find(listLines.begin(), listLines.end(), indexj) != listLines.end()) {
				  cout << "Line dropped, source already transferred: " << sourceLines[indexj].start << "-" << sourceLines[indexj].end << "\n";
				  continue;
				}
				output.push_back(sourceLines[indexj]);
				listLines.push_back(indexj);
				string name = string_match(sourceLines[indexj].attr, "Name");
				cout << "Transferred annotation from range " << sourceLines[indexj].start << "-" << sourceLines[indexj].end << "(" << name << "), replacing the target feature" << "\n";
				
				unsigned int indexk = indexj + 1;
				for (string parent = string_match(sourceLines[indexk].attr, "Parent"); indexk < sourceLines.size() && (parent = string_match(sourceLines[indexk].attr, "Parent")) == name;indexk++) {
				  cout << "Transfer line " << indexk << " from source: child of " << name << "\n";
				  output.push_back(sourceLines[indexk]);
				  listLines.push_back(indexk);
				  
				}
			 }
			 else {
				//if (indexi != lastncRNA) {
				if (!transferred) {
				  lastncRNA = indexi;
				  cout << "Range " << targetLines[indexi].start << "-" << targetLines[indexi].end << " marked as ncRNA from range " << sourceLines[indexj].start << "-" << sourceLines[indexj].end << "\n";
				  string id = string_match(sourceLines[indexj].attr, "ID");

				  featureLine line = targetLines[indexi];
				  line.attr += ";gene=" + id;
				  line.type = "ncRNA";
				  output.push_back(line);
				  line.type = "exon";
				  output.push_back(line);
				}
			 }
			 transferred = true;
		  }
		}
	 }
	 if (!annotated) {
		cout << "No existing annotation found at range " << targetLines[indexi].start << "-" << targetLines[indexi].end  << ", " << targetLines[indexi].attr << ", added to trinotate list\n";
		string id = string_match(targetLines[indexi].attr, "ID");
		listTrinotate(tempSeq, id, targetLines[indexi].start, targetLines[indexi].end, fastamap, targetLines[indexi].seqname);
		featureLine line = targetLines[indexi];
		output.push_back(line);
		namemap.insert(pair<string, int>(id, output.size() - 1));
		//cout << "namemap ID: " << id << "\n" << output[namemap[id]].end << "\n";
		line.type = "exon";
		output.push_back(line);
	 }
  }
  if (force_unsupported_features) {
	cout << "Transferring unsupported features" << "\n";
	for (int i = 0;i < unsupported.size();i++) {
		cout << sourceLines[unsupported[i]].start << "-" << sourceLines[unsupported[i]].end << ":" << sourceLines[unsupported[i]].type << "\n";
		output.push_back(sourceLines[unsupported[i]]);
	}
  }
  if (!in_test_mode) {
	system("./trinotAdd.sh tempSeq.fasta Trinotate.sqlite");
  }
  cout << "importing trinotate report" << "\n";
  map<string, map<string, string> > trinotate;
  readTrinotate("trinotate_annotation_report.xls", trinotate);
  for (map<string, map<string, string> >::iterator i = trinotate.begin();i != trinotate.end();i++) {
	 //cout << i->first << "\n";
	 //cout << namemap[i->first] << "\n";
	 //cout << output[namemap[i->first]].end << "\n";
	 if (namemap.find(i->first) == namemap.end()) continue;
	 for (map<string, string>::iterator j = i->second.begin();j != i->second.end();j++) {
		if (j->second != ".") {
		  output[namemap[i->first]].attr = output[namemap[i->first]].attr + " ; " + j->first + " : " + j->second;
		  output[namemap[i->first] + 1].attr = output[namemap[i->first] + 1].attr + ";" + j->first + ":" + j->second;
		}
	 }
  }
  //writeOutput(outFile, output);
  //map<string, map<string, string> > trinotate;
  //readTrinotate("trinotate_annotation_report.xls", trinotate);
  //return 0;
  writeOutput(outFile, output);
  //writeOutput("source", sourceLines);
  tempSeq.close();
  cout << "Time spent: " << (clock() - timer) / CLOCKS_PER_SEC << endl;
}

void readFile(string filename, vector<featureLine>& lines) {
	ifstream file(filename.c_str());
	string line;
	while (file.good()) {
		getline(file, line);
		//cout << line << "\n";
		if (line.empty()) continue;
		if (line[0] != '#') {
			featureLine feature;
			stringstream ss(line);
			getline(ss, feature.seqname, '\t');
			//cout << feature.seqname << "\n";
			getline(ss, feature.source, '\t');
			//cout << feature.source << "\n";
			getline(ss, feature.type, '\t');
			ss >> feature.start;
			//cout << feature.start << "\n";
			ss >> feature.end;
			ss.ignore(256, '\t');
			//cout << feature.end << "\n";
			getline(ss, feature.score, '\t');
			//cout << feature.score << "\n";
			ss >> feature.strand;
			ss.ignore(256,'\t');
			//cout << feature.strand << "\n";
			getline(ss, feature.phase, '\t');
			//cout << feature.phase << "\n";
			getline(ss, feature.attr);
			//cout << feature.attr << "\n";
			if (ss.bad() || ss.fail()) {
				cout << "Bad format from file " << filename << endl;
				return;
			}
			else {
				lines.push_back(feature);
			}
		}
	}
	file.close();	
}
void overlapByRate(vector<vector<int> >& ret, vector<featureLine>& query, vector<featureLine>& subject, double queryoverlap, double     subjectoverlap) {


	sort(subject.begin(), subject.end(), sortByStart);
	for (vector<featureLine>::const_iterator i = query.cbegin();i != query.cend();i++) {
		ret.push_back(vector<int>());
		int qoverlap = queryoverlap >= 0 ? (i->end - i->start) * queryoverlap : numeric_limits<int>::max();
		bool test = false;
		if (i->start == 132144) {
			cout << "Required overlap: " << qoverlap << "\n";
			test = true;
		}
		for (vector<featureLine>::const_iterator j = subject.cbegin();j != subject.cend();j++) {
			//cout << "i:" << i->start << "-" << i->end << ", j:" << j->start << "-" << j->end << "\n";

			int soverlap = subjectoverlap >= 0 ? (j->end - j->start) * subjectoverlap : numeric_limits<int>::max();
			if (j->start > i->end) break;
			if (j->end <= i->end) {
				if (j->end - max(i->start, j->start) >= min(qoverlap, soverlap)) {
					ret[ret.size() - 1].push_back(j - subject.begin()); 
					
					//cout << "included" << "\n";
					//cout << "i:" << i->start << "-" << i->end << ", j:" << j->start << "-" << j->end << "\n";
					}
			}
			else if (i->end - max(j->start, i->start) >= min(qoverlap, soverlap)){
				ret[ret.size() - 1].push_back(j - subject.begin());
				

				//cout << "i:" << i->start << "-" << i->end << ", j:" << j->start << "-" << j->end << "\n";
			}
		}
	}
}
 
void readFasta(string filename, map<string, string>& output) {
	ifstream file(filename.c_str());
	string name, seq, temp;
	while (file.good()) {
		temp.clear();
		getline(file, temp);
		if (temp.empty()) {
			output.insert(pair<string, string> (name, seq));
			name =  "";
			seq = "";
			continue;
		}
		//cout << temp << "\n";
		if (temp[0] == '>') {
			output.insert(pair<string, string> (name, seq));
			name = temp.substr(1, temp.find_first_of(" \t\n") - 1);
			seq = "";
		}
		else {
			seq += temp;
		}
	}
	
	output.insert(pair<string, string> (name, seq));
	output.erase("");
	if (file.bad()) {
		cout << "Error reading fasta file " << filename << endl;
	}
	file.close();
}
void readTrinotate(string filename, map<string, map<string, string> >& output) {
	ifstream file(filename.c_str());
	string str, id, line;
	// We always skip the first column;
	vector<string> attr;
	file.ignore(256, '\t');
	// second col is the key field
	file.ignore(256, '\t');
	getline(file, line);
	stringstream ss(line);
	while (ss.good()) {
		getline(ss, str, '\t');
		attr.push_back(str);
	}
	while (file.good()) {
		line.clear();
		getline(file,line);
		//cout << line << "\n";
		ss.clear();
		ss.str(line);
		ss.ignore(256, '\t');
		id.clear();
		getline(ss, id, '\t');
		//cout << id << endl;
		map<string, string> values;
		for (vector<string>::iterator i = attr.begin();i != attr.end();i++) {
			str.clear();
			getline(ss, str, '\t');
		//	cout << *i << endl;
		//	cout << str << endl;
			values.insert(pair<string, string>(*i, str));
		}
		output.insert(pair<string, map<string, string> >(id, values));
	}
	if (file.bad()) {
		cout << "Error reading trinotate report from file " << filename << endl;
	}
	file.clear();
	file.close();
}
 
void listTrinotate(ostream& file, string id, int start, int end, map<string, string>& fastamap, string seqname) {
	file << ">" << id;
	string str(fastamap[seqname].substr(start, end - start + 1));
	unsigned int i;
	for (i = 0;i < str.size();i++) {
		if (i % 70 == 0) file << "\n";
		file << str[i];
	}
	file << "\n";
}
 
struct _isSameRep {
        bool operator() (const featureLine& i) {
                static featureLine j = i;
		if (j.seqname != i.seqname) {
			j = i;
			return true;
		}
                return false;
        }
} isSameRep;


void inline sortBySeqname(vector<featureLine>& output) {
	for (vector<featureLine>::iterator i = output.begin();i != output.end();) {
                vector<featureLine>::iterator j = find_if(i, output.end(), isSameRep);
                stable_sort(i, j, sortByStart);
                i = j;
        }

}

void writeOutput(string filename, vector<featureLine>& output) {
	
	ofstream file((filename + ".gff").c_str());

	sortBySeqname(output);
	//sort(output.begin(), output.end(), comp);
	for (vector<featureLine>::const_iterator i = output.cbegin();i != output.cend();i++) {
		file << i->seqname << "\t" << i->source << "\t" << i->type << "\t" << i->start
		<< "\t" << i->end << "\t" << i->score << "\t" << i->strand << "\t" << i->phase
		<< "\t" << i->attr << "\n";
	}
	file.clear();
	file.close();
}

