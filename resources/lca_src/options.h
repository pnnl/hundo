#pragma once
#include "libload.h"


std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);


struct options
{
public:
	options(int argc, char **argv, int);
	~options();
	const string TaxLvl2string();


	//vars
	string RefTaxFile, blastres, outF;
	string input_format;//bl8 , uc
	bool BLfilter;     //do my own blast filter before LCA
	bool calcHighMats; //calculate phylum etc level sum
	bool hitRD;
	bool isReads;
	bool annotateAll; //give out OTU / READ even if not assigned??
	bool nativeSlVdb;
	int numThr; // number of threads
	int taxDepth; //how deep does the taxonomy go?
	double LCAfract;
	vector<double> idThr;
	vector<string> blFiles, refDBs;
	vector<string> Taxlvls;
};

