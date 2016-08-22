#include "options.h"

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

void self_help() {
    cout << "\n";
	cout << "Usage: LCA [options] -i [BLAST,...] -r [taxonomy,...] -o [output file]\n";
    cout << "\n";
	cout << "required:  -i PATH,...   BLAST results in m8 format\n";
	cout << "           -r PATH,...   reference taxonomy file with seq ID[tab]tax assignments;\n";
    cout << "                         seq ID must match BLAST hits\n";
	cout << "           -o PATH       output file containing the sequence name and the assigned\n";
    cout << "                         taxonomy against the ref database\n";
    cout << "\n";
	cout << "optional:  -X            calculate abundance of reads at different taxonomic levels\n";
	cout << "           -S            if a hit can be uniquely assigned to a single entry in the\n";
    cout << "                         ref database, this is reported in the output file\n";
	cout << "           -N            do not perform additional filtering\n";
	cout << "           -C            changes the tags attached to single reads (unsure of functionality)\n";
	cout << "           -L FLOAT      [0-1] the fraction of matching taxonomies required to accept\n";
    cout << "                         this taxonomy on the different levels (default: 0.8)\n";
	cout << "           -P LIST       comma separated list of min % identity to accept a database\n";
    cout << "                         hit as applicable to this taxonomic level, starting from\n";
    cout << "                         Species and going to Kingdom (default: 97,95,93,91,88,78,0)\n";
	cout << endl;
	exit(0);
}

options::options(int argc, char **argv,int defDep):
	RefTaxFile(""), blastres(""), outF(""), input_format("bl8"),
	BLfilter(true), calcHighMats(false), hitRD(false), isReads(false),
	annotateAll(false), nativeSlVdb(false),
	numThr(1), taxDepth(defDep), LCAfract(0.9f), idThr(defDep,0),
	blFiles(0), refDBs(0), Taxlvls(defDep)
{
	idThr[1] = 78; idThr[2] = 88; idThr[3] = 91; idThr[4] = 93;
	idThr[5] = 95; idThr[6] = 97;

	bool newIDthrs = false; string newIDthStr("");

	for (int i = 1; i < argc; i++)
	{
		if (!strcmp(argv[i], "-i"))
			blastres = argv[++i];
		else if (!strcmp(argv[i], "-h"))
			self_help();
		else if (!strcmp(argv[i], "-r"))
			RefTaxFile = argv[++i];
		else if (!strcmp(argv[i], "-o"))
			outF = argv[++i];
		else if (!strcmp(argv[i], "-X"))
			calcHighMats = true;
		else if (!strcmp(argv[i], "-S"))
			hitRD = true;
		else if (!strcmp(argv[i], "-N"))
			BLfilter = false;
		else if (!strcmp(argv[i], "-C"))
			isReads = true;
		else if (!strcmp(argv[i], "-L"))
			LCAfract = atof(argv[++i]);
		else if (!strcmp(argv[i], "-P")) {
			newIDthrs = true; newIDthStr= argv[++i];
		}

	}

	if (hitRD) {//needs to add 1 extra entry to some vectors
		/*taxDepth++;
		idThr.resize(taxDepth, 0);
		Taxlvls.resize(taxDepth,"Read")/*/
	}
	split(blastres, ',', blFiles);
	split(RefTaxFile, ',', refDBs);

	if (blFiles.size() != refDBs.size()) {
		cerr << "Unequal number of blast and refDB files!\n"; exit(24);
	}

	//check that all args were given
	bool isErr(false);
	if (blastres == "") { cerr << "Input file invalid or missing (-i)\n"; isErr = true; }
	if (RefTaxFile == "") { cerr << "RefDb file invalid or missing (-r)\n"; isErr = true; }
	if (outF == "") { cerr << "Output file invalid or missing (-o)\n"; isErr = true; }
	//if (blastres == "") { cerr << "Input file invalid (-f)"; isErr = true; }
	if (isErr) { self_help(); }

	if (newIDthrs) {//todo implement
		//cerr << "ID thresh inplement\n"; exit(55);
		vector<string> idthrsrev;
		split(newIDthStr, ',',idthrsrev);
		if (idthrsrev.size() != 7) {
			cerr << "Wrong number of id threshold levels (needs to be 7).\nAborting..\n"; exit(39);
		}
        cout << "Using tax cutoffs of: " << newIDthStr << "\n";
		for (size_t i = 0; i < 7; i++) {
			idThr[i] = atoi(idthrsrev[6 - i].c_str());
		}
	}

	vector<string> defTLvls(8, "");
	defTLvls[0] = "Domain"; defTLvls[1] = "Phylum"; defTLvls[2] = "Class";
	defTLvls[3] = "Order"; defTLvls[4] = "Family";  defTLvls[5] = "Genus";
	defTLvls[6] = "Species"; defTLvls[7] = "Strain";
	for (size_t i = 0; i < (size_t)taxDepth; i++) {
		if (i >= defTLvls.size()) { break; }
		Taxlvls[i] = defTLvls[i];
	}
}

const string options::TaxLvl2string() {
	string ret = Taxlvls[0];
	for (size_t i = 1; i < Taxlvls.size(); i++) {
		ret += __defaultTaxSep + Taxlvls[i];
	}
	return ret;
}


options::~options() {
}
