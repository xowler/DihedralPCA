#include <iostream>
#include <string>

using namespace std;


#include "ezOptionParser.hpp"



int parseOptions(int argc, const char* argv[]){
    ez::ezOptionParser opt;
	opt.overview = "Calculates PCA on dihedral angles as given by g_rama";
	opt.syntax = "dihedralPCA [-s skip] [-n] [-t] prefix rama0.xvg [rama1.xvg ...]";

    ez::ezOptionValidator* vS4 = new ez::ezOptionValidator(ez::ezOptionValidator::S4);

	opt.add(
		"1", // Default.
		0, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Use every [skip] frames. default: 1", // Help description.
		"-s", // Flag token.
        vS4
	);

	opt.parse(argc, argv);

	if (opt.lastArgs.size() < 2) {
		cout << "ERROR: Expected (at least) 2 arguments, but got " 
             << opt.lastArgs.size() << ".\n\n";
		string usage;
		opt.getUsage(usage);
		cout << usage;
		return 1;
	} 



    

		
	std::cout << "First file: " << *opt.lastArgs[0] << std::endl;
	std::cout << "Second file: " << *opt.lastArgs[1] << std::endl;

	return 0;
}

int main(int argc, const char* argv[]){
    parseOptions(argc,argv);
}
