#include <iostream>
#include <vector>

#include "fa2cg.hpp"

#include "entities/molecule.hpp"
#include "inout/pdbParser.hpp"
#include "inout/opepInputParamWriter.h"
#include "inout/pdbWriter.hpp"
#include "inout/xtcParser.hpp"
#include "geometry/analysis.hpp"
#include "geometry/geometry.hpp"
#include "geometry/vector3d.hpp"
#include "mathematics/stats.hpp"
#include "param/fa2cgParam.hpp"

using namespace std;


int main(int argc, char ** argv)
{
    if(argc < 2 || std::string(argv[1]).find("-h") != std::string::npos  || std::string(argv[1]).find("--help") != std::string::npos)
    {
        std::cout << "A PDB file must be passed as an argument." << std::endl;
	std::cout << "Call: " << argv[0] << "[--circular] [PDBfile]" << std::endl;
        return 0;
    }

	bool circular = false;
	//trjconv -f min.xtc -o traj.pdb -s conf_initiale_RNA.pdb
	string fileName = argv[1];
	if(fileName.find("--circular") != std::string::npos)
	{
	    circular = true;
	    if(argc < 3)
	    {
		std::cout << "A PDB file must be passed as an argument." << std::endl;
		return 0;
	    }
	    fileName = argv[2];
	}

	Molecule molecule = pdbParser::getMoleculeFromPdb(fileName);
	//PdbWriter::write("test.pdb", molecule);
	//Analysis analysis(molecule);
	//analysis.plotAntiSyn();
	//analysis.plotBasesEcarts();

	string name = fileName.substr (0, fileName.size()-4);

	Fa2cg fa = Fa2cg();
	Molecule faMolecule = fa.fa2cg(molecule);
	PdbWriter::write(name+"_CG.pdb", faMolecule);
	PdbWriter::write("conf_initiale_RNA.pdb", faMolecule);
	OpepInputParamWriter opepInputParamWriter(faMolecule);

	opepInputParamWriter.write(name+".top", circular);
	//PdbWriter::write("test.pdb", molecule);


	//string fileName = argv[1];

	//Molecule molecule = pdbParser::getMoleculeFromPdb(fileName);
	//PdbWriter::write("test.pdb", molecule);
	//Analysis analysis(molecule);
	//analysis.plotStrongStopDist();


	return 0;
}
