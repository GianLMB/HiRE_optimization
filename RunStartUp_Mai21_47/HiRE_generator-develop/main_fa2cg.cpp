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
    if(argc < 3 || std::string(argv[1]).find("-h") != std::string::npos  || std::string(argv[1]).find("--help") != std::string::npos)
    {
        std::cout << "Usage: " << argv[0] << " [AAfile] [CGfile]" << std::endl;
        return 0;
    }
    Molecule molecule = pdbParser::getMoleculeFromPdb(argv[1]);


    Fa2cg fa = Fa2cg();
    Molecule faMolecule = fa.fa2cg(molecule);
    PdbWriter::write(argv[2], faMolecule);

    return 0;
}
