#include "entities/molecule.hpp"
#include "entities/chain.hpp"
#include "entities/atom.hpp"
#include "entities/residue.hpp"
#include "geometry/vector3d.hpp"
#include "param/fa2cgParam.hpp"

#include <map>
#include <vector>
#include <algorithm> // for find

class Fa2cg
{
	public:
		Fa2cg();
		Molecule fa2cg(Molecule &molecule);

	private:
		void transformResidue(Residue &residue);
		void normaliseResidue(Residue &residue);
		Vector3d getCoordFromAtomVector(std::vector<Atom>);

};
