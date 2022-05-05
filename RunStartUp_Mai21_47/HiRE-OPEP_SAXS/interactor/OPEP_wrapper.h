#ifndef _OPEP_WRAPPER_H_
#define _OPEP_WRAPPER_H_

extern "C" {
	extern struct energies_t
	{
		double evdw,elec,eph,epa,ethh,etha,ebonh,ebona,
		       enbph,enbpa,eelph,eelpa,ehydro;

	} energies_;

	void interactor_start_(int * nbpos, 
						   double * imd_forcescale,
						   int * imd_wait,
						   int * imd_port,
						   int * imd_debug);

	void interactor_poll_(int * connected);

	void interactor_synchronize_(double * pos, double * forces, int * nbpos,
								 int * nstep, double * total_energy, double * temperature, double * ehb, double * estack, double * eahyd);
	void interactor_stop_();
}

#endif
