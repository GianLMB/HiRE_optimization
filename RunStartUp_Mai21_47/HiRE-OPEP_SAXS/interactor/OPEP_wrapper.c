#include <cstdlib>
#include <cstdio>

#include "OPEP_wrapper.h"
#include "InteractorMDDriver.h"

static InteractorMDDriver * imd = NULL;

extern "C" {

	void interactor_start_(int * nbpos, 
						   double * imd_forcescale,
						   int * imd_wait,
						   int * imd_port,
						   int * imd_debug)
	{
		printf("From OPEP_wrapper\n");
		if(imd == NULL)
			imd = new InteractorMDDriver(*imd_forcescale, 
										 *imd_wait, 
										 *imd_port, 
										 *imd_debug);

		if(imd)
		{
			imd->start(*nbpos);
		}
	}

  // Put true in connected if MDDriver is runnning
	void interactor_poll_(int * connected)
	{
		if(imd)
			*connected = imd->isActive();
    else
      *connected = 0;
	}

	void interactor_synchronize_(double * pos, double * forces, int * nbpos,
								 int * nstep, double * total_energy, double * temperature, double * ehb, double * estack, double * eahyd)
	{
		// printf("VDW: %f\n",energies_.etha);
		 energies_t energies = {*ehb,*estack,*eahyd, 0, 0, 0, 0, 0,
		          0, 0, 0, 0, 0};
		if(imd) {
			imd->synchronize(pos, forces, *nbpos,
							 *nstep, *total_energy, *temperature, energies);
						 }
	}

	void interactor_stop_()
	{
		if(imd)
			imd->stop();
	}
}
