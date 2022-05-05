#ifndef _INTERACTORMDDRIVER_H_
#define _INTERACTORMDDRIVER_H_

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#include "Interactor.h"
#include "imd_interface.h"
#include "OPEP_wrapper.h"

#define FILTERSIZE 256
#define LOGFILENAMESIZE 256



static unsigned _receivedforces = 0;

class InteractorMDDriver : public Interactor
{
	public : 
		InteractorMDDriver();
		InteractorMDDriver(float imd_forcescale, 
											   int imd_wait,
											   int imd_port, 
											   int imd_debug);
		virtual ~InteractorMDDriver();
		void setPort(unsigned port);
		void setWait(unsigned wait);
		void setDebug(unsigned debug);
		void setLog(const char * logfilename);
		void setForceScale(float forcescale);
		
		virtual void start(int nbPos);
		virtual void stop();
		virtual void synchronize(double * pos, double * forces, int nbpos,
								 			 int nstep, double total_energy, 
								 			 double temperature, energies_t energies);
	
	protected : 
		char _IMDlogfilename[LOGFILENAMESIZE];
		int    _IMDdebug ;
		FILE * _IMDlog;
		int    _IMDstop  ;
		int _IMDmode ; 
		int _IMDwait ;
		int _IMDport ;
        float     _IMDforcescale ;

		IMDEnergies _IMDenergies;
		int _nbforces;
		int * _particleforceids; 
		float * _particleforces;
			
		pthread_t _imdthreads;
		pthread_mutex_t mutex;

		static int iterate( InteractorMDDriver * imdl);
		static void init( InteractorMDDriver * imdl);
		static void end( InteractorMDDriver * imdl);	

		static void  * runthread(void * userdata);
		virtual void init(int nbPos);
		
		
		
};

#endif