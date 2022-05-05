#include "InteractorMDDriver.h"

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>



using namespace std;

InteractorMDDriver::InteractorMDDriver() : Interactor()
{
	_IMDdebug  = 0;
	strcpy(_IMDlogfilename,"test.log");
	_IMDlog = NULL;
	_IMDstop   = 0;
	_IMDwait           = 1;
	_IMDport           = 3000;
	_IMDmode = 1; //server
	_IMDforcescale =1.0;
	_IMDenergies.tstep  = 0.0;
	_IMDenergies.T      = 0.0;
	_IMDenergies.Etot   = 0.0;
	_IMDenergies.Epot   = 0.0;
	_IMDenergies.Evdw   = 0.0;
	_IMDenergies.Eelec  = 0.0;
	_IMDenergies.Ebond  = 0.0;
	_IMDenergies.Eangle = 0.0;
	_IMDenergies.Edihe  = 0.0;
	_IMDenergies.Eimpr  = 0.0;
	
	pthread_mutex_t mutextmp1 = PTHREAD_MUTEX_INITIALIZER;
	mutex=mutextmp1;
	
}

InteractorMDDriver::InteractorMDDriver(float imd_forcescale = 1.0, 
									   int imd_wait = 0,
									   int imd_port = 3000,  
									   int imd_debug = 0) : Interactor()
{
	_IMDdebug  = imd_debug;
	strcpy(_IMDlogfilename,"test.log");
	_IMDlog = NULL;
	_IMDstop   = 0;
	_IMDwait           = imd_wait;
	_IMDport           = imd_port;
	_IMDmode = 1; //server
	_IMDforcescale = imd_forcescale;
	_IMDenergies.tstep  = 0.0;
	_IMDenergies.T      = 0.0;
	_IMDenergies.Etot   = 0.0;
	_IMDenergies.Epot   = 0.0;
	_IMDenergies.Evdw   = 0.0;
	_IMDenergies.Eelec  = 0.0;
	_IMDenergies.Ebond  = 0.0;
	_IMDenergies.Eangle = 0.0;
	_IMDenergies.Edihe  = 0.0;
	_IMDenergies.Eimpr  = 0.0;
	
	pthread_mutex_t mutextmp1 = PTHREAD_MUTEX_INITIALIZER;
	mutex=mutextmp1;
	
}

InteractorMDDriver::~InteractorMDDriver()
{
}


void InteractorMDDriver::init(int nbPos)
{
	Interactor::init(nbPos);
	// //_IMDlog =  IIMD_init( "", &(_IMDmode),&(_IMDwait),&(_IMDport), &(_IMDdebug),0 );
	// //cout<<"_IMDforcescale "<<_IMDforcescale<<" kcal.mol-1.A-1"<<endl;
	// _IMDforcescale=_IMDforcescale/ForceField::AVOGADRONUMBER;//kcal.A-1
	// //cout<<"_IMDforcescale "<<_IMDforcescale<<" kcal.A-1"<<endl;
	// _IMDforcescale=_IMDforcescale*ForceField::KCAL2KJOULE*1.0E3;//J.A-1
	// //cout<<"_IMDforcescale "<<_IMDforcescale<<" J.A-1"<<endl;
	// _IMDforcescale=_IMDforcescale/ForceField::ANGSTROM2METER;//J.m-1 ou N ou kg.m.s-2
	// //cout<<"_IMDforcescale "<<_IMDforcescale<<" J.m-1 ou N ou kg.m.s-2"<<endl;
	// _IMDforcescale=_IMDforcescale*ForceField::NEWTON2DALTONANGSTROMPERFEMTOSECOND2;//Da.A.fs-2
	// cout<<"_IMDforcescale "<<_IMDforcescale<<" Da.A.fs-2"<<endl;
	// //_IMDlog =  IIMD_init( "", &_IMDmode,&_IMDwait,&_IMDport, &_IMDdebug,0 );
}




void InteractorMDDriver::start(int nbPos)
{
	init(nbPos);
	
	
	Interactor::start(nbPos);
	int rc;
	rc = pthread_create(&_imdthreads, NULL, runthread, (void *) this);
	
	if (rc)
	{
		exit(-1);
	}
}

void InteractorMDDriver::stop()
{
	Interactor::stop();
	_IMDstop=1;
}




void InteractorMDDriver::setPort(unsigned port)
{
	_IMDport = port;
}


void InteractorMDDriver::setWait(unsigned wait)
{
	_IMDwait = wait ;
}

void InteractorMDDriver::setDebug(unsigned debug)
{
	_IMDdebug = debug;
}

void InteractorMDDriver::setLog(const char * logfilename)
{
	strcpy(_IMDlogfilename,logfilename);
}

void InteractorMDDriver::setForceScale(float forcescale)
{
	_IMDforcescale = forcescale;
}

void InteractorMDDriver::init( InteractorMDDriver * imdl)
{
	static int fp_comm = -1;
	if ( fp_comm == -1) 
	{
		// SpringNetwork * spn=imdl->getSpringNetwork();
		if(imdl->_IMDwait!=0)
		{
		 	imdl->_pause = true;
		}
		imdl->_IMDlog =  IIMD_init( "", &(imdl->_IMDmode),&(imdl->_IMDwait),&(imdl->_IMDport), &(imdl->_IMDdebug),0 );
		
		IIMD_probeconnection();
		IIMD_treatprotocol();
		
		if(imdl->_IMDwait!=0)
		{
			imdl->_pause = false;
		}		
		fp_comm = 1;
	}
}

void  * InteractorMDDriver::runthread(void * userdata)
{
	InteractorMDDriver * imdl=(InteractorMDDriver *)userdata;
	pthread_mutex_lock( &imdl->mutex );
	init(imdl);
	pthread_mutex_unlock( &imdl->mutex );
	while ( imdl->_start && imdl->_IMDstop!=1 ) 
	{
		iterate(imdl);
		usleep(1000);
	}
	IIMD_terminate ( );
	exit(0);
	return NULL;	
}




int InteractorMDDriver::iterate( InteractorMDDriver * imdl)
{
	int ret=0;
	pthread_mutex_lock( &imdl->mutex );	
	init(imdl);
	IIMD_probeconnection();
	IIMD_treatprotocol();
	IIMD_send_coords(&(imdl->_nbpositions), imdl->_positions);
	init(imdl);
	IIMD_probeconnection();
	IIMD_treatprotocol();
	IIMD_send_energies( &(imdl->_IMDenergies) );
	init(imdl);
	IIMD_probeconnection();
	IIMD_treatprotocol();
	IIMD_get_forces( &(imdl->_nbforces), &(imdl->_particleforceids), &(imdl->_particleforces) );
	memset(imdl->_forces,0, sizeof(float)*imdl->_nbpositions*3);
	for(unsigned i=0;i<imdl->_nbforces;i++)
	{
		imdl->_forces[(imdl->_particleforceids)[i]*3]=(imdl->_particleforces)[i*3]*imdl->_IMDforcescale;
		imdl->_forces[(imdl->_particleforceids)[i]*3+1]=(imdl->_particleforces)[i*3+1]*imdl->_IMDforcescale;
		imdl->_forces[(imdl->_particleforceids)[i]*3+2]=(imdl->_particleforces)[i*3+2]*imdl->_IMDforcescale;
	}
	switch( imd_event ) 
	{
		case IMD_KILL:
			imdl->_IMDstop = 1;
			// Important : clear the event
			imd_event = -1;
			break;
		case IMD_TRATE:
			//nstximd = imd_value;
			// Important : clear the event
			imd_event = -1;
			break;
	}
	pthread_mutex_unlock( &imdl->mutex );
	// added return to avoid warning for this non-void function
	return ret;
}

void InteractorMDDriver::synchronize(double * pos, double * forces, int nbpos,
								 	 int nstep, double total_energy, 
								 	 double temperature, energies_t energies)
{
	Interactor::synchronize(pos, forces, nbpos);
	_IMDenergies.tstep = nstep;
	_IMDenergies.T = temperature;          //!< Temperature in degrees Kelvin
	// 	_IMDenergies.Epot=0.0;       //!< Potential energy, in Kcal/mol
	_IMDenergies.Evdw = energies.evdw;       //!< Van der Waals energy, in Kcal/mol
		
	_IMDenergies.Eelec = energies.elec;      //!< Electrostatic energy, in Kcal/mol
		
	_IMDenergies.Ebond = energies.ebona + energies.ebonh;      //!< Bond energy, Kcal/mol
	_IMDenergies.Eangle= energies.etha + energies.ethh;     //!< Angle energy, Kcal/mol
  _IMDenergies.Edihe = energies.epa + energies.eph;
	_IMDenergies.Eimpr = 0.0;
	// 	//if(_probeparticule!=NULL)
	// 	//	_IMDenergies.Edihe=_probeparticule->getStericEnergy();      //!< Dihedral energy, Kcal/mol
	// 	//if(_probeparticule!=NULL)
	// 	//	_IMDenergies.Eimpr=_probeparticule->getForce().norm();      //!< Dihedral energy, Kcal/mol
		
		
	// 	_IMDenergies.Eimpr=0.0;      //!< Improper energy, Kcal/mol
	_IMDenergies.Etot = total_energy;       //!< Total energy, in Kcal/mol
}




