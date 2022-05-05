#ifndef _INTERACTOR_H_
#define _INTERACTOR_H_

#define FILTERSIZE 256
#define LOGFILENAMESIZE 256

// #include "SpringNetwork.h"

class Interactor 
{
	public : 
		Interactor();
		virtual ~Interactor();
		
		// void setSpringNetwork(SpringNetwork * springnetwork);
		// SpringNetwork *  getSpringNetwork() const;
		
		virtual int getNbPositions() const;
		virtual void setPositions(float * positions);
		virtual float * getForces() const;
		
		virtual void start(int nbPos);
		virtual void stop();
		virtual void synchronize(double * pos, double * forces, int nbpositions);

		//Return 0 if pause=true; 1 if pause=false
		virtual int isActive();
		
		

	protected : 
		// SpringNetwork * _springnetwork;
		int _nbpositions;
		float * _positions;
		float * _forces;
		bool _start;
		bool _pause;
		void setNbPositions(int nbpositions);
		virtual void init(int nbpositions);		
		
			
};

#endif

