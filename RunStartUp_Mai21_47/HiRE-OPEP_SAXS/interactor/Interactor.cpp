#include <cstdlib>
#include <cstring>
#include <iostream>
#include "Interactor.h"

using namespace std;

Interactor::Interactor()
{
	_nbpositions	= 0;
	_positions 		= NULL;
	_forces			= NULL;
	_start 			= false;
	_pause 			= false;
}

Interactor::~Interactor()
{
	delete [] _positions;
	delete [] _forces;
}

int Interactor::getNbPositions() const
{
	return _nbpositions;
}

void Interactor::setNbPositions(int  nbpositions)
{
	delete[] _positions;
	delete[] _forces;
	_nbpositions=nbpositions;
	_positions=new float[_nbpositions*3];
	memset(_positions, 0, sizeof(float)*_nbpositions*3);	
	_forces=new float[_nbpositions*3];
	memset(_forces, 0, sizeof(float)*_nbpositions*3);

	cout << "Init: " << _nbpositions << endl;	
}
	
void Interactor::setPositions(float * positions)
{
	memcpy(_positions,positions, sizeof(float)*_nbpositions*3);
}

float * Interactor::getForces() const
{
	return _forces;
}	

void Interactor::start(int nbpositions)
{
	_start=true;
}		

void Interactor::stop()
{
	_start=false;
}	

void Interactor::init(int nbpositions)
{
		setNbPositions(nbpositions);
}

int Interactor::isActive()
{
	return !_pause;
}	

void Interactor::synchronize(double * pos, double * forces, int nbpositions)
{
	// if(_springnetwork!=NULL)	
	{
		//cout<<"Interactor::Synchonise"<<endl;
		float position[3];
		for(unsigned j=0;j<nbpositions;j++)
		{
			// cout << j<<"(" << pos[j*3] << ","<< pos[j*3+1]<<","<< pos[j*3+2]<<") - ";
			// cout << endl;
			//memcpy(&(_positions[j*3]), &(pos[j*3]), sizeof(float)*3);
			_positions[j*3] = pos[j*3];
			_positions[j*3+1] = pos[j*3+1];
			_positions[j*3+2] = pos[j*3+2];

			// cout << j<<"(" << _forces[j*3] << ","<< _forces[j*3+1]<<","<< _forces[j*3+2]<<")" << endl;
			forces[j*3]   = _forces[j*3];
			forces[j*3+1] = _forces[j*3+1];
			forces[j*3+2] = _forces[j*3+2];
			// _springnetwork->setForce(j,&(_forces[j*3]));
			
		}
	}
}

