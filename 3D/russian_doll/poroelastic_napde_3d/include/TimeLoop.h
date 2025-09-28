#ifndef TIMELOOP_H
#define TIMELOOP_H

#include "LinearSystem.h"
#include "DarcyOperators.h"
#include "UsefulFunctions.h"
#include "StringUtility.h"

// per casi tempo dipendenti contiene le info relative all'avanzamento in tempo

class TimeLoop
{
public: 
	 TimeLoop ( const GetPot& dataFile); 
	 
	 inline scalar_type dt()
	 {
	 	return M_DT;
	 }
	 inline scalar_type Tend()
	 {
	 	return M_Tend;
	 }
	 inline scalar_type Nstep()
	 {
	 	return M_nstep;
	 }
	 
	 inline scalar_type time()
	 {
	 	return M_currentTime;
	 }
 	inline size_type timeStep()
	 {
	 	return M_currentTimeStep;
	 }
	 inline void advance()
	 {
	 	M_currentTime+=M_DT;
		M_currentTimeStep+=1;
	 }
	
        
 private:
    
    scalar_type M_Tend;
    scalar_type M_DT;
    size_type M_nstep;
    scalar_type M_currentTime;
    size_type M_currentTimeStep;
    

    mutable LifeV::Parser M_parser;
};


#endif
