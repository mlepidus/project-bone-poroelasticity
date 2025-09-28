#include "../include/TimeLoop.h"

const std::string section="time/";

TimeLoop::TimeLoop (const GetPot& dataFile):
		M_Tend( dataFile ( ( section + "Tend" ).data (), 1.) ),
		M_DT( dataFile ( ( section + "dt" ).data (), 0.01) )

{       M_currentTime=0;
	M_currentTimeStep=0;
	M_nstep= nearbyint(M_Tend/M_DT);
	std::cout <<"N steps:"<<M_nstep<<" end time:"<<M_Tend<<"  delta t:"<<M_DT <<std::endl;
 
}


