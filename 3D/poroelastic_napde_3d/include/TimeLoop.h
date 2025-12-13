// ============================================================================
// TimeLoop.h - Time stepping control for transient problems
// ============================================================================
#ifndef TIMELOOP_H
#define TIMELOOP_H

#include "Core.h"
#include "Parser.h"

/**
 * @class TimeLoop
 * @brief Manages time discretization for time-dependent problems
 * 
 * Handles time step size, current time tracking, and advancing through
 * the temporal domain for transient simulations.
 */
class TimeLoop {
public:
    /**
     * @brief Construct from data file
     * @param dataFile GetPot parameter file containing time stepping parameters
     */
    TimeLoop(const GetPot& dataFile);
    
    /// Get time step size
    inline scalar_type dt() { return M_DT; }
    
    /// Get final simulation time
    inline scalar_type Tend() { return M_Tend; }
    
    /// Get total number of time steps
    inline scalar_type Nstep() { return M_nstep; }
    
    /// Get current simulation time
    inline scalar_type time() { return M_currentTime; }
    
    /// Get current time step index
    inline size_type timeStep() { return M_currentTimeStep; }
    
    /// Advance to next time step
    inline void advance() {
        M_currentTime += M_DT;
        M_currentTimeStep += 1;
    }

private:
    scalar_type M_Tend;              ///< Final time
    scalar_type M_DT;                ///< Time step size
    size_type M_nstep;               ///< Number of time steps
    scalar_type M_currentTime;       ///< Current time value
    size_type M_currentTimeStep;     ///< Current time step index
    mutable LifeV::Parser M_parser;  ///< Expression parser for time-dependent functions
};

#endif // TIMELOOP_H