// ============================================================================
// BulkElastData.h - Elasticity problem data for bulk domain
// ============================================================================
#ifndef BULKELASTDATA_H
#define BULKELASTDATA_H

#include "Core.h"
#include "Parser.h"

/**
 * @class BulkElastData
 * @brief Material properties and loading conditions for elasticity in bulk
 * 
 * Manages Lamé parameters, body forces, and exact solutions for elastic
 * deformation problems. Supports spatially variable material properties.
 */
class BulkElastData {
public:
    /**
     * @brief Construct from data file
     * @param dataFile GetPot parameter file
     * @param section Bulk data section name
     * @param sectionElast Elasticity-specific parameters section
     */
    BulkElastData(const GetPot& dataFile,
                  const std::string& section = "bulkData/",
                  const std::string& sectionElast = "mecc/");
    
    /// Evaluate first Lamé parameter (bulk modulus component)
    scalar_type Lambda(const base_node& x);
    
    /// Evaluate second Lamé parameter (shear modulus)
    scalar_type Mu(const base_node& x);
    
    /// Set Lamé parameters at mesh nodes
    void setLambda(std::vector<base_node> nodes);
    void setMu(std::vector<base_node> nodes);
    void setfluidP(std::vector<base_node> nodes);
    
    /// Get Lambda at DOF i
    inline scalar_type getLambda(size_type i) { return (*M_LambdaVector)[i]; }
    
    /// Get Mu at DOF i
    inline scalar_type getMu(size_type i) { return (*M_MuVector)[i]; }
    
    /// Get all Lambda values
    inline std::vector<scalar_type> getLambda() { return *M_LambdaVector; }
    
    /// Get all Mu values
    inline std::vector<scalar_type> getMu() { return *M_MuVector; }
    
    /// Evaluate volumetric body force (e.g., gravity)
    bgeot::base_node bulkLoad(bgeot::base_node x, scalar_type t);
    
    /// Evaluate exact displacement solution (for verification)
    bgeot::base_node uEx(bgeot::base_node x, scalar_type t);
    
    /// Evaluate initial displacement condition
    bgeot::base_node uIni(bgeot::base_node x);
    
    /// Evaluate prescribed fluid pressure field
    scalar_type fluidP(const base_node& x);

private:
    std::string M_section;       ///< Bulk data section name
    std::string M_sectionElast;  ///< Elasticity parameters section name
    
    // Material and loading expressions
    std::string M_mu;       ///< Shear modulus expression
    std::string M_lambda;   ///< First Lamé parameter expression
    std::string M_load;     ///< Body force expression
    std::string M_uEx;      ///< Exact displacement expression
    std::string M_uIni;     ///< Initial displacement expression
    std::string M_fluidP;   ///< Fluid pressure expression
    std::string M_Gstring;  ///< Gravity vector expression
    
    scalar_type M_rhoR;          ///< Solid density
    bgeot::base_node M_G;        ///< Gravity acceleration vector
    
    // Material property storage (nodal values)
    scalarVectorPtr_Type M_LambdaVector;
    scalarVectorPtr_Type M_MuVector;
    scalarVectorPtr_Type M_LoadX;        ///< X-component of body force
    scalarVectorPtr_Type M_LoadY;        ///< Y-component of body force
    scalarVectorPtr_Type M_fluidPVector; ///< Fluid pressure field
    
    LifeV::Parser M_parser;  ///< Expression parser
};

#endif // BULKELASTDATA_H