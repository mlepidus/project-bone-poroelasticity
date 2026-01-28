// ============================================================================
// BulkDarcyData.h - Darcy flow problem data for bulk domain
// ============================================================================
#ifndef BULKDARCYDATA_H
#define BULKDARCYDATA_H

#include "Core.h"
#include "Parser.h"

/**
 * @class BulkDarcyData
 * @brief Physical parameters and source terms for Darcy flow in bulk domain
 * 
 * Manages permeability tensor, fluid properties, source terms, and
 * exact solutions for verification. Supports spatially variable properties.
 */
class BulkDarcyData {
public:
    /**
     * @brief Construct from data file
     * @param dataFile GetPot parameter file
     * @param section Bulk data section name
     * @param sectionDarcy Darcy-specific parameters section
     */
    BulkDarcyData(const GetPot& dataFile,
                  const std::string& section = "bulkData/",
                  const std::string& sectionDarcy = "darcy/");
    
    // Permeability tensor components (can be spatially variable)
    scalar_type Kxx(const base_node& x);  ///< K_xx component
    scalar_type Kyy(const base_node& x);  ///< K_yy component
    scalar_type Kxy(const base_node& x);  ///< K_xy component (off-diagonal)
    scalar_type Kzz(const base_node& x);  ///< K_zz component
    scalar_type Kxz(const base_node& x);  ///< K_xz component (off-diagonal)
    scalar_type Kyz(const base_node& x);  ///< K_yz component (off-diagonal)
    
    /// Set permeability values at mesh nodes
    void setKxx(std::vector<base_node> nodes);
    void setKxy(std::vector<base_node> nodes);
    void setKyy(std::vector<base_node> nodes);
    void setKzz(std::vector<base_node> nodes);
    void setKxz(std::vector<base_node> nodes);
    void setKyz(std::vector<base_node> nodes);
    
    /// Evaluate spatially distributed permeability
    bgeot::base_node spaceDistrK(const base_node& x);
    
    // Getters for permeability at DOF i
    inline scalar_type getKxx(size_type i) { return (*M_KxxVector)[i]; }
    inline scalar_type getKxy(size_type i) { return (*M_KxyVector)[i]; }
    inline scalar_type getKyy(size_type i) { return (*M_KyyVector)[i]; }
    inline scalar_type getKzz(size_type i) { return (*M_KzzVector)[i]; }
    inline scalar_type getKxz(size_type i) { return (*M_KxzVector)[i]; }
    inline scalar_type getKyz(size_type i) { return (*M_KyzVector)[i]; }
    inline scalar_type getBiotAlpha() { return M_biotAlpha; }
    /// Get fluid density
    inline scalar_type rhoF() { return M_rhoF; }
    
    /// Get gravity vector
    inline bgeot::base_node gravity() { return M_G; }
    
    /// Evaluate source term
    scalar_type source(const base_node& x, const scalar_type t);
    
    /// Evaluate initial pressure condition
    scalar_type pIni(const base_node& x);
    
    /// Evaluate exact pressure solution (for verification)
    scalar_type pEx(const base_node& x, const scalar_type t);
    
    /// Evaluate exact velocity solution (for verification)
    bgeot::base_node uEx(const base_node& x, const scalar_type t);
    
    /// Get Biot modulus M (poroelasticity parameter)
    inline scalar_type M() { return M_biotM; }
    
    /// Load layered material properties
    void getLayers();
    
    bool hasExactSolution() { return !M_pEx.empty() && !M_uEx.empty(); };

    /// Get exact pressure expression string
    inline std::string getPexpr() { return M_pEx; }
    
    /// Get exact velocity expression string
    inline std::string getUexpr() { return M_uEx; }
    
    /// Get leakage coefficient
    inline scalar_type getLeakage() { return M_leakage; }

private:
    std::string M_section;        ///< Bulk data section name
    std::string M_sectionDarcy;   ///< Darcy parameters section name
    
    // Permeability expression strings
    std::string M_Kxx, M_Kyy, M_Kxy;
    std::string M_Kzz, M_Kxz, M_Kyz;
    
    std::string M_source;  ///< Source term expression
    std::string M_pIni;    ///< Initial pressure expression
    std::string M_pEx;     ///< Exact pressure expression
    std::string M_uEx;     ///< Exact velocity expression
    std::string M_Gstring; ///< Gravity vector expression
    
    // Permeability tensor storage (nodal values)
    scalarVectorPtr_Type M_KxxVector;
    scalarVectorPtr_Type M_KyyVector;
    scalarVectorPtr_Type M_KxyVector;
    scalarVectorPtr_Type M_KzzVector;
    scalarVectorPtr_Type M_KxzVector;
    scalarVectorPtr_Type M_KyzVector;
    
    bgeot::base_node M_G;     ///< Gravity acceleration vector
    scalar_type M_rhoF;       ///< Fluid density
    scalar_type M_muF;        ///< Fluid dynamic viscosity
    scalar_type M_biotM;      ///< Biot modulus (poroelasticity)
    scalar_type M_leakage;    ///< Leakage coefficient
    scalar_type M_biotAlpha;  ///< Biot coefficient
    mutable LifeV::Parser M_parser;  ///< Expression parser
};

#endif // BULKDARCYDATA_H