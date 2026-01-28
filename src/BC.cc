// ============================================================================
// BC.cc - Boundary condition implementation with polynomial callback support
// ============================================================================

#include "../include/BC.h"
#include <algorithm>
#include <cmath>
#include <set>

// ============================================================================
// Constructor
// ============================================================================
BC::BC(const GetPot& dataFile,
       const std::string& problem,
       const std::string& section) :
        M_section(section + problem),
        M_nBoundaries(dataFile((M_section + "nBoundaries").data(), 4)),
        M_BCstring(dataFile((M_section + "bcflag").data(), "[1,1,1,1]")),
        M_BCNeum(dataFile((M_section + "p_BC").data(), "0")),
        M_BCNeumVec(dataFile((M_section + "bdLoad").data(), "[0,0,0]")),        
        M_BCDiri(dataFile((M_section + "v_BC").data(), "0")),
        M_BCDiriVec(dataFile((M_section + "bdDisp").data(), "[0,0,0]")),
        M_BCDiriVel(dataFile((M_section + "bdVel").data(), "[0,0,0]"))
{   
    M_BC.resize(M_nBoundaries, 0); 
    M_parser.setString(M_BCstring);

    std::cout << "[BC] Parsing boundary conditions for " << M_section << std::endl;
    std::cout << "  BC flags: ";
    
    for (size_type i = 0; i < M_nBoundaries; ++i)
    {
        M_BC[i] = M_parser.evaluate(i);
        std::cout << M_BC[i] << " ";

        if (M_BC[i] == 0)
        {
            M_DiriRG.push_back(i);
        }
        else if (M_BC[i] == 1)
        {
            M_NeumRG.push_back(i);
        }
        else if (M_BC[i] == 2)
        {
            M_MixedRG.push_back(i);
        }
    }
    std::cout << std::endl;
    
    std::cout << "  Dirichlet regions: ";
    for (auto r : M_DiriRG) std::cout << r << " ";
    std::cout << std::endl;
    
    std::cout << "  Neumann regions: ";
    for (auto r : M_NeumRG) std::cout << r << " ";
    std::cout << std::endl;
}

// ============================================================================
// Polynomial Callback Methods
// ============================================================================

void BC::setPressureCallback(size_type region, ScalarCallback callback) {
    M_pressureCallbacks[region] = callback;
    std::cout << "[BC] Pressure callback set for region " << region << std::endl;
}

void BC::setDisplacementCallback(size_type region, VectorCallback callback) {
    M_displacementCallbacks[region] = callback;
    std::cout << "[BC] Displacement (vector) callback set for region " << region << std::endl;
}

void BC::setDisplacementCallbacks(size_type region,
                                  ScalarCallback callback_x,
                                  ScalarCallback callback_y,
                                  ScalarCallback callback_z) {
    M_dispXCallbacks[region] = callback_x;
    M_dispYCallbacks[region] = callback_y;
    M_dispZCallbacks[region] = callback_z;
    
    // Also create a combined vector callback
    M_displacementCallbacks[region] = [callback_x, callback_y, callback_z](scalar_type z) -> bgeot::base_node {
        bgeot::base_node result(3);
        result[0] = callback_x ? callback_x(z) : 0.0;
        result[1] = callback_y ? callback_y(z) : 0.0;
        result[2] = callback_z ? callback_z(z) : 0.0;
        return result;
    };
    
    std::cout << "[BC] Displacement (component) callbacks set for region " << region << std::endl;
}

void BC::setPressurePolynomial(size_type region,
                               const std::vector<scalar_type>& coefficients,
                               scalar_type z_min,
                               scalar_type z_max) {
    // Store polynomial data
    M_pressurePolynomials[region] = {coefficients, z_min, z_max};
    
    // Create callback that evaluates the polynomial
    M_pressureCallbacks[region] = [this, region](scalar_type z) -> scalar_type {
        const auto& poly = this->M_pressurePolynomials.at(region);
        return this->evaluatePolynomial(z, poly.coefficients, poly.z_min, poly.z_max);
    };
    
    std::cout << "[BC] Pressure polynomial set for region " << region << std::endl;
    #ifdef VERBOSE
    std::cout << "  Z range: [" << z_min << ", " << z_max << "]" << std::endl;
    std::cout << "  Polynomial order: " << (coefficients.size() - 1) << std::endl;
    std::cout << "  Coefficients: [";
    for (size_t i = 0; i < coefficients.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << coefficients[i];
    }
    std::cout << "]" << std::endl;
    #endif
}

void BC::setDisplacementPolynomial(size_type region,
                                   const std::vector<scalar_type>& coeffs_x,
                                   const std::vector<scalar_type>& coeffs_y,
                                   const std::vector<scalar_type>& coeffs_z,
                                   scalar_type z_min,
                                   scalar_type z_max) {
    // Store polynomial data
    M_dispXPolynomials[region] = {coeffs_x, z_min, z_max};
    M_dispYPolynomials[region] = {coeffs_y, z_min, z_max};
    M_dispZPolynomials[region] = {coeffs_z, z_min, z_max};
    
    // Create component callbacks
    M_dispXCallbacks[region] = [this, region](scalar_type z) -> scalar_type {
        const auto& poly = this->M_dispXPolynomials.at(region);
        return this->evaluatePolynomial(z, poly.coefficients, poly.z_min, poly.z_max);
    };
    
    M_dispYCallbacks[region] = [this, region](scalar_type z) -> scalar_type {
        const auto& poly = this->M_dispYPolynomials.at(region);
        return this->evaluatePolynomial(z, poly.coefficients, poly.z_min, poly.z_max);
    };
    
    M_dispZCallbacks[region] = [this, region](scalar_type z) -> scalar_type {
        const auto& poly = this->M_dispZPolynomials.at(region);
        return this->evaluatePolynomial(z, poly.coefficients, poly.z_min, poly.z_max);
    };
    
    // Create combined vector callback
    M_displacementCallbacks[region] = [this, region](scalar_type z) -> bgeot::base_node {
        bgeot::base_node result(3);
        result[0] = this->M_dispXCallbacks.at(region)(z);
        result[1] = this->M_dispYCallbacks.at(region)(z);
        result[2] = this->M_dispZCallbacks.at(region)(z);
        return result;
    };
    
    std::cout << "[BC] Displacement polynomial set for region " << region << std::endl;
    #ifdef VERBOSE
    std::cout << "  Z range: [" << z_min << ", " << z_max << "]" << std::endl;
    #endif
}

void BC::clearPressureCallback(size_type region) {
    M_pressureCallbacks.erase(region);
    M_pressurePolynomials.erase(region);
    std::cout << "[BC] Pressure callback cleared for region " << region << std::endl;
}

void BC::clearDisplacementCallback(size_type region) {
    M_displacementCallbacks.erase(region);
    M_dispXCallbacks.erase(region);
    M_dispYCallbacks.erase(region);
    M_dispZCallbacks.erase(region);
    M_dispXPolynomials.erase(region);
    M_dispYPolynomials.erase(region);
    M_dispZPolynomials.erase(region);
    std::cout << "[BC] Displacement callback cleared for region " << region << std::endl;
}

bool BC::hasPressureCallback(size_type region) const {
    return M_pressureCallbacks.find(region) != M_pressureCallbacks.end();
}

bool BC::hasDisplacementCallback(size_type region) const {
    return M_displacementCallbacks.find(region) != M_displacementCallbacks.end();
}

scalar_type BC::evaluatePolynomial(scalar_type z,
                                   const std::vector<scalar_type>& coeffs,
                                   scalar_type z_min,
                                   scalar_type z_max) const {
    if (coeffs.empty()) {
        return 0.0;
    }
    
    // Normalize z to [0, 1]
    scalar_type t = 0.0;
    if (std::abs(z_max - z_min) > 1e-15) {
        t = (z - z_min) / (z_max - z_min);
    }
    
    // Clamp to [0, 1]
    t = std::max(0.0, std::min(1.0, t));
    
    // Evaluate polynomial: p(t) = c0 + c1*t + c2*t^2 + ...
    scalar_type result = 0.0;
    scalar_type t_power = 1.0;
    
    for (size_type i = 0; i < coeffs.size(); ++i) {
        result += coeffs[i] * t_power;
        t_power *= t;
    }
    
    return result;
}

// ============================================================================
// BC Evaluation Functions
// ============================================================================

scalar_type BC::BCNeum(const base_node& x, const size_type& flag, const scalar_type t)
{   
    // Check if we have a pressure callback for this region
    auto it = M_pressureCallbacks.find(flag);
    if (it != M_pressureCallbacks.end() && it->second) {
        // Use polynomial callback based on z-coordinate
        return it->second(x[2]);
    }
    
    // Otherwise use standard parser-based evaluation
    size_type dim = x.size();

    M_parser.setString(M_BCNeum);
    M_parser.setVariable("x", x[0]);
    M_parser.setVariable("y", x[1]);
    if (dim >= 3) 
        M_parser.setVariable("z", x[2]);
    else 
        M_parser.setVariable("z", 0.0);
    M_parser.setVariable("n", flag);
    M_parser.setVariable("t", t);
    return M_parser.evaluate();
}

scalar_type BC::BCDiri(const base_node& x, const size_type& flag)
{   
    size_type dim = x.size();

    M_parser.setString(M_BCDiri);
    M_parser.setVariable("x", x[0]);
    M_parser.setVariable("y", x[1]);
    if (dim >= 3) 
        M_parser.setVariable("z", x[2]);
    else 
        M_parser.setVariable("z", 0.0);
    M_parser.setVariable("n", flag);
    return M_parser.evaluate();
}

bgeot::base_node BC::BCDiriVec(const base_node& x, const size_type& flag, const scalar_type t)
{   
    size_type dim = x.size();
    bgeot::base_node sol(dim, 0.0);
    
    // Check if we have a displacement callback for this region
    auto it = M_displacementCallbacks.find(flag);
    if (it != M_displacementCallbacks.end() && it->second) {
        // Use polynomial callback based on z-coordinate
        bgeot::base_node callback_result = it->second(x[2]);
        for (size_type i = 0; i < std::min(dim, (size_type)3); ++i) {
            sol[i] = callback_result[i];
        }
        return sol;
    }
    
    // Otherwise use standard parser-based evaluation
    for (size_type i = 0; i < dim; ++i)
    {
        M_parser.setString(M_BCDiriVec);
        M_parser.setVariable("x", x[0]);
        M_parser.setVariable("y", x[1]);
        if (dim >= 3) 
            M_parser.setVariable("z", x[2]);
        else          
            M_parser.setVariable("z", 0.0);
        M_parser.setVariable("n", flag);
        M_parser.setVariable("t", t);
        sol[i] = M_parser.evaluate(i);
    }
    return sol;
}

bgeot::base_node BC::BCNeumVec(const base_node& x, const size_type& flag, const scalar_type t)
{
    size_type dim = x.size();
    bgeot::base_node sol(dim, 0.0);
       
    for (size_type i = 0; i < dim; ++i)
    {
        M_parser.setString(M_BCNeumVec);
        M_parser.setVariable("x", x[0]);
        M_parser.setVariable("y", x[1]);
        if (dim >= 3) 
            M_parser.setVariable("z", x[2]);
        else          
            M_parser.setVariable("z", 0.0);
        M_parser.setVariable("n", flag);
        M_parser.setVariable("t", t);
        sol[i] = M_parser.evaluate(i);
    }
    return sol;
}

bgeot::base_node BC::BCDiriVel(const base_node& x, const size_type& flag, const scalar_type t)
{
    size_type dim = x.size();
    bgeot::base_node sol(dim, 0.0);
    
    for (size_type i = 0; i < dim; ++i)
    {
        M_parser.setString(M_BCDiriVel);
        M_parser.setVariable("x", x[0]);
        M_parser.setVariable("y", x[1]);
        if (dim >= 3) 
            M_parser.setVariable("z", x[2]);
        else          
            M_parser.setVariable("z", 0.0);
        M_parser.setVariable("n", flag);
        M_parser.setVariable("t", t);
        sol[i] = M_parser.evaluate(i);
    }
    return sol;
}

// ============================================================================
// Boundary Region Setup from Gmsh Tags (by Name)
// ============================================================================

void BC::setBoundariesFromTagsName(getfem::mesh* meshPtr, 
                               const std::map<std::string, size_type>& regmap)
{
    std::cout << "\n=== BC::setBoundariesFromTagsName ===" << std::endl;
    std::cout << "Setting boundaries from Gmsh physical tags..." << std::endl;
    
    // Get all regions that exist in the mesh
    dal::bit_vector mesh_regions = meshPtr->regions_index();
    
    std::cout << "Available regions in mesh:" << std::endl;
    for (dal::bv_visitor i(mesh_regions); !i.finished(); ++i) {
        getfem::mesh_region region = meshPtr->region(i);
        size_t count = 0;
        for (getfem::mr_visitor it(region); !it.finished(); ++it) {
            count++;
        }
        std::cout << "  Region " << i << ": " << count << " faces" << std::endl;
    }

    #ifdef VERBOSE
    std::cout << "\nPhysical names from Gmsh:" << std::endl;
    for (const auto& pair : regmap) {
        std::cout << "  '" << pair.first << "' -> region ID " << pair.second << std::endl;
    }
    #endif

    // Map physical names to internal region numbering
    // Standard naming: outer/ext -> 0, inner/int -> 1, bottom/lower -> 2, top/upper -> 3
    std::map<std::string, size_type> name_to_internal_id;
    
    for (const auto& pair : regmap) {
        std::string name_lower = pair.first;
        std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(), ::tolower);
        
        // Mapping for OUTER or LEFT -> ID 0
        if (name_lower.find("outer") != std::string::npos || 
            name_lower.find("external") != std::string::npos ||
            name_lower.find("ext") != std::string::npos ||
            name_lower.find("left") != std::string::npos) {
            
            name_to_internal_id[pair.first] = 0;
            #ifdef VERBOSE
                 std::cout << "  Mapping '" << pair.first << "' -> internal region 0 (Outer/Left)" << std::endl;
            #endif

        }
        // Mapping for INNER or RIGHT -> ID 1
        else if (name_lower.find("inner") != std::string::npos || 
                 name_lower.find("internal") != std::string::npos ||
                 name_lower.find("int") != std::string::npos ||
                 name_lower.find("right") != std::string::npos) {
            
            name_to_internal_id[pair.first] = 1;
            #ifdef VERBOSE
                std::cout << "  Mapping '" << pair.first << "' -> internal region 1 (Inner/Right)" << std::endl;
            #endif
        }
        // Mapping for BOTTOM -> ID 2
        else if (name_lower.find("bottom") != std::string::npos || 
                 name_lower.find("lower") != std::string::npos ||
                 name_lower.find("low") != std::string::npos ||
                 name_lower.find("bot") != std::string::npos) {
            
            name_to_internal_id[pair.first] = 2;
            #ifdef VERBOSE
                std::cout << "  Mapping '" << pair.first << "' -> internal region 2 (Bottom)" << std::endl;
            #endif
        }
        // Mapping for TOP -> ID 3
        else if (name_lower.find("top") != std::string::npos || 
                 name_lower.find("upper") != std::string::npos ||
                 name_lower.find("up") != std::string::npos) {
            
            name_to_internal_id[pair.first] = 3;
                #ifdef VERBOSE
                    std::cout << "  Mapping '" << pair.first << "' -> internal region 3 (Top)" << std::endl;
                #endif
        }
        else {
            std::cout << "  Warning: No mapping found for '" << pair.first << "'" << std::endl;
        }
    }
    
    // Assign Gmsh regions to internal regions
    std::cout << "\nAssigning Gmsh regions to internal regions..." << std::endl;
    
    for (const auto& mapping : name_to_internal_id) {
        const std::string& physical_name = mapping.first;
        size_type internal_id = mapping.second;
        
        auto it = regmap.find(physical_name);
        if (it != regmap.end()) {
            size_type gmsh_region_id = it->second;
            
            if (mesh_regions.is_in(gmsh_region_id)) {
                getfem::mesh_region gmsh_region = meshPtr->region(gmsh_region_id);
                
                size_t face_count = 0;
                for (getfem::mr_visitor face_it(gmsh_region); !face_it.finished(); ++face_it) {
                    meshPtr->region(internal_id).add(face_it.cv(), face_it.f());
                    face_count++;
                }
                #ifdef VERBOSE
                std::cout << "  ✓ Assigned " << face_count << " faces from '" 
                          << physical_name << "' (Gmsh region " << gmsh_region_id 
                          << ") to internal region " << internal_id << std::endl;
                #endif
            }
    #ifdef VERBOSE
            else {
                std::cout << "  ✗ Warning: Gmsh region " << gmsh_region_id 
                          << " for '" << physical_name << "' not found!" << std::endl;
            }
    #endif
        
        }
    }
    

    #ifdef VERBOSE
    // Verify assignment
    std::cout << "\nVerification - Internal regions after assignment:" << std::endl;
    for (size_type i = 0; i < M_nBoundaries; ++i) {
        size_t count = 0;
        getfem::mesh_region region = meshPtr->region(i);
        for (getfem::mr_visitor it(region); !it.finished(); ++it) {
            count++;
        }
        
        std::string bc_type = "Unknown";
        if (M_BC[i] == 0) bc_type = "Dirichlet";
        else if (M_BC[i] == 1) bc_type = "Neumann";
        else if (M_BC[i] == 2) bc_type = "Mixed";
        
        std::cout << "  Region " << i << ": " << count 
                  << " faces (BC type: " << bc_type << ")" << std::endl;
        
        if (count == 0) {
            std::cout << "    WARNING: No faces assigned!" << std::endl;
        }
    }
    #endif

    std::cout << "=== Boundary assignment complete ===\n" << std::endl;
}

// ============================================================================
// Helper Methods for Tag Number Selection
// ============================================================================

std::vector<size_type> BC::extractSortedTagNumbers(
    const std::map<std::string, size_type>& regmap) const 
{
    // Use a set to get unique sorted values
    std::set<size_type> unique_tags;
    for (const auto& pair : regmap) {
        unique_tags.insert(pair.second);
    }
    
    // Convert to sorted vector
    return std::vector<size_type>(unique_tags.begin(), unique_tags.end());
}

std::vector<size_type> BC::selectTagNumbers(
    const std::vector<size_type>& sortedTags,
    size_type n,
    bool largest) const 
{
    if (sortedTags.size() <= n) {
        // Return all tags if we don't have enough
        return sortedTags;
    }
    
    std::vector<size_type> selected;
    
    if (largest) {
        // Select the N largest (from the end of sorted list)
        size_type start_idx = sortedTags.size() - n;
        for (size_type i = start_idx; i < sortedTags.size(); ++i) {
            selected.push_back(sortedTags[i]);
        }
    } else {
        // Select the N smallest (from the beginning of sorted list)
        for (size_type i = 0; i < n; ++i) {
            selected.push_back(sortedTags[i]);
        }
    }
    
    return selected;  // Already sorted in ascending order
}

// ============================================================================
// Boundary Region Setup from Tag Numbers (Regmap Version)
// ============================================================================

void BC::setBoundariesFromTagNumbers(getfem::mesh* meshPtr, 
                                     const std::map<std::string, size_type>& regmap,
                                     bool largest)
{
    std::cout << "\n=== BC::setBoundariesFromTagNumbers ===" << std::endl;
    std::cout << "Selection mode: " << (largest ? "LARGEST" : "SMALLEST") 
              << " tag numbers" << std::endl;
    std::cout << "Number of BCs required: " << M_nBoundaries << std::endl;
    
    // Clear previous mappings
    M_internalToGmsh.clear();
    M_gmshToInternal.clear();
    
    // Get all regions that exist in the mesh
    dal::bit_vector mesh_regions = meshPtr->regions_index();
    
    std::cout << "\nAvailable regions in mesh:" << std::endl;
    for (dal::bv_visitor i(mesh_regions); !i.finished(); ++i) {
        getfem::mesh_region region = meshPtr->region(i);
        size_t count = 0;
        for (getfem::mr_visitor it(region); !it.finished(); ++it) {
            count++;
        }
        std::cout << "  Region " << i << ": " << count << " faces" << std::endl;
    }
    
    // Extract and sort all unique tag numbers from regmap
    std::vector<size_type> all_tags = extractSortedTagNumbers(regmap);
    
    std::cout << "\nAll tags from regmap (sorted): ";
    for (auto t : all_tags) std::cout << t << " ";
    std::cout << std::endl;
    
    // Filter to only tags that exist in the mesh
    std::vector<size_type> existing_tags;
    for (auto tag : all_tags) {
        if (mesh_regions.is_in(tag)) {
            existing_tags.push_back(tag);
        }
    }
    
    std::cout << "Tags existing in mesh: ";
    for (auto t : existing_tags) std::cout << t << " ";
    std::cout << std::endl;
    
    // Select the N tags based on mode
    std::vector<size_type> selected_tags = selectTagNumbers(existing_tags, M_nBoundaries, largest);
    
    if (selected_tags.size() < M_nBoundaries) {
        std::cout << "WARNING: Only " << selected_tags.size() 
                  << " tags available, but " << M_nBoundaries << " BCs required!" << std::endl;
    }
    
    std::cout << "\nSelected tags: ";
    for (auto t : selected_tags) std::cout << t << " ";
    std::cout << std::endl;
    
    // Assign selected tags to internal regions 0, 1, 2, ...
    std::cout << "\nAssigning tags to internal regions:" << std::endl;
    
    for (size_type internal_id = 0; internal_id < selected_tags.size(); ++internal_id) {
        size_type gmsh_tag = selected_tags[internal_id];
        
        // Store the mapping
        M_internalToGmsh[internal_id] = gmsh_tag;
        M_gmshToInternal[gmsh_tag] = internal_id;
        
        // Copy faces from Gmsh region to internal region
        getfem::mesh_region gmsh_region = meshPtr->region(gmsh_tag);
        
        size_t face_count = 0;
        for (getfem::mr_visitor face_it(gmsh_region); !face_it.finished(); ++face_it) {
            meshPtr->region(internal_id).add(face_it.cv(), face_it.f());
            face_count++;
        }
        
        std::string bc_type = "Unknown";
        if (internal_id < M_BC.size()) {
            if (M_BC[internal_id] == 0) bc_type = "Dirichlet";
            else if (M_BC[internal_id] == 1) bc_type = "Neumann";
            else if (M_BC[internal_id] == 2) bc_type = "Mixed";
        }
        
        std::cout << "  Internal region " << internal_id << " <- Gmsh tag " << gmsh_tag
                  << " (" << face_count << " faces, BC type: " << bc_type << ")" << std::endl;
    }
    
    std::cout << "\n=== Tag number boundary assignment complete ===\n" << std::endl;
}

// ============================================================================
// Boundary Region Setup from Tag Numbers (Direct Mesh Version)
// ============================================================================

void BC::setBoundariesFromTagNumbersDirect(getfem::mesh* meshPtr,
                                           bool largest)
{
    std::cout << "\n=== BC::setBoundariesFromTagNumbersDirect ===" << std::endl;
    std::cout << "Selection mode: " << (largest ? "LARGEST" : "SMALLEST") 
              << " tag numbers" << std::endl;
    std::cout << "Number of BCs required: " << M_nBoundaries << std::endl;
    
    // Clear previous mappings
    M_internalToGmsh.clear();
    M_gmshToInternal.clear();
    
    // Get all non-empty regions from the mesh
    dal::bit_vector mesh_regions = meshPtr->regions_index();
    
    // Collect all region IDs that have faces
    std::vector<size_type> all_tags;
    
    std::cout << "\nScanning mesh regions:" << std::endl;
    for (dal::bv_visitor i(mesh_regions); !i.finished(); ++i) {
        getfem::mesh_region region = meshPtr->region(i);
        size_t count = 0;
        for (getfem::mr_visitor it(region); !it.finished(); ++it) {
            count++;
        }
        if (count > 0) {
            all_tags.push_back(i);
            std::cout << "  Region " << i << ": " << count << " faces" << std::endl;
        }
    }
    
    // Sort the tags
    std::sort(all_tags.begin(), all_tags.end());
    
    std::cout << "\nAll non-empty region tags (sorted): ";
    for (auto t : all_tags) std::cout << t << " ";
    std::cout << std::endl;
    
    // Select the N tags based on mode
    std::vector<size_type> selected_tags = selectTagNumbers(all_tags, M_nBoundaries, largest);
    
    if (selected_tags.size() < M_nBoundaries) {
        std::cout << "WARNING: Only " << selected_tags.size() 
                  << " tags available, but " << M_nBoundaries << " BCs required!" << std::endl;
    }
    
    std::cout << "\nSelected tags: ";
    for (auto t : selected_tags) std::cout << t << " ";
    std::cout << std::endl;
    
    // Store face data temporarily to avoid conflicts during reassignment
    struct FaceData {
        size_type cv;   // Convex index
        short int f;   // Face index
        FaceData(size_type c, short int face) : cv(c), f(face) {}
    };
    std::vector<std::vector<FaceData>> face_storage(selected_tags.size());
    
    // First pass: collect all faces
    for (size_type i = 0; i < selected_tags.size(); ++i) {
        size_type gmsh_tag = selected_tags[i];
        getfem::mesh_region gmsh_region = meshPtr->region(gmsh_tag);
        
        for (getfem::mr_visitor face_it(gmsh_region); !face_it.finished(); ++face_it) {
            face_storage[i].push_back(FaceData(face_it.cv(), face_it.f()));
        }
    }
    
    // Clear internal regions that might overlap with selected tags
    for (size_type internal_id = 0; internal_id < M_nBoundaries; ++internal_id) {
        meshPtr->region(internal_id).clear();
    }
    
    // Second pass: assign faces to internal regions
    std::cout << "\nAssigning tags to internal regions:" << std::endl;
    
    for (size_type internal_id = 0; internal_id < selected_tags.size(); ++internal_id) {
        size_type gmsh_tag = selected_tags[internal_id];
        
        // Store the mapping
        M_internalToGmsh[internal_id] = gmsh_tag;
        M_gmshToInternal[gmsh_tag] = internal_id;
        
        // Add faces to internal region
        for (const auto& face : face_storage[internal_id]) {
            meshPtr->region(internal_id).add(face.cv, face.f);
        }
        
        std::string bc_type = "Unknown";
        if (internal_id < M_BC.size()) {
            if (M_BC[internal_id] == 0) bc_type = "Dirichlet";
            else if (M_BC[internal_id] == 1) bc_type = "Neumann";
            else if (M_BC[internal_id] == 2) bc_type = "Mixed";
        }
        
        std::cout << "  Internal region " << internal_id << " <- Gmsh tag " << gmsh_tag
                  << " (" << face_storage[internal_id].size() << " faces, BC type: " 
                  << bc_type << ")" << std::endl;
    }
    
    std::cout << "\n=== Direct tag number boundary assignment complete ===\n" << std::endl;
}

// ============================================================================
// Geometric Boundary Detection (Fallback)
// ============================================================================

void BC::setBoundaries(getfem::mesh* meshPtr)
{
    std::cout << "\n=== BC::setBoundaries (geometric detection) ===" << std::endl;
    
    getfem::mesh_region border_faces;
    getfem::outer_faces_of_mesh(*meshPtr, border_faces);
    
    int face_count[4] = {0};
    
    for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) 
    {
        assert(i.is_face());
    
        base_node un = meshPtr->normal_of_face_of_convex(i.cv(), i.f());
        un /= gmm::vect_norm2(un);
        
        // Cylindrical case: distinguish inner/outer surfaces based on radial direction
        if (gmm::abs(un[2]) < 0.1)
        {
            base_node bary = gmm::mean_value(meshPtr->points_of_face_of_convex(i.cv(), i.f()));
            
            // Assuming cylinder centered at (1, 1)
            double rx = bary[0] - 1.0;
            double ry = bary[1] - 1.0;
            double s = un[0]*rx + un[1]*ry;
            
            if (s > 0) {
                meshPtr->region(0).add(i.cv(), i.f()); // Outer surface
                face_count[0]++;
            } else {
                meshPtr->region(1).add(i.cv(), i.f()); // Inner surface
                face_count[1]++;
            }
        }
        
        // Bottom surface (z ~ 0)
        if (gmm::abs(un[0]) < 1.0e-7 && gmm::abs(un[1]) < 1.0e-7 && gmm::abs(un[2] + 1) < 1.0e-7)
        {
            meshPtr->region(2).add(i.cv(), i.f());
            face_count[2]++;
        }
        
        // Top surface (z ~ 1)
        if (gmm::abs(un[0]) < 1.0e-7 && gmm::abs(un[1]) < 1.0e-7 && gmm::abs(un[2] - 1) < 1.0e-7)
        {
            meshPtr->region(3).add(i.cv(), i.f());
            face_count[3]++;
        }
    }
    
    std::cout << "Boundary detection complete:" << std::endl;
    for (size_type i = 0; i < M_nBoundaries; ++i) {
        std::cout << "  Region " << i << ": " << face_count[i] << " faces" << std::endl;
    }
    std::cout << "===========================================\n" << std::endl;
}

// ============================================================================
// Getters
// ============================================================================

std::vector<size_type> BC::getNeumBD()
{
    return M_NeumRG;
}

std::vector<size_type> BC::getDiriBD()
{
    return M_DiriRG;
}

std::vector<size_type> BC::getMixedBD()
{
    return M_MixedRG;
}