#include "../include/BC.h"

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
    
}

// ============================================================================
// Boundary Assignment from Gmsh Tags
// ============================================================================

void BC::setBoundariesFromTags(getfem::mesh* meshPtr, 
                               const std::map<std::string, size_type>& regmap)
{
    std::cout << "\n=== BC::setBoundariesFromTags ===" << std::endl;
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
    
    std::cout << "\nPhysical names from Gmsh:" << std::endl;
    for (const auto& pair : regmap) {
        std::cout << "  '" << pair.first << "' -> region ID " << pair.second << std::endl;
    }
    
    // Map physical names to internal region numbering
    // Standard naming: outer/ext -> 0, inner/int -> 1, bottom/lower -> 2, top/upper -> 3
    std::map<std::string, size_type> name_to_internal_id;
    
    for (const auto& pair : regmap) {
        std::string name_lower = pair.first;
        std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(), ::tolower);
        
        if (name_lower.find("outer") != std::string::npos || 
            name_lower.find("external") != std::string::npos ||
            name_lower.find("ext") != std::string::npos) {
            name_to_internal_id[pair.first] = 0;
            std::cout << "  Mapping '" << pair.first << "' -> internal region 0 (outer)" << std::endl;
        }
        else if (name_lower.find("inner") != std::string::npos || 
                 name_lower.find("internal") != std::string::npos ||
                 name_lower.find("int") != std::string::npos) {
            name_to_internal_id[pair.first] = 1;
            std::cout << "  Mapping '" << pair.first << "' -> internal region 1 (inner)" << std::endl;
        }
        else if (name_lower.find("bottom") != std::string::npos || 
                 name_lower.find("lower") != std::string::npos ||
                 name_lower.find("low") != std::string::npos ||
                 name_lower.find("bot") != std::string::npos) {
            name_to_internal_id[pair.first] = 2;
            std::cout << "  Mapping '" << pair.first << "' -> internal region 2 (bottom)" << std::endl;
        }
        else if (name_lower.find("top") != std::string::npos || 
                 name_lower.find("upper") != std::string::npos ||
                 name_lower.find("up") != std::string::npos) {
            name_to_internal_id[pair.first] = 3;
            std::cout << "  Mapping '" << pair.first << "' -> internal region 3 (top)" << std::endl;
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
                
                std::cout << "  ✓ Assigned " << face_count << " faces from '" 
                          << physical_name << "' (Gmsh region " << gmsh_region_id 
                          << ") to internal region " << internal_id << std::endl;
            }
            else {
                std::cout << "  ✗ Warning: Gmsh region " << gmsh_region_id 
                          << " for '" << physical_name << "' not found!" << std::endl;
            }
        }
    }
    
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
    
    std::cout << "=== Boundary assignment complete ===\n" << std::endl;
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
// BC Evaluation Functions
// ============================================================================

scalar_type BC::BCNeum(const base_node& x, const size_type& flag, const scalar_type t)
{
    M_parser.setString(M_BCNeum);
    M_parser.setVariable("x", x[0]);
    M_parser.setVariable("y", x[1]);
    M_parser.setVariable("z", x[2]);
    M_parser.setVariable("n", flag);
    M_parser.setVariable("t", t);
    return M_parser.evaluate();
}

scalar_type BC::BCDiri(const base_node& x, const size_type& flag)
{
    M_parser.setString(M_BCDiri);
    M_parser.setVariable("x", x[0]);
    M_parser.setVariable("y", x[1]);
    M_parser.setVariable("z", x[2]);
    M_parser.setVariable("n", flag);
    return M_parser.evaluate();
}

bgeot::base_node BC::BCDiriVec(const base_node& x, const size_type& flag, const scalar_type t)
{
    bgeot::base_node sol(3, 0);
       
    for (size_type i = 0; i < 3; ++i)
    {
        M_parser.setString(M_BCDiriVec);
        M_parser.setVariable("x", x[0]);
        M_parser.setVariable("y", x[1]);
        M_parser.setVariable("z", x[2]);
        M_parser.setVariable("n", flag);
        M_parser.setVariable("t", t);
        sol[i] = M_parser.evaluate(i);
    }
    return sol;
}

bgeot::base_node BC::BCNeumVec(const base_node& x, const size_type& flag, const scalar_type t)
{
    bgeot::base_node sol(3, 0);
       
    for (size_type i = 0; i < 3; ++i)
    {
        M_parser.setString(M_BCNeumVec);
        M_parser.setVariable("x", x[0]);
        M_parser.setVariable("y", x[1]);
        M_parser.setVariable("z", x[2]);
        M_parser.setVariable("n", flag);
        M_parser.setVariable("t", t);
        sol[i] = M_parser.evaluate(i);
    }
    return sol;
}

bgeot::base_node BC::BCDiriVel(const base_node& x, const size_type& flag, const scalar_type t)
{
    bgeot::base_node sol(3, 0);
    
    for (size_type i = 0; i < 3; ++i)
    {
        M_parser.setString(M_BCDiriVel);
        M_parser.setVariable("x", x[0]);
        M_parser.setVariable("y", x[1]);
        M_parser.setVariable("z", x[2]);
        M_parser.setVariable("n", flag);
        M_parser.setVariable("t", t);
        sol[i] = M_parser.evaluate(i);
    }
    return sol;
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