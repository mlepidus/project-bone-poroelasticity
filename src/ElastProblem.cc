#include "../include/ElastProblem.h"

// ============================================================================
// Static helper to parse boundary type from string
// ============================================================================
BoundaryAssignmentType ElastProblem::parseBoundaryType(const std::string& typeStr) {
    std::string lower = typeStr;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

    if (lower == "geometric_cylinder" || lower == "geo_cylinder" || lower == "auto_cylinder" || lower == "cylinder") {
        return BoundaryAssignmentType::GEOMETRIC_CYLINDER;
    } else if (lower == "geometric_square" || lower == "geo_square" || lower == "auto_square" || lower == "square") {
        return BoundaryAssignmentType::GEOMETRIC_SQUARE;
    } else if (lower == "tagname" || lower == "tag_name" || lower == "name" || lower == "byname") {
        return BoundaryAssignmentType::TAG_NAME;
    } else if (lower == "tagnumber" || lower == "tag_number" || lower == "number" || lower == "bynumber") {
        return BoundaryAssignmentType::TAG_NUMBER;
    }
    
    // Default fallback
    std::cout << "[ElastProblem] WARNING: Unknown boundary type '" << typeStr 
              << "', defaulting to GEOMETRIC_SQUARE" << std::endl;
    return BoundaryAssignmentType::GEOMETRIC_SQUARE;
}

// ============================================================================
// Constructor
// ============================================================================
ElastProblem::ElastProblem(const GetPot& dataFile, Bulk* bulk, 
                          const std::string& basePath):
    M_time(nullptr),
    M_Bulk(bulk),
    M_BC(dataFile, "mecc/", basePath),
    M_boundaryType(BoundaryAssignmentType::GEOMETRIC_SQUARE),
    M_tagNumberLargest(true),
    M_DispFEM(bulk->getMesh(), dataFile, "mecc/", "Displacement", basePath),
    M_CoeffFEM(bulk->getMesh(), dataFile, "mecc/", "Coeff", basePath),
    M_Sys(nullptr),
    M_nbTotDOF(0),
    M_intMethod(*(bulk->getMesh()))
{
    M_nbTotDOF = M_DispFEM.nb_dof();
    
    // Read boundary assignment configuration
    std::string boundaryTypeStr = dataFile((basePath + "domain/boundaryType").c_str(), "geometric");
    M_boundaryType = parseBoundaryType(boundaryTypeStr);
    
    // For TAG_NUMBER mode: read whether to select largest or smallest
    std::string tagSelectionStr = dataFile((basePath + "domain/tagSelection").c_str(), "largest");
    std::transform(tagSelectionStr.begin(), tagSelectionStr.end(), tagSelectionStr.begin(), ::tolower);
    M_tagNumberLargest = (tagSelectionStr != "smallest");
    
    std::cout << "[ElastProblem] Boundary assignment type: ";
    switch (M_boundaryType) {
        case BoundaryAssignmentType::GEOMETRIC_CYLINDER:  std::cout << "GEOMETRIC_CYLINDER"; break;
        case BoundaryAssignmentType::GEOMETRIC_SQUARE:    std::cout << "GEOMETRIC_SQUARE"; break;
        case BoundaryAssignmentType::TAG_NAME:   std::cout << "TAG_NAME"; break;
        case BoundaryAssignmentType::TAG_NUMBER: 
            std::cout << "TAG_NUMBER (" << (M_tagNumberLargest ? "largest" : "smallest") << ")"; 
            break;
    }
    std::cout << std::endl;
    
    // SMART DEFAULT: Check mesh dimension
    size_type meshDim = bulk->getMesh()->dim();
    std::string defaultMethod = (meshDim == 2) ? "IM_TRIANGLE(2)" : "IM_TETRAHEDRON(2)";

    std::string intMethod(dataFile((basePath + "mecc/integrationMethod").c_str(), defaultMethod.c_str()));
    
    M_intMethod.set_integration_method(bulk->getMesh()->convex_index(), 
                                       getfem::int_method_descriptor(intMethod));
    
    // Set material properties
    M_Bulk->getElastData()->setLambda(M_CoeffFEM.getDOFpoints());
    M_Bulk->getElastData()->setMu(M_CoeffFEM.getDOFpoints());
    M_Bulk->getElastData()->setfluidP(M_CoeffFEM.getDOFpoints());
    
    // Setup boundary conditions based on type
    setupBoundaryConditions();
}

// ============================================================================
// Setup Boundary Conditions
// ============================================================================
void ElastProblem::setupBoundaryConditions() {
    switch (M_boundaryType) {
        case BoundaryAssignmentType::GEOMETRIC_CYLINDER:
            std::cout << "[ElastProblem] Using geometric detection for boundary conditions..." << std::endl;
            M_BC.setBoundariesCylinder(M_Bulk->getMesh());
            break;
            
        case BoundaryAssignmentType::GEOMETRIC_SQUARE:
            std::cout << "[ElastProblem] Using geometric detection for boundary conditions..." << std::endl;
            M_BC.setBoundariesSquare(M_Bulk->getMesh());
            break;

        case BoundaryAssignmentType::TAG_NAME:
            if (M_Bulk->hasExternalMesh()) {
                std::cout << "[ElastProblem] Using Gmsh physical tag names for boundary conditions..." << std::endl;
                M_BC.setBoundariesFromTagsName(M_Bulk->getMesh(), M_Bulk->getRegionMap());
            } else {
                std::cout << "[ElastProblem] WARNING: TAG_NAME requested but no external mesh. "
                          << "Falling back to geometric detection." << std::endl;
                M_BC.setBoundariesSquare(M_Bulk->getMesh());
            }
            break;
            
        case BoundaryAssignmentType::TAG_NUMBER:
            if (M_Bulk->hasExternalMesh()) {
                std::cout << "[ElastProblem] Using Gmsh tag numbers for boundary conditions..." << std::endl;
                M_BC.setBoundariesFromTagNumbersDirect(M_Bulk->getMesh(), M_tagNumberLargest);
            } else {
                std::cout << "[ElastProblem] WARNING: TAG_NUMBER requested but no external mesh. "
                          << "Falling back to geometric detection." << std::endl;
                M_BC.setBoundariesSquare(M_Bulk->getMesh());
            }
            break;
    }
}

void ElastProblem::addToSys(LinearSystem* sys) {
    M_Sys = sys;
    M_Sys->addToMatrix(M_nbTotDOF);
}

size_type ElastProblem::getNDOF(std::string variable) {
    if (variable == "Disp" || variable == "all") {
        return M_DispFEM.nb_dof();
    }
    std::cerr << "ERROR: Unknown variable '" << variable << "' in ElastProblem::getNDOF" << std::endl;
    return 0;
}

void ElastProblem::initialize() {
    M_DispSol.reset(new scalarVector_Type(M_DispFEM.nb_dof()));
    M_DispSolOld.reset(new scalarVector_Type(M_DispFEM.nb_dof()));
    gmm::clear(*M_DispSol);
    gmm::clear(*M_DispSolOld);
    
    size_type meshDim = M_Bulk->getMesh()->dim();
    int counter = 0;
    for (size_type i = 0; i < M_DispFEM.nb_dof("base"); i += meshDim) {
        base_node nodo(M_DispFEM.point_of_basic_dof(i));
        
        (*M_DispSol)[counter]     = M_Bulk->getElastData()->uIni(nodo)[0];
        (*M_DispSol)[counter + 1] = M_Bulk->getElastData()->uIni(nodo)[1];
        
        if (meshDim == 3) {
            (*M_DispSol)[counter + 2] = M_Bulk->getElastData()->uIni(nodo)[2];
        }
        
        counter += meshDim;
    }
}

void ElastProblem::assembleMatrix() {
    sparseMatrixPtr_Type A;
    A.reset(new sparseMatrix_Type(M_DispFEM.nb_dof(), M_DispFEM.nb_dof()));
    gmm::clear(*A);
    
    if (!M_dispMassMatrix) {
        M_dispMassMatrix.reset(new sparseMatrix_Type(M_DispFEM.nb_dof(), M_DispFEM.nb_dof()));
        gmm::clear(*M_dispMassMatrix);
        massMatrix(M_dispMassMatrix, M_DispFEM, M_intMethod);
    }
    
    stiffElast(A, M_Bulk, M_DispFEM, M_CoeffFEM, M_intMethod);
    
    M_Sys->addSubMatrix(A, 0, 0);
}

void ElastProblem::assembleRHS() {
    scalarVectorPtr_Type source;
    source.reset(new scalarVector_Type(M_DispFEM.nb_dof()));
    gmm::clear(*source);
    bulkLoad(source, M_Bulk, M_DispFEM, M_CoeffFEM, M_intMethod, M_time->time());
    M_Sys->addSubVector(source, 0);
    
    scalarVectorPtr_Type BCvec;
    BCvec.reset(new scalarVector_Type(M_DispFEM.nb_dof()));
    gmm::clear(*BCvec);
    stressRHS(BCvec, M_Bulk, M_time->time(), &M_BC, M_DispFEM, M_CoeffFEM, M_intMethod);
    M_Sys->addSubVector(BCvec, 0);
}

void ElastProblem::enforceStrongBC(bool firstTime) {
    size_type meshDim = M_Bulk->getMesh()->dim();

    if (firstTime) {
        for (size_type bndID = 0; bndID < M_BC.getDiriBD().size(); bndID++) {
            dal::bit_vector quali = M_DispFEM.getFEM()->dof_on_region(M_BC.getDiriBD()[bndID]);
            for (dal::bv_visitor i(quali); !i.finished(); ++i) {
                M_rowsStrongBC.push_back(i);
                M_rowsStrongBCFlags.push_back(bndID);
            }
        }
        
        for (size_type i = 0; i < M_rowsStrongBC.size(); i += meshDim) {
            
            // --- X COMPONENT ---
            size_type ii = M_rowsStrongBC[i];
            bgeot::base_node where = M_DispFEM.getFEM()->point_of_basic_dof(ii);
            
            M_Sys->setNullRow(ii);
            M_Sys->setMatrixValue(ii, ii, 1);
            scalar_type value = (M_BC.BCDiriVec(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[0] +
                                M_BC.BCDiriVel(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[0] * M_time->time());
            M_Sys->setRHSValue(ii, value);
            
            // --- Y COMPONENT ---
            ii = M_rowsStrongBC[i + 1];
            M_Sys->setNullRow(ii);
            M_Sys->setMatrixValue(ii, ii, 1);
            value = (M_BC.BCDiriVec(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[1] +
                    M_BC.BCDiriVel(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[1] * M_time->time());
            M_Sys->setRHSValue(ii, value);
            
            // --- Z COMPONENT (SOLO 3D) ---
            if (meshDim == 3) {
                ii = M_rowsStrongBC[i + 2];
                M_Sys->setNullRow(ii);
                M_Sys->setMatrixValue(ii, ii, 1);
                value = (M_BC.BCDiriVec(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[2] +
                        M_BC.BCDiriVel(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[2] * M_time->time());
                M_Sys->setRHSValue(ii, value);
            }
        }
    } else {
        for (size_type i = 0; i < M_rowsStrongBC.size(); i += meshDim) {
            
            size_type ii = M_rowsStrongBC[i];
            bgeot::base_node where = M_DispFEM.getFEM()->point_of_basic_dof(ii);
            scalar_type value = (M_BC.BCDiriVec(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[0] +
                                M_BC.BCDiriVel(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[0] * M_time->time());
            M_Sys->setRHSValue(ii, value);
            
            ii = M_rowsStrongBC[i + 1];
            value = (M_BC.BCDiriVec(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[1] +
                    M_BC.BCDiriVel(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[1] * M_time->time());
            M_Sys->setRHSValue(ii, value);
            
            if (meshDim == 3) {
                ii = M_rowsStrongBC[i + 2];
                value = (M_BC.BCDiriVec(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[2] +
                        M_BC.BCDiriVel(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[2] * M_time->time());
                M_Sys->setRHSValue(ii, value);
            }
        }
    }
}

void ElastProblem::solve() {
    M_Sys->solve();
    
    M_DispSol.reset(new scalarVector_Type(M_DispFEM.nb_dof()));
    gmm::clear(*M_DispSol);
    M_Sys->extractSubVector(M_DispSol, 0, "sol");
}

void ElastProblem::extractSol(scalarVectorPtr_Type sol) {
    M_DispSol.reset(new scalarVector_Type(M_DispFEM.nb_dof()));
    gmm::clear(*M_DispSol);
    gmm::copy(gmm::sub_vector(*sol, gmm::sub_interval(0, M_DispFEM.nb_dof())), *M_DispSol);
}

void ElastProblem::updateSol() {
    gmm::clear(*M_DispSolOld);
    gmm::copy(*M_DispSol, *M_DispSolOld);
}

scalar_type ElastProblem::computeError(std::string what, scalar_type time) {
    scalar_type error = 0.0;
    
    if (what == "Displacement") {
        scalarVector_Type loc_err(M_DispFEM.nb_dof());
        
        size_type qdim = M_DispFEM.get_qdim();
        for (size_type i = 0; i < M_DispFEM.nb_dof(); ++i) {
            size_type comp = i % qdim;
            size_type base_dof = i / qdim;
            
            bgeot::base_node where = M_DispFEM.getFEM()->point_of_basic_dof(base_dof);
            base_small_vector exact_val = M_Bulk->getElastData()->uEx(where, time);
            
            loc_err[i] = (*M_DispSol)[i] - exact_val[comp];
        }
        
        error = L2Norm_Elast(M_dispMassMatrix, loc_err);
        std::cout << "at time " << time << " error displacement " << error << std::endl;
    } else {
        std::cerr << "WARNING: Only displacement error computation supported" << std::endl;
    }
    
    return error;
}

void ElastProblem::exportVtk(std::string folder, int frame) {
    std::cout << "Exporting elasticity solution..." << std::endl;
    
    std::string frameNum;
    if (frame >= 0) {
        if (frame < 10) {
            frameNum = "00" + LifeV::number2string(frame);
        } else if (frame < 100) {
            frameNum = "0" + LifeV::number2string(frame);
        } else {
            frameNum = LifeV::number2string(frame);
        }
    }
    
    std::string filename = folder + "Displacement_bulk" + frameNum + ".vtk";
    getfem::vtk_export exp(filename);
    
    exp.exporting(*(M_DispFEM.getFEM()));
    std::vector<scalar_type> dispCut(M_DispFEM.getFEM()->nb_dof(), 0.0);
    gmm::copy(*M_DispSol, dispCut);
    exp.write_mesh();
    exp.write_point_data(*(M_DispFEM.getFEM()), dispCut, "u");
    
    if (M_Bulk->getElastData()->hasExactSolution()) {
        std::string exactFilename = folder + "Displacement_exact" + frameNum + ".vtk";
        getfem::vtk_export expE(exactFilename);
        
        expE.exporting(*(M_DispFEM.getFEM()));
        gmm::clear(dispCut);
        for (size_type i = 0; i < dispCut.size(); i += M_Bulk->getMesh()->dim()) {
            bgeot::base_node where = M_DispFEM.getFEM()->point_of_basic_dof(i);
            dispCut[i] = M_Bulk->getElastData()->uEx(where, M_time->time())[0];
            dispCut[i + 1] = M_Bulk->getElastData()->uEx(where, M_time->time())[1];
            
            if (M_Bulk->getMesh()->dim() == 3) {
                 dispCut[i + 2] = M_Bulk->getElastData()->uEx(where, M_time->time())[2];
            }
        }
        
        expE.write_mesh();
        expE.write_point_data(*(M_DispFEM.getFEM()), dispCut, "u");
    }
}