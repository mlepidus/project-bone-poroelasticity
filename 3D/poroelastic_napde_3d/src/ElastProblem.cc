#include "../include/ElastProblem.h"

ElastProblem::ElastProblem(const GetPot& dataFile, Bulk* bulk):
    M_Bulk(bulk),
    M_BC(dataFile, "mecc/"),
    M_DispFEM(bulk->getMesh(), dataFile, "mecc/", "Displacement", "bulkData/"),
    M_CoeffFEM(bulk->getMesh(), dataFile, "mecc/", "Coeff", "bulkData/"),
    M_Sys(NULL),
    M_intMethod(*(bulk->getMesh()))
{
    M_nbTotDOF = M_DispFEM.nb_dof();
    
    std::string intMethod(dataFile(std::string("bulkData/mecc/integrationMethod").data(), "IM_TETRAHEDRON(2)"));
    M_intMethod.set_integration_method(bulk->getMesh()->convex_index(), 
                                       getfem::int_method_descriptor(intMethod));
    
    // Set material properties
    M_Bulk->getElastData()->setLambda(M_CoeffFEM.getDOFpoints());
    M_Bulk->getElastData()->setMu(M_CoeffFEM.getDOFpoints());
    M_Bulk->getElastData()->setfluidP(M_CoeffFEM.getDOFpoints());
    
    // Set up boundary conditions
    if (M_Bulk->hasExternalMesh()) {
        std::cout << "Using Gmsh physical tags for boundary conditions..." << std::endl;
        M_BC.setBoundariesFromTags(M_Bulk->getMesh(), M_Bulk->getRegionMap());
    } else {
        std::cout << "Using geometric detection for boundary conditions..." << std::endl;
        M_BC.setBoundaries(M_Bulk->getMesh());
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
    
    // Initialize with initial condition if provided
    int counter = 0;
    for (size_type i = 0; i < M_DispFEM.nb_dof("base"); i += 3) {
        base_node nodo(M_DispFEM.point_of_basic_dof(i));
        (*M_DispSol)[counter] = M_Bulk->getElastData()->uIni(nodo)[0];
        (*M_DispSol)[counter + 1] = M_Bulk->getElastData()->uIni(nodo)[1];
        (*M_DispSol)[counter + 2] = M_Bulk->getElastData()->uIni(nodo)[2];
        counter += 3;
    }
}

void ElastProblem::assembleMatrix() {
    sparseMatrixPtr_Type A;
    A.reset(new sparseMatrix_Type(M_DispFEM.nb_dof(), M_DispFEM.nb_dof()));
    gmm::clear(*A);
    
    // Build displacement mass matrix for error computation (only once)
    if (!M_dispMassMatrix) {
        M_dispMassMatrix.reset(new sparseMatrix_Type(M_DispFEM.nb_dof(), M_DispFEM.nb_dof()));
        gmm::clear(*M_dispMassMatrix);
        massMatrix(M_dispMassMatrix, M_DispFEM, M_intMethod);
    }
    
    // Assemble stiffness matrix
    stiffElast(A, M_Bulk, M_DispFEM, M_CoeffFEM, M_intMethod);
    
    // Add to global system
    M_Sys->addSubMatrix(A, 0, 0);
}

void ElastProblem::assembleRHS() {
    // Body forces
    scalarVectorPtr_Type source;
    source.reset(new scalarVector_Type(M_DispFEM.nb_dof()));
    gmm::clear(*source);
    bulkLoad(source, M_Bulk, M_DispFEM, M_CoeffFEM, M_intMethod, M_time->time());
    M_Sys->addSubVector(source, 0);
    
    // Neumann boundary conditions (tractions)
    scalarVectorPtr_Type BCvec;
    BCvec.reset(new scalarVector_Type(M_DispFEM.nb_dof()));
    gmm::clear(*BCvec);
    stressRHS(BCvec, M_Bulk, M_time->time(), &M_BC, M_DispFEM, M_CoeffFEM, M_intMethod);
    M_Sys->addSubVector(BCvec, 0);
}

void ElastProblem::enforceStrongBC(bool firstTime) {
    if (firstTime) {
        // Collect DOFs with Dirichlet BC
        for (size_type bndID = 0; bndID < M_BC.getDiriBD().size(); bndID++) {
            dal::bit_vector quali = M_DispFEM.getFEM()->dof_on_region(M_BC.getDiriBD()[bndID]);
            for (dal::bv_visitor i(quali); !i.finished(); ++i) {
                M_rowsStrongBC.push_back(i);
                M_rowsStrongBCFlags.push_back(bndID);
            }
        }
        
        // Enforce Dirichlet BC in matrix and RHS
        for (size_type i = 0; i < M_rowsStrongBC.size(); i += 3) {
            size_type ii = M_rowsStrongBC[i];
            bgeot::base_node where = M_DispFEM.getFEM()->point_of_basic_dof(ii);
            
            // Set row to identity
            M_Sys->setNullRow(ii);
            M_Sys->setMatrixValue(ii, ii, 1);
            scalar_type value = (M_BC.BCDiriVec(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[0] +
                                M_BC.BCDiriVel(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[0] * M_time->time());
            M_Sys->setRHSValue(ii, value);
            
            // y-component
            ii = M_rowsStrongBC[i + 1];
            M_Sys->setNullRow(ii);
            M_Sys->setMatrixValue(ii, ii, 1);
            value = (M_BC.BCDiriVec(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[1] +
                    M_BC.BCDiriVel(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[1] * M_time->time());
            M_Sys->setRHSValue(ii, value);
            
            // z-component
            ii = M_rowsStrongBC[i + 2];
            M_Sys->setNullRow(ii);
            M_Sys->setMatrixValue(ii, ii, 1);
            value = (M_BC.BCDiriVec(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[2] +
                    M_BC.BCDiriVel(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[2] * M_time->time());
            M_Sys->setRHSValue(ii, value);
        }
    } else {
        // Only update RHS values for time-dependent BC
        for (size_type i = 0; i < M_rowsStrongBC.size(); i += 3) {
            size_type ii = M_rowsStrongBC[i];
            bgeot::base_node where = M_DispFEM.getFEM()->point_of_basic_dof(ii);
            
            scalar_type value = (M_BC.BCDiriVec(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[0] +
                                M_BC.BCDiriVel(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[0] * M_time->time());
            M_Sys->setRHSValue(ii, value);
            
            ii = M_rowsStrongBC[i + 1];
            value = (M_BC.BCDiriVec(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[1] +
                    M_BC.BCDiriVel(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[1] * M_time->time());
            M_Sys->setRHSValue(ii, value);
            
            ii = M_rowsStrongBC[i + 2];
            value = (M_BC.BCDiriVec(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[2] +
                    M_BC.BCDiriVel(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[2] * M_time->time());
            M_Sys->setRHSValue(ii, value);
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
        // Create local error vector
        scalarVector_Type loc_err(M_DispFEM.nb_dof());
        
        // Compute pointwise error
        size_type qdim = M_DispFEM.get_qdim(); // Get dimension of vector field
        for (size_type i = 0; i < M_DispFEM.nb_dof(); ++i) {
            size_type comp = i % qdim; // Component (0=x, 1=y, 2=z)
            size_type base_dof = i / qdim; // Base dof index
            
            bgeot::base_node where = M_DispFEM.getFEM()->point_of_basic_dof(base_dof);
            base_small_vector exact_val = M_Bulk->getElastData()->uEx(where, time);
            
            loc_err[i] = (*M_DispSol)[i] - exact_val[comp];
        }
        
        // Compute L2 norm of error using mass matrix
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
    
    // Export displacement
    std::string filename = folder + "Displacement_bulk" + frameNum + ".vtk";
    getfem::vtk_export exp(filename);
    
    exp.exporting(*(M_DispFEM.getFEM()));
    std::vector<scalar_type> dispCut(M_DispFEM.getFEM()->nb_dof(), 0.0);
    gmm::copy(*M_DispSol, dispCut);
    exp.write_mesh();
    exp.write_point_data(*(M_DispFEM.getFEM()), dispCut, "u");
    
    // Export exact solution for comparison
    if (M_Bulk->getElastData()->hasExactSolution()) {
        std::string exactFilename = folder + "Displacement_exact" + frameNum + ".vtk";
        getfem::vtk_export expE(exactFilename);
        
        expE.exporting(*(M_DispFEM.getFEM()));
        gmm::clear(dispCut);
        for (size_type i = 0; i < dispCut.size(); i += 3) {
            bgeot::base_node where = M_DispFEM.getFEM()->point_of_basic_dof(i);
            dispCut[i] = M_Bulk->getElastData()->uEx(where, M_time->time())[0];
            dispCut[i + 1] = M_Bulk->getElastData()->uEx(where, M_time->time())[1];
            dispCut[i + 2] = M_Bulk->getElastData()->uEx(where, M_time->time())[2];
        }
        
        expE.write_mesh();
        expE.write_point_data(*(M_DispFEM.getFEM()), dispCut, "u");
    }
}