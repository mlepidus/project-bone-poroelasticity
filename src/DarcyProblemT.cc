#include "../include/DarcyProblemT.h"



// ============================================================================
// Static helper to parse boundary type from string
// ============================================================================
BoundaryAssignmentType DarcyProblemT::parseBoundaryType(const std::string& typeStr) {
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
    std::cout << "[DarcyProblemT] WARNING: Unknown boundary type '" << typeStr 
              << "', defaulting to GEOMETRIC_SQUARE" << std::endl;
    return BoundaryAssignmentType::GEOMETRIC_SQUARE;
}

// ============================================================================
// Constructor
// ============================================================================
DarcyProblemT::DarcyProblemT (const GetPot& dataFile, Bulk* bulk,
                              const std::string& basePath ):
    M_timeLoop(nullptr),
    M_Bulk(bulk),
    M_BC(dataFile, "darcy/", basePath),
    M_boundaryType(BoundaryAssignmentType::GEOMETRIC_SQUARE),
    M_tagNumberLargest(true),
    M_PressureFEM( bulk->getMesh(), dataFile, "darcy/", "Pressure", basePath),
    M_CoeffFEM( bulk->getMesh(), dataFile, "darcy/", "Coeff", basePath),
    M_VelocityFEM( bulk->getMesh(), dataFile, "darcy/", "Velocity", basePath),
    M_visualizationVFEM( bulk->getMesh(), dataFile, "darcy/", "VelocityVis", basePath),
    M_Sys(nullptr),
    M_nbTotDOF(0),			
    M_intMethod(*(bulk->getMesh()))
{  
    M_nbTotDOF = M_PressureFEM.nb_dof() + M_VelocityFEM.nb_dof();
  
    // Read boundary assignment configuration
    std::string boundaryTypeStr = dataFile((basePath + "domain/boundaryType").c_str(), "geometric");
    M_boundaryType = parseBoundaryType(boundaryTypeStr);
    
    // For TAG_NUMBER mode: read whether to select largest or smallest
    std::string tagSelectionStr = dataFile((basePath + "domain/tagSelection").c_str(), "largest");
    std::transform(tagSelectionStr.begin(), tagSelectionStr.end(), tagSelectionStr.begin(), ::tolower);
    M_tagNumberLargest = (tagSelectionStr != "smallest");
    
    std::cout << "[DarcyProblemT] Boundary assignment type: ";
    switch (M_boundaryType) {
        case BoundaryAssignmentType::GEOMETRIC_CYLINDER:  std::cout << "GEOMETRIC_CYLINDER"; break;
        case BoundaryAssignmentType::GEOMETRIC_SQUARE:    std::cout << "GEOMETRIC_SQUARE"; break;
        case BoundaryAssignmentType::TAG_NAME:   std::cout << "TAG_NAME"; break;
        case BoundaryAssignmentType::TAG_NUMBER: 
            std::cout << "TAG_NUMBER (" << (M_tagNumberLargest ? "largest" : "smallest") << ")"; 
            break;
    }
    std::cout << std::endl;
  
    // SMART DEFAULT for integration method
    size_type meshDim = bulk->getMesh()->dim();
    std::string defaultMethod = (meshDim == 2) ? "IM_TRIANGLE(2)" : "IM_TETRAHEDRON(2)";

    std::string intMethod(dataFile((basePath + "darcy/integrationMethod").c_str(), defaultMethod.c_str()));
   
    M_intMethod.set_integration_method(bulk->getMesh()->convex_index(), getfem::int_method_descriptor(intMethod));
   
    M_Bulk->getDarcyData()->setKxx(M_CoeffFEM.getDOFpoints());
    M_Bulk->getDarcyData()->setKyy(M_CoeffFEM.getDOFpoints());
    M_Bulk->getDarcyData()->setKxy(M_CoeffFEM.getDOFpoints());
    M_Bulk->getDarcyData()->setKzz(M_CoeffFEM.getDOFpoints());
    M_Bulk->getDarcyData()->setKxz(M_CoeffFEM.getDOFpoints());
    M_Bulk->getDarcyData()->setKyz(M_CoeffFEM.getDOFpoints());

    // Setup boundary conditions based on type
    setupBoundaryConditions();
}

// ============================================================================
// Setup Boundary Conditions
// ============================================================================
void DarcyProblemT::setupBoundaryConditions() {
    switch (M_boundaryType) {
        case BoundaryAssignmentType::GEOMETRIC_CYLINDER:
            std::cout << "[DarcyProblemT] Using geometric detection for boundary conditions..." << std::endl;
            M_BC.setBoundariesCylinder(M_Bulk->getMesh());
            break;
        case BoundaryAssignmentType::GEOMETRIC_SQUARE:
            std::cout << "[DarcyProblemT] Using geometric detection for boundary conditions..." << std::endl;
            M_BC.setBoundariesSquare(M_Bulk->getMesh());
            break;
        case BoundaryAssignmentType::TAG_NAME:
            if (M_Bulk->hasExternalMesh()) {
                std::cout << "[DarcyProblemT] Using Gmsh physical tag names for boundary conditions..." << std::endl;
                M_BC.setBoundariesFromTagsName(M_Bulk->getMesh(), M_Bulk->getRegionMap());
            } else {
                std::cout << "[DarcyProblemT] WARNING: TAG_NAME requested but no external mesh. "
                          << "Falling back to geometric detection." << std::endl;
                M_BC.setBoundariesCylinder(M_Bulk->getMesh());
            }
            break;
            
        case BoundaryAssignmentType::TAG_NUMBER:
            if (M_Bulk->hasExternalMesh()) {
                std::cout << "[DarcyProblemT] Using Gmsh tag numbers for boundary conditions..." << std::endl;
                M_BC.setBoundariesFromTagNumbersDirect(M_Bulk->getMesh(), M_tagNumberLargest);
            } else {
                std::cout << "[DarcyProblemT] WARNING: TAG_NUMBER requested but no external mesh. "
                          << "Falling back to geometric detection." << std::endl;
                M_BC.setBoundariesCylinder(M_Bulk->getMesh());
            }
            break;
    }
}

FEM* DarcyProblemT::getFEM(const std::string variable)
{
    if (variable=="Pressure")
    {
        return &M_PressureFEM;
    }	
    if (variable=="Velocity")
    {
        return &M_VelocityFEM;
    }
    std::cerr << "ERROR: Unknown variable '" << variable << "' in DarcyProblemT::getFEM" << std::endl;
    return &M_PressureFEM;  
}

size_type DarcyProblemT::getNDOF(std::string variable)
{
    if (variable=="Pressure")
    {
        return M_PressureFEM.nb_dof();
    }

    if (variable=="Velocity")
    {
        return M_VelocityFEM.nb_dof();
    }

    if (variable=="all")
    {
        return M_VelocityFEM.nb_dof()+M_PressureFEM.nb_dof();
    }
    std::cerr << "ERROR: Unknown variable '" << variable << "' in DarcyProblemT::getNDOF" << std::endl;
    return 0;
}

void DarcyProblemT::addToSys(LinearSystem* sys)
{
    M_Sys=sys;
    M_Sys->addToMatrix(M_nbTotDOF);
}

void DarcyProblemT::initialize()
{
    M_pressureSolOld.reset(new scalarVector_Type (M_PressureFEM.nb_dof()));
    M_pressureSol.reset(new scalarVector_Type (M_PressureFEM.nb_dof()));
    
    gmm::clear(*M_pressureSol);
    
    M_velocitySol.reset(new scalarVector_Type (M_VelocityFEM.nb_dof()));
    
    gmm::clear(*M_velocitySol);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (size_type i=0; i<M_PressureFEM.nb_dof("base");++i)
    {
        base_node nodo(M_PressureFEM.point_of_basic_dof(i));
        (*M_pressureSol)[i]=M_Bulk->getDarcyData()->pIni(nodo);
        (*M_pressureSolOld)[i]=(*M_pressureSol)[i];
    }
}

void DarcyProblemT::updateSol()
{	
    gmm::clear(*M_pressureSolOld);
    gmm::copy(*M_pressureSol, *M_pressureSolOld);
}

//A11=∫ u·v K⁻¹
//A12=∫ div(u) φ
//A22=∫ φψ * (1/dt)
void DarcyProblemT::assembleMatrix()
{
    sparseMatrixPtr_Type A11;
    A11.reset(new sparseMatrix_Type (M_VelocityFEM.nb_dof(),M_VelocityFEM.nb_dof()));
    gmm::clear(*A11);
    sparseMatrixPtr_Type A12;
    A12.reset(new sparseMatrix_Type (M_VelocityFEM.nb_dof(),M_PressureFEM.nb_dof()));
    gmm::clear(*A12);

    if (!M_pressureMass) {
        M_pressureMass.reset(new sparseMatrix_Type(M_PressureFEM.nb_dof(), M_PressureFEM.nb_dof()));
        gmm::clear(*M_pressureMass);
        massL2Standard(M_pressureMass, M_PressureFEM, M_CoeffFEM, M_intMethod);
    }

    massHdiv( A11, M_Bulk, M_VelocityFEM, M_CoeffFEM, M_intMethod);
            
    divHdiv( A12, M_VelocityFEM, M_PressureFEM, M_intMethod);

    essentialWNitsche( A11, M_Bulk, &M_BC, M_VelocityFEM, M_CoeffFEM, M_intMethod);

    sparseMatrixPtr_Type A22;
    A22.reset(new sparseMatrix_Type(M_PressureFEM.nb_dof(), M_PressureFEM.nb_dof()));
    gmm::clear(*A22);
    
    scalar_type M_v = M_Bulk->getDarcyData()->M();
    scalar_type gamma = M_Bulk->getDarcyData()->getLeakage();
    if (M_v<=0.0){
        std::cerr << "ERROR: Biot modulus must be positive. M_v = " 
                  << M_v << std::endl;
        throw std::runtime_error("Non-positive Biot modulus");
    }
    
    // Compute: (1/(M_v*Δt) + γ) * M_standard
    gmm::copy(*M_pressureMass, *A22);
    scalar_type scaling = 1.0/(M_v * M_dt) + gamma;
    gmm::scale(*A22, scaling);

    M_Sys->addSubMatrix(A11, 0,0);
        
    M_Sys->addSubMatrix(A12, 0, gmm::mat_ncols(*A11));
    M_Sys->addSubMatrix(A12, gmm::mat_nrows(*A11),0, -1.0, true);

    M_Sys->addSubMatrix(A22, M_VelocityFEM.nb_dof(),M_VelocityFEM.nb_dof());

    M_Sys->saveMatrix("mat_darcy.mm");
}

void DarcyProblemT::assembleRHS()
{
    scalarVectorPtr_Type  source;
    source.reset(new scalarVector_Type (M_PressureFEM.getFEM()->nb_dof()));
    gmm::clear(*source);

    scalarSource(source, M_Bulk, M_PressureFEM, M_CoeffFEM, M_intMethod, M_timeLoop->time() );
        
    M_Sys->addSubVector(source,M_VelocityFEM.nb_dof());
        
    scalarVectorPtr_Type  BCvec;
    BCvec.reset(new scalarVector_Type (M_VelocityFEM.getFEM()->nb_dof()));
    gmm::clear(*BCvec);
    essentialWNitscheRHS( BCvec, M_Bulk, &M_BC, M_VelocityFEM, M_CoeffFEM, M_intMethod);
    M_Sys->addSubVector(BCvec, 0);		
        
    gmm::clear(*BCvec);
        
    naturalRHS( BCvec, M_Bulk, &M_BC, M_VelocityFEM, M_CoeffFEM, M_intMethod, M_timeLoop->time() );

    M_Sys->addSubVector(BCvec, 0);	

    gmm::clear(*BCvec);

    vectorSource( BCvec, M_Bulk, M_VelocityFEM, M_CoeffFEM, M_intMethod);
        
    M_Sys->addSubVector(BCvec, 0);

    if (M_Bulk->getDarcyData()->M() > 0) {
        scalar_type scaling = 1.0 / (M_Bulk->getDarcyData()->M() * M_dt);
        M_Sys->multAddToRHS(M_pressureMass, M_pressureSolOld, 0,
                            M_VelocityFEM.nb_dof(), scaling);
    }
}

void DarcyProblemT::solve()
{
    M_Sys->solve();
    M_pressureSol.reset(new scalarVector_Type (M_PressureFEM.nb_dof()));
    gmm::clear(*M_pressureSol);
    M_Sys->extractSubVector(M_pressureSol, M_VelocityFEM.nb_dof(), "sol");

    gmm::clear(*M_velocitySol);
    M_Sys->extractSubVector(M_velocitySol, 0, "sol");

    if (M_timeLoop->time()==0)
    {
        M_pressureSolIni.reset(new scalarVector_Type (M_PressureFEM.nb_dof()));
        M_Sys->extractSubVector(M_pressureSolIni, M_VelocityFEM.nb_dof(), "sol");
    }
}

void DarcyProblemT::extractSol(scalarVectorPtr_Type sol)
{
    M_pressureSol.reset(new scalarVector_Type (M_PressureFEM.nb_dof()));
    gmm::clear(*M_pressureSol);
    gmm::copy ( gmm::sub_vector (*sol, gmm::sub_interval ( M_VelocityFEM.nb_dof(),  M_PressureFEM.nb_dof()) ), *M_pressureSol);

    M_velocitySol.reset(new scalarVector_Type (M_VelocityFEM.nb_dof()));
    gmm::clear(*M_velocitySol);
    gmm::copy ( gmm::sub_vector (*sol, gmm::sub_interval ( 0, M_VelocityFEM.nb_dof()) ), *M_velocitySol);
}

scalar_type DarcyProblemT::computeError(std::string what, scalar_type time)
{
    if (M_Bulk->getDarcyData()->hasExactSolution()){
        std::cout << "compute error"<<std::endl;
        scalar_type error;
        scalarVector_Type loc_err(M_PressureFEM.getFEM()->nb_dof());
        if (what=="Pressure" || what=="all")
        {

            for (size_type i=0;i<loc_err.size();++i)
            {
                bgeot::base_node where=M_PressureFEM.getFEM()->point_of_basic_dof(i);
                loc_err[i]=gmm::abs(((*M_pressureSol)[i]) - M_Bulk->getDarcyData()->pEx(where,time) );
            }
        }	
        error=L2Norm ( M_pressureMass,loc_err, M_PressureFEM, M_intMethod);

        std::cout <<"at time   "<<time<< "   error pressure   "<<error<<std::endl; 
        return error;
    }
    else return -1;
}

void DarcyProblemT::exportVtk(std::string folder, std::string what, int frame)
{
    std::cout << "export Darcy solution"<<std::endl;
    std::string frameNum;
    
    if (frame<0)
    {
        frameNum="";
    }
    else
    {
        if (frame<10)
        {
            frameNum.append("00"+LifeV::number2string(frame));
        }
        else
        {
            if (frame<100)
            {
                frameNum.append("0"+LifeV::number2string(frame));
            }
            else
            {
                frameNum.append(LifeV::number2string(frame));
            }
        }
    }		

    if (what=="Pressure" || what=="all")
    {
        getfem::vtk_export exp(folder + "Pressure_bulk"+ frameNum + ".vtk" );
        getfem::vtk_export expE(folder + "Pressure_ex"+ frameNum + ".vtk" );
        getfem::vtk_export expD(folder + "dP"+ frameNum + ".vtk" );
        
        exp.exporting( *(M_PressureFEM.getFEM()));
        exp.write_mesh();
        expE.exporting( *(M_PressureFEM.getFEM()));
        expE.write_mesh();
        expD.exporting( *(M_PressureFEM.getFEM()));
        expD.write_mesh();
        std::vector<scalar_type> diff(M_PressureFEM.getFEM()->nb_dof());
        std::vector<scalar_type> pe(M_PressureFEM.getFEM()->nb_dof());
        for (size_type i=0; i<diff.size(); ++i)
        {
            bgeot::base_node where=M_PressureFEM.getFEM()->point_of_basic_dof(i);
            pe[i]=M_Bulk->getDarcyData()->pEx(where,M_timeLoop->time() );
        }
        if (M_pressureSolIni!=NULL)
        {
            for (size_type i=0; i<diff.size(); ++i)
            {
                diff[i]=(*M_pressureSol)[i] - (*M_pressureSolIni)[i] ;
            }
        }
        if (true)
        {
            exp.write_point_data( *(M_PressureFEM.getFEM()), *M_pressureSol, "p");
            expE.write_point_data( *(M_PressureFEM.getFEM()), pe, "pex");
            expD.write_point_data( *(M_PressureFEM.getFEM()), diff, "pIni");
        }
        else
        {
            exp.write_point_data( *(M_PressureFEM.getFEM()), *M_pressureSolOld, "p");
        }
    }
    
    if (what=="Velocity" || what=="all")
    {
        getfem::vtk_export  exp(folder + "Velocity_bulk"+ frameNum + ".vtk" );
        
        std::vector<scalar_type> velInt(M_visualizationVFEM.getFEM()->nb_dof(),0.0);
        getfem::interpolation(*(M_VelocityFEM.getFEM()), *(M_visualizationVFEM.getFEM()),*M_velocitySol, velInt );
        exp.exporting(*(M_visualizationVFEM.getFEM()));
        exp.write_mesh();
        exp.write_point_data( *(M_visualizationVFEM.getFEM()), velInt, "v");
    }
}