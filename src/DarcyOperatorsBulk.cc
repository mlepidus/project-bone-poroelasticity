#include "../include/DarcyOperatorsBulk.h"


void massL2Standard(sparseMatrixPtr_Type M,  FEM& femP, 
                    FEM& femC, getfem::mesh_im& im)
{
    // Same as massL2 but with dt=1 and leakage=0
    getfem::mesh_fem femP_mf = *(femP.getFEM());
    getfem::mesh_fem femC_mf = *(femC.getFEM());
    
    std::vector<scalar_type> ones(femC_mf.nb_dof(), 1.0);
    
    getfem::generic_assembly assem;
    assem.set("x=data(#2);"
              "a=comp(vBase(#1).vBase(#1).Base(#2));"
              "M(#1,#1)+=a(:,i,:,i,j).x(j)");
    
    assem.push_mi(im);
    assem.push_mf(femP_mf);
    assem.push_mf(femC_mf);
    assem.push_data(ones);  // Coefficient = 1 everywhere
    assem.push_mat(*M);
    assem.assembly(-1);
}


//without providing mass matrix, recompute it
scalar_type L2Norm ( scalarVectorPtr_Type V,  FEM& femP, getfem::mesh_im& im)
{
	
 	sparseMatrixPtr_Type M;
	M.reset(new sparseMatrix_Type (femP.nb_dof(),femP.nb_dof()));
			gmm::clear(*M);
	massL2Standard(M, femP,  femP, im);
	scalarVector_Type VV(V->size(),0);
	gmm::mult(*M, *V,VV);
	scalar_type norm=gmm::vect_sp(*V, VV);
	return pow(norm,0.5);

}


scalar_type L2Norm ( scalarVector_Type V, FEM& femP, getfem::mesh_im& im)
{
 	sparseMatrixPtr_Type M;
	M.reset(new sparseMatrix_Type (femP.nb_dof(),femP.nb_dof()));
			gmm::clear(*M);
		massL2Standard(M, femP,  femP, im);
	scalarVector_Type VV(V.size(),0);
	gmm::mult(*M, V,VV);
	scalar_type norm=gmm::vect_sp(V, VV);
	return pow(norm,0.5);

}


scalar_type L2Norm ( sparseMatrixPtr_Type M, scalarVector_Type V, FEM& femP, getfem::mesh_im& im)
{
	if (M==nullptr){
        sparseMatrixPtr_Type M;
        M.reset(new sparseMatrix_Type (femP.nb_dof(),femP.nb_dof()));
                gmm::clear(*M);
       	massL2Standard(M, femP,  femP, im);
    }
	scalarVector_Type VV(V.size(),0);
	gmm::mult(*M, V,VV);
	scalar_type norm=gmm::vect_sp(V, VV);
	return pow(norm,0.5);
}

void massHdiv(sparseMatrixPtr_Type M,
              Bulk* medium, FEM& FemV, FEM& FemP, getfem::mesh_im& im)
{
    getfem::mesh_fem femV(*(FemV.getFEM()));
    getfem::mesh_fem femP((*FemP.getFEM()));

    size_type dim = femV.linked_mesh().dim();
    
    // Prepare coefficient vectors
    std::vector<scalar_type> iKxx(femP.nb_dof(), 0.0);
    std::vector<scalar_type> iKyy(femP.nb_dof(), 0.0);
    std::vector<scalar_type> iKzz(femP.nb_dof(), 0.0);
    std::vector<scalar_type> iKxy(femP.nb_dof(), 0.0);
    std::vector<scalar_type> iKxz(femP.nb_dof(), 0.0);
    std::vector<scalar_type> iKyz(femP.nb_dof(), 0.0);


    // Compute inverse of permeability tensor at each DOF
    for (size_type i = 0; i < femP.nb_dof(); ++i)
    {   
        if (dim == 2) {
            scalar_type Kxx = medium->getDarcyData()->getKxx(i);
            scalar_type Kyy = medium->getDarcyData()->getKyy(i);
            scalar_type Kxy = medium->getDarcyData()->getKxy(i);
            
            scalar_type det = Kxx * Kyy - Kxy * Kxy;
            if (std::abs(det) < 1e-14) det = 1.0;
            
            iKxx[i] = Kyy / det;
            iKyy[i] = Kxx / det;
            iKxy[i] = -Kxy / det;
        }
        else if (dim == 3) {
            scalar_type Kxx = medium->getDarcyData()->getKxx(i);
            scalar_type Kyy = medium->getDarcyData()->getKyy(i);
            scalar_type Kzz = medium->getDarcyData()->getKzz(i);
            scalar_type Kxy = medium->getDarcyData()->getKxy(i);
            scalar_type Kxz = medium->getDarcyData()->getKxz(i);
            scalar_type Kyz = medium->getDarcyData()->getKyz(i);
            
            scalar_type det = Kxx*(Kyy*Kzz - Kyz*Kyz) 
                            - Kxy*(Kxy*Kzz - Kyz*Kxz) 
                            + Kxz*(Kxy*Kyz - Kyy*Kxz);
            
            if (std::abs(det) < 1e-14) {
                iKxx[i] = 1.0;
                iKyy[i] = 1.0;
                iKzz[i] = 1.0;
                iKxy[i] = 0.0;
                iKxz[i] = 0.0;
                iKyz[i] = 0.0;
                continue;
            }
            
            iKxx[i] = (Kyy*Kzz - Kyz*Kyz) / det;
            iKyy[i] = (Kxx*Kzz - Kxz*Kxz) / det;
            iKzz[i] = (Kxx*Kyy - Kxy*Kxy) / det;
            iKxy[i] = (Kxz*Kyz - Kxy*Kzz) / det;
            iKxz[i] = (Kxy*Kyz - Kxz*Kyy) / det;
            iKyz[i] = (Kxy*Kxz - Kxx*Kyz) / det;
        }
    }
    
    getfem::generic_assembly assem;

    if (dim == 2) {
        assem.set("xx=data(#2); yy=data$2(#2); xy=data$3(#2);"
                  "t=comp(vBase(#1).vBase(#1).Base(#2));"
                  "M(#1,#1)+=t(:,1,:,1,j).xx(j) + t(:,2,:,2,j).yy(j) + "
                  "t(:,1,:,2,j).xy(j) + t(:,2,:,1,j).xy(j);");
                  
        assem.push_mi(im);
        assem.push_mf(femV);
        assem.push_mf(femP);
        assem.push_data(iKxx);
        assem.push_data(iKyy);
        assem.push_data(iKxy);
        assem.push_mat(*M);
        assem.assembly(-1);
    }
    else if (dim == 3) {
        // For 3D, assemble diagonal terms first, then off-diagonal
        assem.set("xx=data(#2);"
                  "t=comp(vBase(#1).vBase(#1).Base(#2));"
                  "M(#1,#1)+=t(:,1,:,1,j).xx(j);");
        assem.push_mi(im);
        assem.push_mf(femV);
        assem.push_mf(femP);
        assem.push_data(iKxx);
        assem.push_mat(*M);
        assem.assembly(-1);
        
        // Assemble Kyy term
        getfem::generic_assembly assem_yy;
        assem_yy.set("yy=data(#2);"
                     "t=comp(vBase(#1).vBase(#1).Base(#2));"
                     "M(#1,#1)+=t(:,2,:,2,j).yy(j);");
        assem_yy.push_mi(im);
        assem_yy.push_mf(femV);
        assem_yy.push_mf(femP);
        assem_yy.push_data(iKyy);
        assem_yy.push_mat(*M);
        assem_yy.assembly(-1);
        
        // Assemble Kzz term
        getfem::generic_assembly assem_zz;
        assem_zz.set("zz=data(#2);"
                     "t=comp(vBase(#1).vBase(#1).Base(#2));"
                     "M(#1,#1)+=t(:,3,:,3,j).zz(j);");
        assem_zz.push_mi(im);
        assem_zz.push_mf(femV);
        assem_zz.push_mf(femP);
        assem_zz.push_data(iKzz);
        assem_zz.push_mat(*M);
        assem_zz.assembly(-1);
        
        // Assemble Kxy off-diagonal terms
        getfem::generic_assembly assem_xy;
        assem_xy.set("xy=data(#2);"
                     "t=comp(vBase(#1).vBase(#1).Base(#2));"
                     "M(#1,#1)+=t(:,1,:,2,j).xy(j) + t(:,2,:,1,j).xy(j);");
        assem_xy.push_mi(im);
        assem_xy.push_mf(femV);
        assem_xy.push_mf(femP);
        assem_xy.push_data(iKxy);
        assem_xy.push_mat(*M);
        assem_xy.assembly(-1);
        
        // Assemble Kxz off-diagonal terms
        getfem::generic_assembly assem_xz;
        assem_xz.set("xz=data(#2);"
                     "t=comp(vBase(#1).vBase(#1).Base(#2));"
                     "M(#1,#1)+=t(:,1,:,3,j).xz(j) + t(:,3,:,1,j).xz(j);");
        assem_xz.push_mi(im);
        assem_xz.push_mf(femV);
        assem_xz.push_mf(femP);
        assem_xz.push_data(iKxz);
        assem_xz.push_mat(*M);
        assem_xz.assembly(-1);
        
        // Assemble Kyz off-diagonal terms
        getfem::generic_assembly assem_yz;
        assem_yz.set("yz=data(#2);"
                     "t=comp(vBase(#1).vBase(#1).Base(#2));"
                     "M(#1,#1)+=t(:,2,:,3,j).yz(j) + t(:,3,:,2,j).yz(j);");
        assem_yz.push_mi(im);
        assem_yz.push_mf(femV);
        assem_yz.push_mf(femP);
        assem_yz.push_data(iKyz);
        assem_yz.push_mat(*M);
        assem_yz.assembly(-1);
    }
    else {
        throw std::runtime_error("massHdiv: unsupported dimension");
    }

    std::cout << "[DARCY] operator a(volume) assembled - " 
              << gmm::nnz(*M) << " non-zeros" << std::endl;
}


void divHdiv( sparseMatrixPtr_Type M, FEM& FemV, FEM& FemP, getfem::mesh_im& im)
{

    getfem::mesh_fem femV(*(FemV.getFEM()));
    getfem::mesh_fem femP((*FemP.getFEM()));

   
    // H(div) mass matrix: integral of (K^{-1} u) · v over the domain

    getfem::generic_assembly assem;

     assem.set("M(#1,#2)+=-comp(vGrad(#1).Base(#2))"
        "(:,i,i, :);");
                
    // Assign the mesh integration method
    assem.push_mi(im);
    // Assign the mesh finite element space
    assem.push_mf(femV);
    // Assign the mesh finite element space
    
    // Assign the mesh finite element space for the coefficients
    assem.push_mf(femP);


    // Set the matrices to save the evaluations
    assem.push_mat(*M);

    // Computes the matrices
    assem.assembly(-1);


    std::cout << "[DARCY] operator b(volume) assembled" << std::endl;

}

void scalarSource (scalarVectorPtr_Type V, Bulk* medium,  FEM& FemP, FEM& FemC, getfem::mesh_im& im, const scalar_type t=0)
{

	  getfem::mesh_fem femP(*(FemP.getFEM()));
          getfem::mesh_fem femC((*FemC.getFEM()));
    
	  getfem::generic_assembly assem;
	  assem.set("w=data(#2);"
		       "a=comp(Base(#1).Base(#2));"
		       "V(#1)+=a(:, k).w(k)");
	  std::vector<scalar_type>  sourceVal(femC.nb_dof(),0.0);
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
	  for (size_type i=0; i<femC.nb_dof();++i)
	  {
	  	sourceVal[i]=medium->getDarcyData()->source(femC.point_of_basic_dof(i),t);
	  }  	     
	    
          // Assign the mesh integration method
    	  assem.push_mi(im);

	  // Assign the mesh finite element space
          assem.push_mf(femP);
          // Assign the mesh finite element space for the coefficients
    	  assem.push_mf(femC);
	  assem.push_data(sourceVal);
	  
          // Set the matrices to save the evaluations
          assem.push_vec(*V);

          // Computes the matrices
          assem.assembly(-1);      

}

void vectorSource (scalarVectorPtr_Type V, Bulk* medium,  FEM& FemV, FEM& FemC, getfem::mesh_im& im)
{

	getfem::mesh_fem femV(*(FemV.getFEM()));
    getfem::mesh_fem femC((*FemC.getFEM()));
    
    size_type dim = medium->getDim();

    std::vector<scalar_type> sourceValX(femC.nb_dof(),0.0);
    std::vector<scalar_type> sourceValY(femC.nb_dof(),0.0);
    std::vector<scalar_type> sourceValZ;

	getfem::generic_assembly assem;
    if (dim == 2) {
        assem.set("vx=data$1(#2);" "vy=data$2(#2);"
                  "a=comp(vBase(#1).Base(#2));"
                  "V(#1)+=a(:,1,k).vx(k)+a(:,2,k).vy(k)");
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (size_type i=0; i<femC.nb_dof();++i) {
            sourceValX[i]=medium->getDarcyData()->rhoF()*medium->getDarcyData()->gravity()[0];
            sourceValY[i]=medium->getDarcyData()->rhoF()*medium->getDarcyData()->gravity()[1];
        }
        
        assem.push_mi(im);
        assem.push_mf(femV);
        assem.push_mf(femC);
        assem.push_data(sourceValX);
        assem.push_data(sourceValY);
        assem.push_vec(*V);
        assem.assembly(-1);
    }
    else if (dim == 3) {
        assem.set("vx=data$1(#2);" "vy=data$2(#2);" "vz=data$3(#2);"
                  "a=comp(vBase(#1).Base(#2));"
                  "V(#1)+=a(:,1,k).vx(k)+a(:,2,k).vy(k)+a(:,3,k).vz(k)");
        
        sourceValZ.resize(femC.nb_dof(), 0.0);  
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (size_type i=0; i<femC.nb_dof();++i) {
            sourceValX[i]=medium->getDarcyData()->rhoF()*medium->getDarcyData()->gravity()[0];
            sourceValY[i]=medium->getDarcyData()->rhoF()*medium->getDarcyData()->gravity()[1];
            sourceValZ[i]=medium->getDarcyData()->rhoF()*medium->getDarcyData()->gravity()[2];  
        }
        
        assem.push_mi(im);
        assem.push_mf(femV);
        assem.push_mf(femC);
        assem.push_data(sourceValX);
        assem.push_data(sourceValY);
        assem.push_data(sourceValZ);  
        assem.push_vec(*V);
        assem.assembly(-1);
    }

	else {
	throw std::runtime_error("vectorSource: unsupported dimension " + std::to_string(dim));
	  }
  
}

void vectorSource (scalarVectorPtr_Type V, FEM& FemV, FEM& FemC, scalarVectorPtr_Type data, getfem::mesh_im& im)
{

	  getfem::mesh_fem femV(*(FemV.getFEM()));
          getfem::mesh_fem femC((*FemC.getFEM()));
    
	  getfem::generic_assembly assem;
	  assem.set("v=data$1(#2);" 
		       "a=comp(vBase(#1).vBase(#2));"
		       "V(#1)+=a(:,i,k,i).v(k)");

	 

          // Assign the mesh integration method
    	  assem.push_mi(im);

	  // Assign the mesh finite element space
          assem.push_mf(femV);
          // Assign the mesh finite element space for the coefficients
    	  assem.push_mf(femC);
	  assem.push_data(*data);
	  
          // Set the matrices to save the evaluations
          assem.push_vec(*V);

          // Computes the matrices
          assem.assembly(-1);      

}

