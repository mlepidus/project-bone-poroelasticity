#include "../include/ElastOperatorsBulk.h"
#include <stdexcept>

void stiffElast( sparseMatrixPtr_Type M,
               Bulk* medium, FEM& FemD, FEM& FemC, getfem::mesh_im& im)
{
    getfem::mesh_fem femD(*(FemD.getFEM()));
    getfem::mesh_fem femC((*FemC.getFEM()));

   
    getfem::generic_assembly assem;
    
    assem.set("lambda=data$1(#2); mu=data$2(#2);"
			   "t=comp(vGrad(#1).vGrad(#1).Base(#2));"
			   "M(#1,#1)+= sym(t(:,i,j,:,i,j,k).mu(k)"
			   "+ t(:,j,i,:,i,j,k).mu(k)"
			   "+ t(:,i,i,:,j,j,k).lambda(k))");
     // Assign the mesh integration method
    assem.push_mi(im);
    // Assign the mesh finite element space
    assem.push_mf(femD);
 
    // Assign the mesh finite element space for the coefficients
    assem.push_mf(femC);
    
    std::vector<scalar_type> lambda(femC.nb_dof(),0.0);
    std::vector<scalar_type> mu(femC.nb_dof(),0.0);
   
    for (size_type i=0; i<femC.nb_dof();++i)
    {   
    	lambda[i]=medium->getElastData()->getLambda(i);
    	mu[i]=medium->getElastData()->getMu(i); 
    }

    assem.push_data(lambda);

    assem.push_data(mu);
    
    // Set the matrices to save the evaluations
    assem.push_mat(*M);

    assem.assembly(-1);
  

    std::cout << "elast :: operator a(volume)      [OK]" << std::endl;

}

void bulkLoad( scalarVectorPtr_Type V,
               Bulk* medium, FEM& FemD, FEM& FemC, getfem::mesh_im& im, scalar_type time)
{
    getfem::mesh_fem femD(*(FemD.getFEM()));
    getfem::mesh_fem femC((*FemC.getFEM()));
    
    size_type dim = femD.linked_mesh().dim();
    
    getfem::generic_assembly assem;
    
    if (dim == 2) {
        assem.set("datax=data$1(#2);""datay=data$2(#2);"
            "t=comp(vBase(#1).Base(#2));"
            "V(#1)+=t(:,1,k).datax(k)+t(:,2,k).datay(k);");
    }
    else if (dim == 3) {
        assem.set("datax=data$1(#2);""datay=data$2(#2);""dataz=data$3(#2);"
            "t=comp(vBase(#1).Base(#2));"
            "V(#1)+=t(:,1,k).datax(k)+t(:,2,k).datay(k)+t(:,3,k).dataz(k);");
    }
    else {
        throw std::runtime_error("bulkLoad: unsupported dimension " + std::to_string(dim));
    }
        
    // Assign the mesh integration method
    assem.push_mi(im);
    // Assign the mesh finite element space
    assem.push_mf(femD);
 
    // Assign the mesh finite element space for the coefficients
    assem.push_mf(femC);
    
    // Declare all data vectors at same scope to ensure they remain valid until assembly
    scalarVector_Type datax(femC.nb_dof());
    scalarVector_Type datay(femC.nb_dof());
    scalarVector_Type dataz;  // Only populated for 3D
    
    for (size_type i=0; i<femC.nb_dof();++i)
    {
    	datax [ i ] = medium->getElastData()->bulkLoad(femC.point_of_basic_dof(i),time)[0];
    	datay [ i ] = medium->getElastData()->bulkLoad(femC.point_of_basic_dof(i),time)[1];
    }
    
    assem.push_data(datax);
    assem.push_data(datay);
    
    if (dim == 3) {
        dataz.resize(femC.nb_dof());
        for (size_type i=0; i<femC.nb_dof();++i)
        {
            dataz [ i ] = medium->getElastData()->bulkLoad(femC.point_of_basic_dof(i),time)[2];
        }
        assem.push_data(dataz);
    }
    
    // Set the matrices to save the evaluations
    assem.push_vec(*V);
    assem.assembly(-1);

    std::cout << "elast :: bulk load (" << dim << "D)     [OK]" << std::endl;
}


void givenFluidP( scalarVectorPtr_Type V,
               Bulk* medium, FEM& FemD, FEM& FemC, getfem::mesh_im& im)
{
    getfem::mesh_fem femD(*(FemD.getFEM()));
    getfem::mesh_fem femC((*FemC.getFEM()));
    
    
    getfem::generic_assembly assem;
    
    assem.set("datax=data(#2);"
        "t=comp(vGrad(#1).Base(#2));"
        "V(#1)+=t(:,i,i,k).datax(k);");
    
      // Assign the mesh integration method
    assem.push_mi(im);
    
    // Assign the mesh finite element space
    assem.push_mf(femD);
    
    // Assign the mesh finite element space for the coefficients
    assem.push_mf(femC);
    
    scalarVector_Type datax(femC.nb_dof());
    
    for (size_type i=0; i<femC.nb_dof();++i)
    {
    	datax [ i ] = medium->getElastData()->fluidP(femC.point_of_basic_dof(i));
    }
  
    assem.push_data(datax);
    
    // Set the matrices to save the evaluations
    assem.push_vec(*V);
    assem.assembly(-1); 

    

}


void givenFluidP( scalarVectorPtr_Type V, scalarVectorPtr_Type pressure,
               FEM& FemD, FEM& FemC, getfem::mesh_im& im)
{
    getfem::mesh_fem femD(*(FemD.getFEM()));
    getfem::mesh_fem femC((*FemC.getFEM()));
    
   
    getfem::generic_assembly assem;
    
    assem.set("datax=data(#2);"
        "t=comp(vGrad(#1).Base(#2));"
        "V(#1)+=t(:,i,i,k).datax(k);");
    
      // Assign the mesh integration method
    assem.push_mi(im);
    
    // Assign the mesh finite element space
    assem.push_mf(femD);
    
    // Assign the mesh finite element space for the coefficients
    assem.push_mf(femC);
    
     assem.push_data(*pressure);
    
    // Set the matrices to save the evaluations
    assem.push_vec(*V);
    assem.assembly(-1); 

}

void massMatrix( sparseMatrixPtr_Type M, 
                FEM& FemD, getfem::mesh_im& im)
{
    getfem::mesh_fem femD(*(FemD.getFEM()));
        
    getfem::generic_assembly assem;

    if (femD.get_qdim()>1)
    {    
   	 assem.set("t=comp(vBase(#1).vBase(#1));"
        "M(#1,#1)+=t(:,i,:,i);");
    }
    else
    {
	 assem.set("t=comp(Base(#1).Base(#1));"
        "M(#1,#1)+=t(:,:);");
    }
      // Assign the mesh integration method
    assem.push_mi(im);
    
    // Assign the mesh finite element space
    assem.push_mf(femD);
     
     
    // Set the matrices to save the evaluations
    assem.push_mat(*M);
    assem.assembly(-1); 
   				
 
}

void matrixFluidP( sparseMatrixPtr_Type M, 
                FEM& FemD, FEM& FemS, getfem::mesh_im& im)
{
    getfem::mesh_fem femD(*(FemD.getFEM()));
    getfem::mesh_fem femS((*FemS.getFEM()));
    
   
    getfem::generic_assembly assem;
    
    assem.set("t=-comp(vGrad(#1).Base(#2));"
        "M(#1,#2)+=t(:,i,i,:);");
    
      // Assign the mesh integration method
    assem.push_mi(im);
    
    // Assign the mesh finite element space
    assem.push_mf(femD);
    assem.push_mf(femS);
     
     
    // Set the matrices to save the evaluations
    assem.push_mat(*M);
    assem.assembly(-1); 
  
}

scalar_type L2Norm_Elast (sparseMatrixPtr_Type M, scalarVector_Type V)
{
    scalarVector_Type VV(V.size(), 0);
    gmm::mult(*M, V, VV);
    scalar_type norm = gmm::vect_sp(V, VV);
    return pow(norm, 0.5);
}