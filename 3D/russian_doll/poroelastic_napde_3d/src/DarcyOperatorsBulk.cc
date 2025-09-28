#include "../include/DarcyOperatorsBulk.h"

void massL2( sparseMatrixPtr_Type M,
               Bulk* medium, FEM& FemP, FEM& FemC, getfem::mesh_im& im, scalar_type dt)
{

    getfem::mesh_fem femP(*(FemP.getFEM()));
    getfem::mesh_fem femC((*FemC.getFEM()));

   
    massL2(  M,medium,  femP,  femC, im,  dt,-1);
}

void massL2Local( sparseMatrixPtr_Type M,
               Bulk* medium, FEM& FemP, FEM& FemC, getfem::mesh_im& im, scalar_type dt, size_type icv)
{

    getfem::mesh_fem femP(*(FemP.getFEM()));
    getfem::mesh_fem femC((*FemC.getFEM()));

    sparseMatrix_Type M_(femP.nb_dof(),femP.nb_dof());
    
    //questo è il termine tau tau (matrice di "massa" per le velocità)

    getfem::generic_assembly assem;

    assem.set ( "x=data(#2);"
                "a=comp(vBase(#1).vBase(#1).Base(#2));"
                "M(#1,#1)+=a(:,i,:,i,j).x(j)");
      
                
    // Assign the mesh integration method
    assem.push_mi(im);
    // Assign the mesh finite element space
    assem.push_mf(femP);
    // Assign the mesh finite element space
    
    // Assign the mesh finite element space for the coefficients
    assem.push_mf(femC);
    std::vector<scalar_type> iKxx(femC.nb_dof(),0.0);
    
    scalar_type leakage = medium->getDarcyData()->getLeakage();
    size_type ii=femC.ind_basic_dof_of_element(icv)[0];
    	iKxx[ii]=1./dt;
   
    
    
    // Assign the coefficients
    
    assem.push_data(iKxx);

    // Set the matrices to save the evaluations
    assem.push_mat(M_);

    // Computes the matrices
    assem.assembly(-1);

   for ( size_type i = 0; i < femP.nb_dof(); ++i )
    {
        for ( size_type j = 0; j < femP.nb_dof(); ++j )
        {
            (*M)(i, j) += M_(i, j);
	    
        }
    }

    std::cout << "DARCY :: operator a(volume)      [OK]" << std::endl;

}

void massL2( sparseMatrixPtr_Type M,
               Bulk* medium, getfem::mesh_fem& femP, getfem::mesh_fem& femC, getfem::mesh_im& im, scalar_type dt, int region)
{

   
    //questo è il termine tau tau (matrice di "massa" per le velocità)

    getfem::generic_assembly assem;

       assem.set ( "x=data(#2);"
                "a=comp(vBase(#1).vBase(#1).Base(#2));"
                "M(#1,#1)+=a(:,i,:,i,j).x(j)");
                
    // Assign the mesh integration method
    assem.push_mi(im);
    // Assign the mesh finite element space
    assem.push_mf(femP);
    // Assign the mesh finite element space
    
    // Assign the mesh finite element space for the coefficients
    assem.push_mf(femC);
    std::vector<scalar_type> iKxx(femC.nb_dof(),0.0);
    
    scalar_type leakage = medium->getDarcyData()->getLeakage();
    for (int i=0; i<femC.nb_dof();++i)
    {   
    	iKxx[i]=1./dt + leakage;
   
    }
    
    // Assign the coefficients
    
    assem.push_data(iKxx);

    // Set the matrices to save the evaluations
    assem.push_mat(*M);

    // Computes the matrices
    assem.assembly(region);

    std::cout << "DARCY :: operator mass(volume)      [OK]" << std::endl;

}

scalar_type L2Norm ( scalarVectorPtr_Type V, Bulk* medium, FEM& femP, getfem::mesh_im& im)
{
	
 	sparseMatrixPtr_Type M;
	M.reset(new sparseMatrix_Type (femP.nb_dof(),femP.nb_dof()));
			gmm::clear(*M);
	massL2(M,  medium, femP,  femP, im, 1);
	scalarVector_Type VV(V->size(),0);
	gmm::mult(*M, *V,VV);
	scalar_type norm=gmm::vect_sp(*V, VV);
	return pow(norm,0.5);

}

scalar_type L2Norm ( scalarVector_Type V, Bulk* medium, FEM& femP, getfem::mesh_im& im)
{
	
 	sparseMatrixPtr_Type M;
	M.reset(new sparseMatrix_Type (femP.nb_dof(),femP.nb_dof()));
			gmm::clear(*M);
	massL2(M,  medium, femP,  femP, im, 1);
	scalarVector_Type VV(V.size(),0);
	gmm::mult(*M, V,VV);
	scalar_type norm=gmm::vect_sp(V, VV);
	return pow(norm,0.5);

}
scalar_type L2Norm ( scalarVector_Type V, Bulk* medium, FEM& femV,FEM& femP, getfem::mesh_im& im)
{
	
 	sparseMatrixPtr_Type M;
	M.reset(new sparseMatrix_Type (femV.nb_dof(),femV.nb_dof()));
			gmm::clear(*M);
	massL2(M,  medium, femV,  femP, im, 1);
	scalarVector_Type VV(V.size(),0);
	gmm::mult(*M, V,VV);
	scalar_type norm=gmm::vect_sp(V, VV);
	return pow(norm,0.5);

}

scalar_type L2Norm ( sparseMatrixPtr_Type M, scalarVector_Type V, Bulk* medium, getfem::mesh_fem& femP, getfem::mesh_im& im, int region)
{
	
 	/*sparseMatrixPtr_Type M;
	M.reset(new sparseMatrix_Type (femP.nb_dof(),femP.nb_dof()));
			gmm::clear(*M);*/
//	massL2(M,  medium, femP,  femP, im, 1, region);
	scalarVector_Type VV(V.size(),0);
	gmm::mult(*M, V,VV);
	scalar_type norm=gmm::vect_sp(V, VV);
	return pow(norm,0.5);

}

void massHdiv( sparseMatrixPtr_Type M,
               Bulk* medium, FEM& FemV, FEM& FemP, getfem::mesh_im& im)
{

    getfem::mesh_fem femV(*(FemV.getFEM()));
    getfem::mesh_fem femP((*FemP.getFEM()));

    
    //questo è il termine tau tau (matrice di "massa" per le velocità)

    getfem::generic_assembly assem;

    assem.set ( "xx=data(#2);""yy=data$2(#2);""xy=data$3(#2);""xz=data$4(#2);""yz=data$5(#2);""zz=data$6(#2);"
                "a=comp(vBase(#1).vBase(#1).Base(#2));"
                "M(#1,#1)+=a(:,1,:,1,j).xx(j)+a(:,2,:,2,j).yy(j)+a(:,3,:,3,j).zz(j)+a(:,1,:,2,j).xy(j)+a(:,2,:,1,j).xy(j)+ a(:,1,:,3,j).xz(j)+a(:,3,:,1,j).xz(j)+a(:,2,:,3,j).yz(j)+a(:,3,:,2,j).yz(j);");
      
                
    // Assign the mesh integration method
    assem.push_mi(im);
    // Assign the mesh finite element space
    assem.push_mf(femV);
    // Assign the mesh finite element space
    
    // Assign the mesh finite element space for the coefficients
    assem.push_mf(femP);
    std::vector<scalar_type> iKxx(femP.nb_dof(),0.0);
    std::vector<scalar_type> iKxy(femP.nb_dof(),0.0);
    std::vector<scalar_type> iKyy(femP.nb_dof(),0.0);
    std::vector<scalar_type> iKxz(femP.nb_dof(),0.0);
    std::vector<scalar_type> iKyz(femP.nb_dof(),0.0);
    std::vector<scalar_type> iKzz(femP.nb_dof(),0.0);

    for (int i=0; i<femP.nb_dof();++i)
    {   
        scalar_type xx=medium->getDarcyData()->getKxx(i);
        scalar_type yy=medium->getDarcyData()->getKyy(i);
        scalar_type zz=medium->getDarcyData()->getKzz(i);
        scalar_type xy=medium->getDarcyData()->getKxy(i);
        scalar_type xz=medium->getDarcyData()->getKxz(i);
        scalar_type yz=medium->getDarcyData()->getKyz(i);
        scalar_type det= -zz*xy*xy + 2*xy*xz*yz - yy*xz*xz - xx*yz*yz + xx*yy*zz;
    	iKxx[i]=(- yz*yz + yy*zz)/det;
        iKxy[i]=(xz*yz - xy*zz)/det;
        iKxz[i]=(xy*yz - xz*yy)/det;
        iKyy[i]=(- xz*xz + xx*zz)/det;
        iKyz[i]=(xy*xz - xx*yz)/det;
        iKzz[i]=(- xy*xy + xx*yy)/det;
    }
    
    // Assign the coefficients
    
    assem.push_data(iKxx);
    assem.push_data(iKyy);
    assem.push_data(iKxy);
    assem.push_data(iKxy);
    assem.push_data(iKyz);
    assem.push_data(iKzz);

    // Set the matrices to save the evaluations
    assem.push_mat(*M);

    // Computes the matrices
    assem.assembly(-1);

    std::cout << "DARCY :: operator a(volume)      [OK]" << std::endl;

}

void divHdiv( sparseMatrixPtr_Type M,
               Bulk* medium, FEM& FemV, FEM& FemP, getfem::mesh_im& im)
{

    getfem::mesh_fem femV(*(FemV.getFEM()));
    getfem::mesh_fem femP((*FemP.getFEM()));

   
    //questo è il termine tau tau (matrice di "massa" per le velocità)

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


    std::cout << "DARCY :: operator b(volume)      [OK]" << std::endl;

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
    
	  getfem::generic_assembly assem;
	  assem.set("vx=data$1(#2);" "vy=data$2(#2);"
		       "a=comp(vBase(#1).Base(#2));"
		       "V(#1)+=a(:,1,k).vx(k)+a(:,2,k).vy(k)");

	  std::vector<scalar_type>  sourceValX(femC.nb_dof(),0.0);
	  std::vector<scalar_type>  sourceValY(femC.nb_dof(),0.0);
	
	  for (size_type i=0; i<femC.nb_dof();++i)
	  {
	  	sourceValX[i]=medium->getDarcyData()->rhoF()*medium->getDarcyData()->gravity()[0];
	  	sourceValY[i]=medium->getDarcyData()->rhoF()*medium->getDarcyData()->gravity()[1];
	  }  	     

          // Assign the mesh integration method
    	  assem.push_mi(im);

	  // Assign the mesh finite element space
          assem.push_mf(femV);
          // Assign the mesh finite element space for the coefficients
    	  assem.push_mf(femC);
	  assem.push_data(sourceValX);
	  assem.push_data(sourceValY);
	  
          // Set the matrices to save the evaluations
          assem.push_vec(*V);

          // Computes the matrices
          assem.assembly(-1);      

}

void vectorSource (scalarVectorPtr_Type V, Bulk* medium,  FEM& FemV, FEM& FemC, scalarVectorPtr_Type data, getfem::mesh_im& im)
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

