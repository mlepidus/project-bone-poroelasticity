#include "../include/ElastOperatorsBD.h"

void essentialWNitscheVec( sparseMatrixPtr_Type M,
               Bulk* medium, BC* bcPtr,  FEM& FemD, FEM& FemC, getfem::mesh_im& im)
               
{
    getfem::mesh_fem femD(*(FemD.getFEM()));
    getfem::mesh_fem femC((*FemC.getFEM()));
    
    scalarVector_Type etaGammaUinvh(femC.nb_dof());
    std::vector<scalar_type> lambda(femC.nb_dof(),0.0);
    std::vector<scalar_type> mu(femC.nb_dof(),0.0);

  for ( size_type i = 0; i < femC.nb_dof(); ++i )
    {
        scalar_type h_el= medium->getMesh()->convex_radius_estimate(femC.first_convex_of_basic_dof(i));
        scalar_type inv_h_mesh = 1.0 / h_el;
	lambda[i]=medium->getElastData()->getLambda(i);
    	mu[i]=medium->getElastData()->getMu(i);
        etaGammaUinvh [ i ] = 5*mu[i] *inv_h_mesh;
    }

 
    getfem::generic_assembly assem_surf, assem_surf2;
    getfem::generic_assembly assem_cons, assem_cons2;

    assem_surf.set("gamma=data$1(#2);"
        "t=comp(vBase(#1).vBase(#1).Base(#2));"
        "M$1(#1,#1)+=(t(:,i, :,i, k).gamma(k));");
    assem_cons.set("lambda=data$1(#2); mu=data$2(#2);"
			   "t=-comp(vBase(#1).vGrad(#1).Normal().Base(#2));"
			   "M(#1,#1)+=t(:,i,:,i,j,j,k).mu(k)+t(:,i,:,j,i,j,k).mu(k)+t(:,i,:,j,j,i,k).lambda(k)");

    assem_surf2.set("gamma=data$1(#2);"
        "t=comp(vBase(#1).Normal().vBase(#1).Normal().Base(#2));"
        "M$1(#1,#1)+=(t(:,i,i, :,j,j, k).gamma(k));");

    // Assign the M_mediumMesh integration method
    assem_surf.push_mi(im);
    assem_surf2.push_mi(im);
    assem_cons.push_mi(im);

    // Assign the M_mediumMesh finite element space
    assem_surf.push_mf(femD);
    assem_surf2.push_mf(femD);
    assem_cons.push_mf(femD);

    // Assign the M_mediumMesh finite element space for the coefficients
    assem_surf.push_mf(femC);
    assem_surf2.push_mf(femC);
    assem_cons.push_mf(femC);

    // Assign the coefficients
    assem_surf.push_data(etaGammaUinvh);
    assem_surf2.push_data(etaGammaUinvh);
    assem_cons.push_data(lambda);
    assem_cons.push_data(mu);

    // Set the matrices to save the evaluations
    assem_surf.push_mat(*M);
    assem_surf2.push_mat(*M);
    assem_cons.push_mat(*M);

    // Assemble in each sub region
    for ( size_type bndID = 0; bndID < bcPtr->getDiriBD().size(); bndID++ )
    {
             assem_surf.assembly( medium->getMesh()->region(bcPtr->getDiriBD()[bndID]));
	     assem_surf2.assembly( medium->getMesh()->region(bcPtr->getDiriBD()[bndID]));
	     assem_cons.assembly( medium->getMesh()->region(bcPtr->getDiriBD()[bndID]));
    }

    for ( size_type bndID = 0; bndID < bcPtr->getMixedBD().size(); bndID++ )
    {
             assem_surf2.assembly( medium->getMesh()->region(bcPtr->getMixedBD()[bndID]));
    }
    
 
 //   std::cout << "elast :: operator nitsche    [OK]" << std::endl;


}


void essentialWNitscheRHSVec( scalarVectorPtr_Type V,
               Bulk* medium, scalar_type time, BC* bcPtr,  FEM& FemD, FEM& FemC, getfem::mesh_im& im)
               
{
    getfem::mesh_fem femD(*(FemD.getFEM()));
    getfem::mesh_fem femC((*FemC.getFEM()));
    
    size_type dim=medium->getDim();
   
    scalarVector_Type etaGammaUinvh(femC.nb_dof());
    std::vector<scalar_type> lambda(femC.nb_dof(),0.0);
    std::vector<scalar_type> mu(femC.nb_dof(),0.0);

    for ( size_type i = 0; i < femC.nb_dof(); ++i )
    {
        scalar_type h_el= medium->getMesh()->convex_radius_estimate(femC.first_convex_of_basic_dof(i));
        scalar_type inv_h_mesh = 1.0 / h_el;
	    lambda[i]=medium->getElastData()->getLambda(i);
    	mu[i]=medium->getElastData()->getMu(i);
        etaGammaUinvh [ i ] = 5*mu[i] *inv_h_mesh;
    }

    for ( size_type bndID = 0; bndID < bcPtr->getDiriBD().size(); bndID++ )
    {
    getfem::generic_assembly assem_surf, assem_surf2, assem_cons;


        if (dim == 2) {
            assem_surf.set("datax=data$1(#2);""datay=data$2(#2);"
                "t=comp(vBase(#1).Base(#2));"
                "V(#1)+=t(:,1,k).datax(k)+t(:,2,k).datay(k);");

            assem_surf2.set("datax=data$1(#2);""datay=data$2(#2);"
                "t=comp(vBase(#1).Normal().Normal().Base(#2));"
                "V(#1)+=t(:,1,1,1,k).datax(k)+t(:,2,2,2,k).datay(k);");

            assem_cons.set("lambda=data$1(#3); mu=data$2(#3);" "datax=data$3(#2);""datay=data$4(#2);"
                           "t=-comp(Base(#2).vGrad(#1).Normal().Base(#3));"
             "V(#1)+=t(l,:,1,j,j,k).mu(k).datax(l)+t(l,:,i,1,i,k).mu(k).datax(l)+t(l,:,i,i,1,k).lambda(k).datax(l)+t(l,:,2,j,j,k).mu(k).datay(l)+t(l,:,i,2,i,k).mu(k).datay(l)+t(l,:,i,i,2,k).lambda(k).datay(l)");
        }
        else if (dim == 3) {
            assem_surf.set("datax=data$1(#2);""datay=data$2(#2);""dataz=data$3(#2);"
                "t=comp(vBase(#1).Base(#2));"
                "V(#1)+=t(:,1,k).datax(k)+t(:,2,k).datay(k)+t(:,3,k).dataz(k);");

            assem_surf2.set("datax=data$1(#2);""datay=data$2(#2);""dataz=data$3(#2);"
                "t=comp(vBase(#1).Normal().Normal().Base(#2));"
                "V(#1)+=t(:,1,1,1,k).datax(k)+t(:,2,2,2,k).datay(k)+t(:,3,3,3,k).dataz(k);");

            assem_cons.set("lambda=data$1(#3); mu=data$2(#3);" "datax=data$3(#2);""datay=data$4(#2);""dataz=data$5(#2);"
                           "t=-comp(Base(#2).vGrad(#1).Normal().Base(#3));"
             "V(#1)+=t(l,:,1,j,j,k).mu(k).datax(l)+t(l,:,i,1,i,k).mu(k).datax(l)+t(l,:,i,i,1,k).lambda(k).datax(l)"
                   "+t(l,:,2,j,j,k).mu(k).datay(l)+t(l,:,i,2,i,k).mu(k).datay(l)+t(l,:,i,i,2,k).lambda(k).datay(l)"
                   "+t(l,:,3,j,j,k).mu(k).dataz(l)+t(l,:,i,3,i,k).mu(k).dataz(l)+t(l,:,i,i,3,k).lambda(k).dataz(l)");
        }
        else {
            throw std::runtime_error("essentialWNitscheRHSVec: unsupported dimension " + std::to_string(dim));
        }

    // Assign the M_mediumMesh integration method
    assem_surf.push_mi(im);
    assem_surf2.push_mi(im);
    assem_cons.push_mi(im);

    // Assign the M_mediumMesh finite element space
    assem_surf.push_mf(femD);
    assem_surf2.push_mf(femD);
    assem_cons.push_mf(femD);
    assem_cons.push_mf(femC); 
    assem_cons.push_mf(femC);
   

    // Assign the M_mediumMesh finite element space for the coefficients
    assem_surf.push_mf(femC);
    assem_surf2.push_mf(femC);

    // Assemble in each sub region
    
     	    assem_cons.push_data(lambda);
  	    assem_cons.push_data(mu);
	    scalarVector_Type datax(femC.nb_dof());
	    scalarVector_Type datax2(femC.nb_dof());
	    for (size_type i=0; i<femC.nb_dof();++i)
	    {
	    	datax [ i ] = (bcPtr->BCDiriVec(femC.point_of_basic_dof(i), bcPtr->getDiriBD()[bndID], time)[0]+
				bcPtr->BCDiriVel(femC.point_of_basic_dof(i), bcPtr->getDiriBD()[bndID], time)[0]*time);
	    }
    
	    for (size_type i=0; i<femC.nb_dof();++i)
	    {
	    	datax2 [ i ] = etaGammaUinvh [ i ]*datax[i];
	    }

	
   	    
   	scalarVector_Type datay(femC.nb_dof());
    scalarVector_Type datay2(femC.nb_dof());
    
	    for (size_type i=0; i<femC.nb_dof();++i)
	    {
	    	datay [ i ] = (bcPtr->BCDiriVec(femC.point_of_basic_dof(i), bcPtr->getDiriBD()[bndID],time)[1] + 
				bcPtr->BCDiriVel(femC.point_of_basic_dof(i), bcPtr->getDiriBD()[bndID],time)[1]*time);
	    }

        for (size_type i=0; i<femC.nb_dof();++i)
	    {
	    	datay2 [ i ] = etaGammaUinvh [ i ]*datay[i];
	    }
    	    // Assign the coefficients
    assem_cons.push_data(datax);
    
   	    assem_surf.push_data(datax2);
    assem_surf2.push_data(datax2);
    assem_cons.push_data(datay);
   	    assem_surf.push_data(datay2);
    assem_surf2.push_data(datay2);

   	     if (dim == 3) {
            scalarVector_Type dataz(femC.nb_dof());
            scalarVector_Type dataz2(femC.nb_dof());
            
            for (size_type i=0; i<femC.nb_dof();++i)
            {
                dataz [ i ] = (bcPtr->BCDiriVec(femC.point_of_basic_dof(i), bcPtr->getDiriBD()[bndID],time)[2] + 
                        bcPtr->BCDiriVel(femC.point_of_basic_dof(i), bcPtr->getDiriBD()[bndID],time)[2]*time);
            }
            for (size_type i=0; i<femC.nb_dof();++i)
            {
                dataz2 [ i ] = etaGammaUinvh [ i ]*dataz[i];
            }
            
            assem_cons.push_data(dataz);
            assem_surf.push_data(dataz2);
            assem_surf2.push_data(dataz2);
        }
	    // Set the matrices to save the evaluations
   	    assem_surf.push_vec(*V);
    assem_surf2.push_vec(*V);
	    assem_cons.push_vec(*V);

            assem_surf.assembly(medium->getMesh()->region(bcPtr->getDiriBD()[bndID]));
    assem_surf2.assembly(medium->getMesh()->region(bcPtr->getDiriBD()[bndID]));
	    assem_cons.assembly(medium->getMesh()->region(bcPtr->getDiriBD()[bndID]));

    }
//    std::cout << "elast :: operator nitsche    [OK]" << std::endl;

}


void stressRHS( scalarVectorPtr_Type V,
               Bulk* medium, scalar_type time, BC* bcPtr,  FEM& FemD, FEM& FemC, getfem::mesh_im& im)
               
{
    getfem::mesh_fem femD(*(FemD.getFEM()));
    getfem::mesh_fem femC((*FemC.getFEM()));
    
    size_type dim = medium->getDim();

    // Assemble in each sub region
    for ( size_type bndID = 0; bndID < bcPtr->getNeumBD().size(); bndID++ )
    {
        getfem::generic_assembly assem_surf;

        if (dim == 2) {
            assem_surf.set("datax=data$1(#2);""datay=data$2(#2);"
                "t=comp(vBase(#1).Base(#2));"
                "V(#1)+=t(:,1,k).datax(k)+t(:,2,k).datay(k);");
        }
        else if (dim == 3) {
            assem_surf.set("datax=data$1(#2);""datay=data$2(#2);""dataz=data$3(#2);"
                "t=comp(vBase(#1).Base(#2));"
                "V(#1)+=t(:,1,k).datax(k)+t(:,2,k).datay(k)+t(:,3,k).dataz(k);");
        }
        else {
            throw std::runtime_error("stressRHS: unsupported dimension " + std::to_string(dim));
        }

        // Assign the mesh integration method
        assem_surf.push_mi(im);

        // Assign the mesh finite element space
        assem_surf.push_mf(femD);

        // Assign the mesh finite element space for the coefficients
        assem_surf.push_mf(femC);
  
        // Declare ALL data vectors at same scope to ensure they remain valid until assembly
        scalarVector_Type datax(femC.nb_dof());
        scalarVector_Type datay(femC.nb_dof());
        scalarVector_Type dataz;  // Only populated for 3D
    
        for (size_type i=0; i<femC.nb_dof();++i)
        {
            datax [ i ] = bcPtr->BCNeumVec(femC.point_of_basic_dof(i), bcPtr->getNeumBD()[bndID], time)[0];
        }
    
        for (size_type i=0; i<femC.nb_dof();++i)
        {
            datay [ i ] = bcPtr->BCNeumVec(femC.point_of_basic_dof(i), bcPtr->getNeumBD()[bndID], time)[1];
        }
     
        // Assign the coefficients
        assem_surf.push_data(datax);
        assem_surf.push_data(datay);
        
        if (dim == 3) {
            dataz.resize(femC.nb_dof());
        
            for (size_type i=0; i<femC.nb_dof();++i)
            {
                dataz [ i ] = bcPtr->BCNeumVec(femC.point_of_basic_dof(i), bcPtr->getNeumBD()[bndID], time)[2];
            }
            assem_surf.push_data(dataz);
        }
        
        // Set the matrices to save the evaluations
        assem_surf.push_vec(*V);
       
        assem_surf.assembly(medium->getMesh()->region(bcPtr->getNeumBD()[bndID]));

    }
}