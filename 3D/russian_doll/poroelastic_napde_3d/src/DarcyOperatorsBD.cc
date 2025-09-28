#include "../include/DarcyOperatorsBD.h"
/*
scalar_BD_data::scalar_BD_data ( const getfem::mesh_fem& mf_) :
    mf(mf_), U(mf_.nb_basic_dof()), N(mf_.linked_mesh().dim()), gradU(1, N)
{
    sizes_.resize(1);
    sizes_ [ 0 ] = bgeot::short_type(N);

 }

void scalar_BD_data::compute ( getfem::fem_interpolation_context& ctx,
                                      bgeot::base_tensor& t )
{
    size_type cv = ctx.convex_num();
    coeff.resize(mf.nb_basic_dof_of_element(cv));
   
    bgeot::base_node where=ctx.xreal();
    //[-pi*cos(pi*x)*sin(pi*y),-pi*cos(pi*y)*sin(pi*x)]
scalar_type pi=std::acos(-1);
    t [ 0 ] = sin(pi*where[0])*sin(pi*where[1]);  

}
*/
void essentialWNitsche( sparseMatrixPtr_Type M,
               Bulk* medium, BC* bcPtr,  FEM& FemV, FEM& FemC, getfem::mesh_im& im)
               
{
    getfem::mesh_fem femV(*(FemV.getFEM()));
    getfem::mesh_fem femC((*FemC.getFEM()));
    
    scalarVector_Type etaGammaUinvh(femC.nb_dof());
    scalar_type penaltyParameterVelocity=5*medium->Lx()*medium->Ly();

    for ( size_type i = 0; i < femC.nb_dof(); ++i )
    {
        const scalar_type firstRow = medium->getDarcyData()->getKxx(i) + medium->getDarcyData()->getKxy(i) ;
        const scalar_type secondRow = medium->getDarcyData()->getKyy(i) + medium->getDarcyData()->getKxy(i) ;
        const scalar_type norminf = 1./std::max ( firstRow, secondRow );
        scalar_type h_el= medium->getMesh()->convex_radius_estimate(femC.first_convex_of_basic_dof(i));
        scalar_type inv_h_mesh = 1.0 / h_el;
        etaGammaUinvh [ i ] = norminf * penaltyParameterVelocity *inv_h_mesh;
    }
    
    getfem::generic_assembly assem_surf;

    assem_surf.set("gamma=data$1(#2);"
        "t=comp(vBase(#1).Normal().vBase(#1).Normal().Base(#2));"
        "M$1(#1,#1)+=(t(:,i, i, :,j, j, k).gamma(k));");

    // Assign the M_mediumMesh integration method
    assem_surf.push_mi(im);

    // Assign the M_mediumMesh finite element space
    assem_surf.push_mf(femV);

    // Assign the M_mediumMesh finite element space for the coefficients
    assem_surf.push_mf(femC);

    // Assign the coefficients
    assem_surf.push_data(etaGammaUinvh);

    // Set the matrices to save the evaluations
    assem_surf.push_mat(*M);

    // Assemble in each sub region
    for ( size_type bndID = 0; bndID < bcPtr->getDiriBD().size(); bndID++ )
    {
             assem_surf.assembly(
                medium->getMesh()->region(bcPtr->getDiriBD()[bndID]));
    }
    

}

void essentialWNitscheRHS( scalarVectorPtr_Type V,
               Bulk* medium, BC* bcPtr,  FEM& FemV, FEM& FemC, getfem::mesh_im& im)
{
    getfem::mesh_fem femV(*(FemV.getFEM()));
    getfem::mesh_fem femC((*FemC.getFEM()));
    
    scalarVector_Type etaGammaUinvh(femC.nb_dof());
    scalar_type penaltyParameterVelocity=5*medium->Lx()*medium->Ly();

    for ( size_type i = 0; i < femC.nb_dof(); ++i )
    {
        const scalar_type firstRow = medium->getDarcyData()->getKxx(i) + medium->getDarcyData()->getKxy(i) ;
        const scalar_type secondRow = medium->getDarcyData()->getKyy(i) + medium->getDarcyData()->getKxy(i) ;
        const scalar_type norminf = 1./std::max ( firstRow, secondRow );
        scalar_type h_el= medium->getMesh()->convex_radius_estimate(femC.first_convex_of_basic_dof(i));
        scalar_type inv_h_mesh = 1.0 / h_el;
        etaGammaUinvh [ i ] = norminf * penaltyParameterVelocity *inv_h_mesh;
    }
    

    getfem::generic_assembly assem_surf;

    assem_surf.set("gammavel=data$1(#2);"
        "t=comp(vBase(#1).Normal().Base(#2));"
        "V$1(#1)+=(t(:,i, i, k).gammavel(k));");

    // Assign the M_mediumMesh integration method
    assem_surf.push_mi(im);

    // Assign the M_mediumMesh finite element space
    assem_surf.push_mf(femV);

    // Assign the M_mediumMesh finite element space for the coefficients
    assem_surf.push_mf(femC);

  // Assemble in each sub region
    for ( size_type bndID = 0; bndID < bcPtr->getDiriBD().size(); bndID++ )
    {
    
	    scalarVector_Type data(femC.nb_dof());
    
	    for (size_type i=0; i<femC.nb_dof();++i)
	    {
	    	data [ i ] = etaGammaUinvh [ i ]*bcPtr->BCDiri(femC.point_of_basic_dof(i), bcPtr->getDiriBD()[bndID]);
	    }
    
    	    // Assign the coefficients
   	    assem_surf.push_data(data);

	    // Set the matrices to save the evaluations
   	    assem_surf.push_vec(*V);
	   
            assem_surf.assembly(
            medium->getMesh()->region(bcPtr->getDiriBD()[bndID]));
    
    }
    
    
}

void naturalRHS( scalarVectorPtr_Type V,
               Bulk* medium, BC* bcPtr,  FEM& FemV, FEM& FemC, getfem::mesh_im& im, scalar_type time)
{

    getfem::mesh_fem femV(*(FemV.getFEM()));
    getfem::mesh_fem femC((*FemC.getFEM()));
    

    getfem::generic_assembly assem_surf;

    assem_surf.set("p=data$1(#2);"
        "t=-comp(vBase(#1).Normal().Base(#2));"
        "V$1(#1)+=(t(:,i, i, k).p(k));");

    // Assign the M_mediumMesh integration method
    assem_surf.push_mi(im);

    // Assign the M_mediumMesh finite element space
    assem_surf.push_mf(femV);

    // Assign the M_mediumMesh finite element space for the coefficients
    assem_surf.push_mf(femC);

  // Assemble in each sub region
    for ( size_type bndID = 0; bndID < bcPtr->getNeumBD().size(); bndID++ )
    {

	    scalarVector_Type data(femC.nb_dof());
    
	    for (size_type i=0; i<femC.nb_dof();++i)
	    {

	    	data [ i ] = bcPtr->BCNeum(femC.point_of_basic_dof(i), bcPtr->getNeumBD()[bndID], time);
            
	    }
    
    	    // Assign the coefficients
   	    assem_surf.push_data(data);

	    // Set the matrices to save the evaluations
   	    assem_surf.push_vec(*V);
	   
            assem_surf.assembly(
            medium->getMesh()->region(bcPtr->getNeumBD()[bndID]));
         	   
    }
    
    
}
