#include "../include/FEM.h"

FEM::FEM ( const getfem::mesh* mesh,
	   const GetPot& dataFile,
	   const std::string& problem,
	   const std::string& variable,
	   const std::string& section) :
           M_section ( section + problem ),
           M_femType ( ),
           M_SpaceDim( ),
           M_FEM(*mesh),
           M_meshPtr(mesh)
{

 
   if (section=="fractureData/")
   {
  	 M_femType= dataFile ( ( M_section+ "FEMType"+variable ).data (), "FEM_PK(1,1)" ) ;
         M_SpaceDim= dataFile ( ( M_section+ "spaceDimension" ).data (), 1 );
   }
   else
   {
   	 M_femType = dataFile ( ( M_section+ "FEMType"+variable ).data (), "FEM_PK(2,1)" );
         M_SpaceDim= dataFile ( ( M_section+ "spaceDimension" ).data (), 2 );
   }
   getfem::pfem pf_v;
   pf_v = getfem::fem_descriptor(M_femType);
   if (variable=="Velocity"||variable=="VelocityVis"|| variable=="Displacement" || variable=="Stress")
   {
   	M_FEM.set_qdim(M_SpaceDim);

   }
 
   M_FEM.set_finite_element(mesh->convex_index(), pf_v);

   for (size_type i=0; i<M_FEM.nb_dof();++i)
   {
   	M_DOFpoints.push_back(M_FEM.point_of_basic_dof(i));
   }

   for (size_type i=0; i<M_extended.size();++i)
   {
   	M_DOFpoints.push_back(M_FEM.point_of_basic_dof(M_extended[i]));
   }
}
// questo secondo costruttore non si basa sul file di input

FEM::FEM (const getfem::mesh* mesh,
	         std::string femType,
   		 size_type spaceDim
 		
              ):
		M_femType ( ),
          	 M_SpaceDim( ),
         	  M_FEM(*mesh),
         	  M_meshPtr(mesh)
{
   getfem::pfem pf_v;
   pf_v = getfem::fem_descriptor(femType);
   M_FEM.set_qdim(spaceDim);
   M_SpaceDim=spaceDim;
   M_femType=femType;
   
 
   M_FEM.set_finite_element(mesh->convex_index(), pf_v);

   for (size_type i=0; i<M_FEM.nb_dof();++i)
   {
   	M_DOFpoints.push_back(M_FEM.point_of_basic_dof(i));
   }

   for (size_type i=0; i<M_extended.size();++i)
   {
   	M_DOFpoints.push_back(M_FEM.point_of_basic_dof(M_extended[i]));
   }
}



size_type FEM::nb_dof(std::string which)
{
	if (which=="base")
	{
		return M_FEM.nb_dof();
	}
	else
	{
		if (which=="extended")
		{
			return M_extended.size();
		}
		else
		{
			return M_FEM.nb_dof()+M_extended.size();
		}
	}
}     
         
