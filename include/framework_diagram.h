/**
 * @page framework_diagram Poroelasticity Framework Architecture
 * 
 * @dot
 * digraph PoroelasticityFramework {
 *   rankdir=TB;
 *   node [fontname="Helvetica", fontsize=10, shape=box];
 *   edge [fontname="Helvetica", fontsize=9];
 *   
 *   // Top level: Russian Doll (extension)
 *   subgraph cluster_russian {
 *     label="Russian Doll Upscaling (Extension)";
 *     style=filled;
 *     color=lightblue;
 *     
 *     RussianDollProblem [label="RussianDollProblem", URL="\ref RussianDollProblem"];
 *     PLCProblem [label="PLCProblem", URL="\ref PLCProblem"];
 *     InterpolationManager [label="InterpolationManager", URL="\ref InterpolationManager"];
 *   }
 *   
 *   // Core poroelasticity
 *   subgraph cluster_core {
 *     label="Core Poroelasticity";
 *     style=filled;
 *     color=lightyellow;
 *     
 *     DarcyProblemT [label="DarcyProblemT", URL="\ref DarcyProblemT"];
 *     ElastProblem [label="ElastProblem", URL="\ref ElastProblem"];
 *     CoupledProblem [label="CoupledProblem", URL="\ref CoupledProblem"];
 *   }
 *   
 *   // FEM Infrastructure
 *   subgraph cluster_fem {
 *     label="FEM Infrastructure";
 *     style=filled;
 *     color=lightgreen;
 *     
 *     FEM [label="FEM", URL="\ref FEM"];
 *     TimeLoop [label="TimeLoop", URL="\ref TimeLoop"];
 *     BC [label="BC", URL="\ref BC"];
 *     LinearSystem [label="LinearSystem", URL="\ref LinearSystem"];
 *    
 *     
 *     // Operators sub-cluster
 *     subgraph cluster_operators {
 *       label="Operators";
 *       style=dashed;
 *       color=gray;
 *       
 *       DarcyOperators [label="DarcyOperators", URL="\ref DarcyOperators"];
 *       ElastOperators [label="ElastOperators", URL="\ref ElastOperators"];
 *       DarcyOperatorsBulk [label="DarcyOperatorsBulk"];
 *       ElastOperatorsBulk [label="ElastOperatorsBulk"];
 *       DarcyOperatorsBD [label="DarcyOperatorsBD"];
 *       ElastOperatorsBD [label="ElastOperatorsBD"];
 *     }
 *   }
 *   
 *   // Helper utilities
 *   subgraph cluster_helper {
 *     label="Utilities";
 *     style=filled;
 *     color=orange;
 *    
 *     GetPot [label="GetPot", URL="\ref GetPot"];
 *     Parser [label="Parser", URL="\ref Parser"];
 *     StringUtility [label="StringUtility", URL="\ref StringUtility"];
 *   }
 *   
 *   // Data structures
 *   subgraph cluster_data {
 *     label="Data Structures";
 *     style=filled;
 *     color=purple;
 *     Core [label="Core", URL="\ref Core"];
 *     Bulk [label="Bulk", URL="\ref Bulk"];
 *     BulkDarcyData [label="BulkDarcyData", URL="\ref BulkDarcyData"];
 *     BulkElastData [label="BulkElastData", URL="\ref BulkElastData"];
 *   }
 *   
 *   // Key dependencies
 *   RussianDollProblem -> CoupledProblem [label="couples"];
 *   RussianDollProblem -> PLCProblem [label="couples"];
 *   RussianDollProblem -> InterpolationManager [label="uses"];
 * 
 *   CoupledProblem -> {DarcyProblemT, ElastProblem} [label="couples",color=orange];
 *   CoupledProblem -> LinearSystem [style=dashed,color=orange];
 *   CoupledProblem -> TimeLoop [style=dashed,color=orange];
 *   CoupledProblem -> Bulk[style=dashed,color=orange];   
 * 
 *   DarcyProblemT -> FEM [style=dashed, color=blue];
 *   DarcyProblemT -> DarcyOperators[color=blue];
 *   DarcyProblemT -> BC [style=dashed, color=blue];
 *   DarcyProblemT -> TimeLoop [style=dashed, color=blue];
 *   DarcyProblemT -> LinearSystem [style=dashed, color=blue];
 *   DarcyProblemT -> Bulk [style=dashed, color=blue];
 *   DarcyOperators -> DarcyOperatorsBD [label="includes"];
 *   DarcyOperators -> DarcyOperatorsBulk [label="includes"];
 *   
 *   ElastProblem -> BC [style=dashed, color=red];
 *   ElastProblem -> LinearSystem [style=dashed, color=red];
 *   ElastProblem -> TimeLoop [style=dashed, color=red];
 *   ElastProblem -> FEM [style=dashed, color=red];
 *   ElastProblem -> Bulk [style=dashed, color=red];
 *   ElastProblem -> ElastOperators [color=red];
 *   ElastOperators -> ElastOperatorsBD [label="includes"];
 *   ElastOperators -> ElastOperatorsBulk [label="includes"];
 *   
 *   Bulk -> BulkDarcyData [label="contains"];
 *   Bulk -> BulkElastData [label="contains"];
 *   
 *   // Utilities used by everyone
 *   {DarcyProblemT, ElastProblem, FEM} -> GetPot [style=dashed];
 *   GetPot -> Parser [label="uses"];
 * }
 * @enddot
 */