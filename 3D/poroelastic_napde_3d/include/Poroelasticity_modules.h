/**
 * @mainpage Poroelasticity Solver Framework
 * 
 * A comprehensive framework for poroelasticity simulation with Russian Doll upscaling.
 * 
 * @image html poroelasticity_structure.png "Framework Architecture"
 */

/**
 * @defgroup RussianDollGroup Russian Doll Upscaling
 * @brief Advanced upscaling techniques for multi-scale poroelasticity
 * 
 * This module implements the Russian Doll algorithm for efficient
 * multi-scale computation with Proper Localization Conditions (PLC).
 * 
 * @{
 */

// Forward declarations for clean documentation
/** @class RussianDollProblem */
/** @class PLCProblem */
/** @class InterpolationManager */

/** @} */ // End of RussianDollGroup

/**
 * @defgroup PoroelasticityGroup Core Poroelasticity
 * @brief Base poroelasticity problem formulation
 * 
 * Core poroelastic coupling implementation including Darcy flow
 * and linear elasticity with various boundary conditions.
 * 
 * @extends RussianDollGroup
 * @{
 */

/** @class DarcyProblemT */
/** @class ElastProblem */
/** @class CoupledProblem */

/** @} */ // End of PoroelasticityGroup

/**
 * @defgroup FEMGroup Finite Element Infrastructure
 * @ingroup PoroelasticityGroup
 * @brief FEM assembly, time integration, and linear algebra
 * 
 * Provides the computational infrastructure for discretization
 * and solution of the poroelasticity equations.
 * @{
 */

/** @class FEM */
/** @class TimeLoop */
/** @class BC */
/** @class LinearSystem */
/** @class Core */
/** @class Core_modern */

// Operators are part of FEM infrastructure
/** @class DarcyOperators */
/** @class DarcyOperatorsBD */
/** @class DarcyOperatorsBulk */
/** @class ElastOperators */
/** @class ElastOperatorsBD */
/** @class ElastOperatorsBulk */

/** @} */ // End of FEMGroup

/**
 * @defgroup HelperGroup Utility & Infrastructure
 * @ingroup PoroelasticityGroup
 * @brief Parser, input handling, and utility functions
 * 
 * Support classes for input parsing, string manipulation,
 * and general utility functions.
 * @{
 */

/** @class GetPot */
/** @class Parser */
/** @class ParserDefinitions */
/** @class ParserSpiritGrammar */
/** @class StringUtility */
/** @class UsefulFunctions */

/** @} */ // End of HelperGroup

/**
 * @defgroup DataGroup Data Structures
 * @ingroup PoroelasticityGroup
 * @brief Data containers for bulk properties
 * 
 * Data structures for storing and managing material properties
 * and simulation data.
 * @{
 */

/** @class Bulk */
/** @class BulkDarcyData */
/** @class BulkElastData */
/** @class LinePressure */

/** @} */ // End of DataGroup