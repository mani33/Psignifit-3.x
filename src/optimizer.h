/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include <vector>
#include "psychometric.h"
#include "data.h"

/** \brief Simplex optimization */
class PsiOptimizer
{
	private:
		// some variables for the internal processing of the optimization process
		int nparameters;                             // Number of parameters
		std::vector< std::vector<double> > simplex;  // data of the simplex
		std::vector<double> fx;                      // function values at the simplex nodes
		std::vector<double> centroid;                // controid of the simplex
		std::vector<double> newnode;                 // new, suggested simplex node
		std::vector<double> start;                   // starting values
		std::vector<bool>   modified;                // bookkeeping vector to indicate which simplex nodes have changed, i.e. which function values need to be updated
		void initialize_simplex (
				const PsiPsychometric * model,
				const PsiData * data,
				const std::vector<double>* startingvalue=NULL
				);
		void reevaluate_simplex ( const PsiPsychometric * model,
				const PsiData * data,
				int * maxind,
				int * minind,
				std::vector<double> * workspace
				);
		void calculate_centroid ( int excluded );
		void reflect ( int maxind );
		void expand ( int maxind );
		void shrink ( int minind );
	public:
		PsiOptimizer (
			const PsiPsychometric * model,           ///< model to be fitted (this is needed at this point only to determine the amount of internal memory that is required)
			const PsiData * data                     ///< data to be fitted (this is needed at this point only to determine the amount of internal memory that is required)
			); ///< set up everything
		~PsiOptimizer ( void );                                   ///< clean up everything
		std::vector<double> optimize (
			const PsiPsychometric * model,           ///< model to be fitted
			const PsiData * data,                    ///< data to be fitted
			const std::vector<double>* startingvalue=NULL    ///< starting value for optimization --- if this is longer the the number of parameters in the model, the additional values are used to span the simplex
			); ///< Start the optimization process
};

#endif
