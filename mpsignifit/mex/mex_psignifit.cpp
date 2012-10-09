#include "mex.h"
#include "psipp.h"
#include <string.h>

PsiData * make_dataset ( const mxArray * );

void logmsg ( const char *msg ) {
	mexPrintf ( msg );
	mexEvalString ( "drawnow" );
}

unsigned int get_nafc ( const mxArray * inp ) {
	mxArray *nafc = mxGetField ( inp, 0, "nafc" );
	if ( mxIsDouble ( nafc ) )
		return *( (double*)mxGetData ( nafc ) );
	else
		mexErrMsgTxt ( "(mex) get_nafc: could not interprete nafc field\n" );
}

unsigned int get_nsamples ( const mxArray * inp ) {
	if ( mxGetFieldNumber ( inp, "nsamples" ) <0 )
		mexErrMsgTxt ( "(mex) get_nsamples: 'nsamples' was not specified!\n" );
	mxArray *nsamples = mxGetField ( inp, 0, "nsamples" );
	if ( mxIsDouble ( nsamples ) )
		return *( (double*)mxGetData ( nsamples ) );
	else
		mexErrMsgTxt ( "(mex) get_nsamples: could not interprete nsamples field\n" );
}

bool get_parametric ( const mxArray * inp ) {
	if ( mxGetFieldNumber ( inp, "parametric" ) <0 )
		mexErrMsgTxt ( "(mex) get_parametric: 'parametric' was not specified!\n" );
	mxArray * parametric = mxGetField ( inp, 0, "parametric" );
	if ( mxIsDouble ( parametric ) )
		return *( (double*)mxGetData ( parametric ) );
	else
		mexErrMsgTxt ( "(mex) get_parametric: could not interprete parametric field\n" );
}

bool get_gammaislambda ( const mxArray * inp ) {
	if ( mxGetFieldNumber ( inp, "gammaislambda" ) <0 )
		mexErrMsgTxt ( "(mex) get_gammaislambda: 'gammaislambda' was not specified!\n" );
	mxArray * gammaislambda = mxGetField ( inp, 0, "gammaislambda" );
	if ( mxIsLogical ( gammaislambda ) )
		return *( (double*)mxGetData ( gammaislambda ) );
	else
		mexErrMsgTxt ( "(mex) get_gammaislambda: could not interprete gammaislambda field\n" );
}

PsiSigmoid * get_sigmoid ( const mxArray * inp ) {
	if ( mxGetFieldNumber ( inp, "sigmoid" )  <0 )
		mexErrMsgTxt ( "(mex) get_sigmoid: 'sigmoid' was not specified!\n" );

	char ssigmoid[10];
	mxGetString ( mxGetField ( inp, 0, "sigmoid" ), ssigmoid, 10 );

	if ( !strcmp ( "logistic", ssigmoid ) )
		return new PsiLogistic;
	else if ( !strcmp ( "cauchy", ssigmoid ) )
		return new PsiCauchy;
	else if ( !strcmp ( "exp", ssigmoid ) )
		return new PsiExponential;
	else if ( !strcmp ( "gauss", ssigmoid ) )
		return new PsiGauss;
	else if ( !strcmp ( "gumbel_l", ssigmoid ) )
		return new PsiGumbelL;
	else if ( !strcmp ( "gumbel_r", ssigmoid ) )
		return new PsiGumbelR;
	else
		mexErrMsgTxt ( "(mex) get_sigmoid: unknown sigmoid\n" );
}

void get_core_and_sigmoid ( const mxArray * inp, PsiCore **core, PsiSigmoid **sigmoid ) {
	if ( mxGetFieldNumber ( inp, "core" ) <0 )
		mexErrMsgTxt ( "(mex) get_core_and_sigmoid: 'core' was not specified!\n" );

	char score[10];
	mxGetString ( mxGetField ( inp, 0, "core" ), score, 10 );
	PsiData *pdata;
	double alpha;

	*sigmoid = get_sigmoid ( inp );
	// core = NULL;

	if ( !strcmp ( "ab", score ) ) {
		logmsg ( "Using ab core\n" );
		*core = new abCore;
	} else if ( !strcmp ( "linear", score ) ) {
		logmsg ( "Using linear core\n" );
		*core = new linearCore ;
	} else if ( !strncmp ( "mw", score, 2 ) ) {
		sscanf ( score, "mw%lf", &alpha );
		mexPrintf ( "Using mw%g core (extracted from %s)\n", alpha, score );
		*core = new mwCore ( NULL, (*sigmoid)->getcode (), alpha );
	} else {
		pdata = make_dataset ( inp );
		if ( !strcmp ( "log", score ) ) {
			logmsg ( "Using log core\n" );
			*core = new logCore ( pdata );
		} else if ( !strcmp ( "weibull", score ) ) {
			logmsg ( "Using weibull core\n" );
			*core = new weibullCore ( pdata );
		} else if ( !strcmp ( "poly", score ) ) {
			logmsg ( "Using poly core\n" );
			*core = new polyCore ( pdata );
		}
		delete pdata;
	}

	mexPrintf ( "sigmoid.code = %d\n", (*sigmoid)->getcode() );
	// logmsg ( "core: %s", core->getDescriptor ().c_str() );

	if ( *core==NULL )
		mexErrMsgTxt ( "(mex) get_core_and_sigmoid: unknown core\n" );
	else
		logmsg ( "core and sigmoid ok\n" );
}

std::vector<double> get_cuts ( const mxArray * inp ) {
	if ( mxGetFieldNumber ( inp, "cuts" ) <0 )
		mexErrMsgTxt ( "(mex) get_cuts: 'cuts' was not specified!\n" );

	mxArray *cuts = mxGetField ( inp, 0, "cuts" );
	unsigned int i;
	unsigned int M ( mxGetM ( cuts ) );
	unsigned int N ( mxGetN ( cuts ) );
	if ( ! (M==1||N==1) )
		mexErrMsgTxt ( "(mex) get_cuts: 'cuts' should be one dimensional!\n" );

	std::vector<double> ccuts ( N*M );
	for ( i=0; i<N*M; i++ ) {
		ccuts[i] = ((double*)mxGetData ( cuts ))[i];
	}

	return ccuts;
}

std::vector<double> get_params ( const mxArray * inp ) {
	if ( mxGetFieldNumber ( inp, "params" ) <0 )
		mexErrMsgTxt ( "(mex) get_params: 'params' was not specified!\n" );

	mxArray *params = mxGetField ( inp, 0, "params" );
	unsigned int i;
	unsigned int M ( mxGetM ( params ) );
	unsigned int N ( mxGetN ( params ) );
	if ( !(M==1||N==1) )
		mexErrMsgTxt ( "(mex) get_params: 'params' should be one dimensional!\n" );

	std::vector<double> cparams ( N*M );
	for ( i=0; i<N*M; i++ ) {
		cparams[i] = ((double*)mxGetData ( params ) )[i];
	}

	return cparams;
}

PsiData *make_dataset ( const mxArray * inp ) {
	// Extract data set from input

	// Check the input
	if ( mxGetFieldNumber ( inp, "nafc" ) <0 )
		mexErrMsgTxt ( "(mex) make_dataset: 'nafc' was not specified!\n" );
	if ( mxGetFieldNumber ( inp, "data" ) <0 )
		mexErrMsgTxt ( "(mex) make_dataset: 'data' was not specified!\n" );
	mxArray *data;
	data = mxGetField ( inp, 0, "data" );
	if ( mxGetN ( data ) != 3 )
		mexErrMsgTxt ( "Data should have three columns!" );
	// Put data to a cpointer
	double *cdata = (double*)mxGetData ( data );

	unsigned int nblocks = mxGetM ( data );
	mexPrintf ( "Allocating %d data blocks\n", nblocks );

	std::vector<double> x ( nblocks );
	std::vector<int>    N ( nblocks );
	std::vector<int>    k ( nblocks );

	int nafc;
	unsigned int i;

	// Extract relevant information
	nafc = *((double*)mxGetData ( mxGetField ( inp, 0, "nafc" ) ));

	mexPrintf ( "Dataset from %d-AFC task\n", nafc );
	for ( i=0; i<nblocks; i++ ) {
		x[i] = cdata[i];
		k[i] = cdata[nblocks+i];
		N[i] = cdata[2*nblocks+i];
		mexPrintf ( "%6.3f %3i %3i\n", x[i], k[i], N[i] );
	}

	return new PsiData ( x, N, k, nafc );
}

PsiPrior *get_prior ( const mxArray * inp, unsigned int prm ) {
	if ( mxGetFieldNumber ( inp, "priors" ) <0 )
		mexErrMsgTxt ( "(mex) get_prior: 'priors' was not specified!\n" );

	char sprior[20];
	const mxArray * priors = mxGetField ( inp, 0, "priors" );
	if ( mxGetNumberOfElements ( priors ) < prm+1 )
		mexErrMsgTxt ( "(mex) get_prior: not enough priors specified\n" );
	mxGetString ( mxGetCell ( priors, prm ), sprior, 20 );
	double al,bt;

	if ( !strncmp ( "Gauss", sprior, 5 ) ) {
		sscanf ( sprior, "Gauss(%lf,%lf)", &al, &bt );
		return new GaussPrior ( al, bt );
	} else if ( !strncmp ( "Gamma", sprior, 5 ) ) {
		sscanf ( sprior, "Gamma(%lf,%lf)", &al, &bt );
		return new GammaPrior ( al, bt );
	} else if ( !strncmp ( "nGamma", sprior, 6 ) ) {
		sscanf ( sprior, "nGamma(%lf,%lf)", &al, &bt );
		return new nGammaPrior ( al, bt );
	} else if ( !strncmp ( "Beta", sprior, 4 ) ) {
		sscanf ( sprior, "Beta(%lf,%lf)", &al, &bt );
		return new BetaPrior ( al, bt );
	} else if ( !strncmp ( "Uniform", sprior, 4 ) ) {
		sscanf ( sprior, "UniformPrior(%lf,%lf)", &al, &bt );
		return new UniformPrior ( al, bt );
	} else if ( sprior, !strncmp ( "unconstrained", sprior, 5 ) || !strncmp ( "flat", sprior, 5 ) ) {
		return NULL;
	} else {
		mexErrMsgTxt ( "(mex) get_prior: unknown prior" );
	}
}

PsiPsychometric *make_pmf ( const mxArray * inp ) {
	// Extract a psychometric function model from input

	// Check input first
	if ( mxGetFieldNumber ( inp, "nafc" ) <0 )
		mexErrMsgTxt ( "(mex) make_pmf: 'nafc' was not specified!\n"  );
	if ( mxGetFieldNumber ( inp, "sigmoid" ) <0 )
		mexErrMsgTxt ( "(mex) make_pmf: 'sigmoid' was not specified!\n" );

	unsigned int nafc ( get_nafc ( inp ) ), i;
	PsiSigmoid *sigmoid = NULL;
	PsiCore *core = NULL;
	PsiPrior *prior;
	get_core_and_sigmoid ( inp, &core, &sigmoid );
	mexPrintf ( "HL: sigmoid.code = %d", sigmoid->getcode () );

	PsiPsychometric * pmf = new PsiPsychometric ( nafc, core, sigmoid );
	if ( get_nafc ( inp )==1 && get_gammaislambda ( inp ) )
		pmf->setgammatolambda();
	delete core;
	delete sigmoid;

	for ( i=0; i<pmf->getNparams(); i++ ) {
		prior = get_prior ( inp, i );
		if ( prior != NULL ) {
			pmf->setPrior ( i, prior );
			delete prior;
		}
	}

	return pmf;
}

void make_dataset_and_pmf ( const mxArray * inp, PsiData **dataset, PsiPsychometric **pmf ) {
	// Extract both, dataset and psychometric function
	// in principle this is a point at which we could optimize a bit
	*dataset = make_dataset ( inp );
	*pmf     = make_pmf ( inp );
}

/***************************** Interface methods **********************************/

mxArray *bootstrap ( const mxArray *inp ) {
	PsiData *psidata;
	PsiPsychometric *pmf;
	mxArray *output;
	unsigned int i,j;
	const char *fieldnames[25];
	std::vector<double> cuts ( get_cuts ( inp ) );
	std::vector<int> resp;
	unsigned int ncuts ( cuts.size() );

	make_dataset_and_pmf ( inp, &psidata, &pmf );

	unsigned int nblocks ( psidata->getNblocks() );

	for ( i=0; i<6; i++ ) {
		mexPrintf ( "%6.3f %d %d\n", psidata->getIntensity ( i ), psidata->getNcorrect ( i ), psidata->getNtrials ( i ) );
	}

	logmsg ( "\nRunning simulations ... BS ... " );
	unsigned int nsamples = get_nsamples ( inp );
	BootstrapList bs_list = bootstrap ( nsamples, psidata, pmf, cuts, NULL, true, get_parametric ( inp ) );
	logmsg ( "JK ..." );
	JackKnifeList jk_list = jackknifedata ( psidata, pmf );

	logmsg ( "converting output data" );
	fieldnames[0]  = (char*)"mcsamples";
	fieldnames[1]  = (char*)"mcestimates";
	fieldnames[2]  = (char*)"mcdeviance";
	fieldnames[3]  = (char*)"mcthres";
	fieldnames[4]  = (char*)"thbias";
	fieldnames[5]  = (char*)"thacc";
	fieldnames[6]  = (char*)"mcslope";
	fieldnames[7]  = (char*)"slbias";
	fieldnames[8]  = (char*)"slacc";
	fieldnames[9]  = (char*)"mcRpd";
	fieldnames[10] = (char*)"mcRkd";
	fieldnames[11] = (char*)"outliers";
	fieldnames[12] = (char*)"influential";
	output = mxCreateStructMatrix ( 1, 1, 13, fieldnames );

	mxArray *mcsamples   = mxCreateNumericMatrix ( nsamples, psidata->getNblocks(), mxDOUBLE_CLASS, mxREAL );
	mxArray *mcestimates = mxCreateNumericMatrix ( nsamples, pmf->getNparams(), mxDOUBLE_CLASS, mxREAL );
	mxArray *mcdeviance  = mxCreateNumericMatrix ( nsamples, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *mcthres     = mxCreateNumericMatrix ( nsamples, ncuts, mxDOUBLE_CLASS, mxREAL );
	mxArray *mcslope     = mxCreateNumericMatrix ( nsamples, ncuts, mxDOUBLE_CLASS, mxREAL );
	mxArray *mcRpd       = mxCreateNumericMatrix ( nsamples, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *mcRkd       = mxCreateNumericMatrix ( nsamples, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *outliers    = mxCreateNumericMatrix ( nblocks, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *influential = mxCreateNumericMatrix ( nblocks, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *thbias      = mxCreateNumericMatrix ( ncuts, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *slbias      = mxCreateNumericMatrix ( ncuts, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *thacc       = mxCreateNumericMatrix ( ncuts, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *slacc       = mxCreateNumericMatrix ( ncuts, 1, mxDOUBLE_CLASS, mxREAL );
	double *cmcsamples   = mxGetPr ( mcsamples );
	double *cmcestimates = mxGetPr ( mcestimates );
	double *cmcdeviance  = mxGetPr ( mcdeviance );
	double *cmcthres     = mxGetPr ( mcthres );
	double *cmcslope     = mxGetPr ( mcslope );
	double *cmcRpd       = mxGetPr ( mcRpd );
	double *cmcRkd       = mxGetPr ( mcRkd );
	double *coutliers    = mxGetPr ( outliers );
	double *cinfluential = mxGetPr ( influential );
	double *cthbias      = mxGetPr ( thbias );
	double *cslbias      = mxGetPr ( slbias );
	double *cthacc       = mxGetPr ( thacc );
	double *cslacc       = mxGetPr ( slacc );

	std::vector<double> ci_lower (pmf->getNparams () );
	std::vector<double> ci_upper (pmf->getNparams () );

	for ( i=0; i<nsamples; i++ ) {
		resp = bs_list.getData ( i );
		for ( j=0; j<nblocks; j++ )
			cmcsamples[j*nsamples + i] = resp[j];
		for ( j=0; j<pmf->getNparams(); j++ )
			cmcestimates[j*nsamples +i] = bs_list.getEst ( i, j );
		cmcdeviance[i] = bs_list.getdeviance ( i );
		for ( j=0; j<ncuts; j++ ) {
			cmcthres[j*nsamples + i] = bs_list.getThres_byPos ( i, j );
			cmcslope[j*nsamples + i] = bs_list.getSlope_byPos ( i, j );
		}
		cmcRpd[i] = bs_list.getRpd ( i );
		cmcRkd[i] = bs_list.getRkd ( i );
	}

	for ( i=0; i<pmf->getNparams(); i++ ) {
		ci_lower[i] = bs_list.getPercentile ( 0.025, i );
		ci_upper[i] = bs_list.getPercentile ( 0.975, i );
	}
	for ( i=0; i<nblocks; i++ ) {
		coutliers[i]    = jk_list.outlier ( i );
		cinfluential[i] = jk_list.influential ( i, ci_lower, ci_upper );
	}

	for ( i=0; i<ncuts; i++ ) {
		cthacc[i]  = bs_list.getAcc_t ( i );
		cthbias[i] = bs_list.getBias_t ( i );
		cslacc[i] = bs_list.getAcc_s ( i );
		cthacc[i] = bs_list.getAcc_t ( i );
	}

	mxSetField ( output, 0, "mcsamples",   mcsamples );
	mxSetField ( output, 0, "mcestimates", mcestimates );
	mxSetField ( output, 0, "mcdeviance",  mcdeviance );
	mxSetField ( output, 0, "mcthres",     mcthres );
	mxSetField ( output, 0, "mcslope",     mcslope );
	mxSetField ( output, 0, "mcRpd",       mcRpd );
	mxSetField ( output, 0, "mcRkd",       mcRkd );
	mxSetField ( output, 0, "thbias",      thbias );
	mxSetField ( output, 0, "thacc",       thacc );
	mxSetField ( output, 0, "slbias",      slbias );
	mxSetField ( output, 0, "slacc",       slacc );
	mxSetField ( output, 0, "outliers",    outliers );
	mxSetField ( output, 0, "influential", influential );

	return output;
}

mxArray *mapestimate ( const mxArray *inp ) {
	PsiData *psidata;
	PsiPsychometric *pmf;
	mxArray *output;
	unsigned int i,j;
	const char *fieldnames[25];
	std::vector<double> cuts ( get_cuts ( inp ) );
	unsigned int ncuts ( cuts.size() );

	make_dataset_and_pmf ( inp, &psidata, &pmf );

	unsigned int nblocks ( psidata->getNblocks() ), nparams ( pmf->getNparams() );

	for ( i=0; i<6; i++ ) {
		mexPrintf ( "%6.3f %d %d\n", psidata->getIntensity ( i ), psidata->getNcorrect ( i ), psidata->getNtrials ( i ) );
	}

	PsiOptimizer opt ( pmf, psidata );
	std::vector<double> estimate ( opt.optimize ( pmf, psidata ) ); // NOTE: No starting values...

	Matrix *H = pmf->ddnegllikeli ( estimate, psidata );

	mxArray *thres    = mxCreateNumericMatrix ( ncuts, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *slope    = mxCreateNumericMatrix ( ncuts, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *deviance = mxCreateNumericMatrix ( 1, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *mapest   = mxCreateNumericMatrix ( nparams, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *fisher   = mxCreateNumericMatrix ( nparams, nparams, mxDOUBLE_CLASS, mxREAL );

	fieldnames[0] = "estimate";
	fieldnames[1] = "fisher";
	fieldnames[2] = "thres";
	fieldnames[3] = "slope";
	fieldnames[4] = "deviance";
	output = mxCreateStructMatrix ( 1, 1, 5, fieldnames );

	double *cthres    = mxGetPr ( thres );
	double *cslope    = mxGetPr ( slope );
	double *cestimate = mxGetPr ( mapest );
	double *cfisher   = mxGetPr ( fisher );
	double *cdeviance = mxGetPr ( deviance );

	for ( i=0; i<ncuts; i++ ) {
		cthres[i] = pmf->getThres ( estimate, cuts[i] );
		cslope[i] = pmf->getSlope ( estimate, cthres[i] );
	}
	*cdeviance = pmf->deviance ( estimate, psidata );
	for ( i=0; i<nparams; i++ ) {
		cestimate[i] = estimate[i];
		for ( j=0; j<nparams; j++ ) {
			cfisher[j+nparams*i] = (*H)(i,j);
		}
	}

	mxSetField ( output, 0, "estimate", mapest );
	mxSetField ( output, 0, "fisher", fisher );
	mxSetField ( output, 0, "thres", thres );
	mxSetField ( output, 0, "slope", slope );
	mxSetField ( output, 0, "deviance", deviance );

	delete H;
	delete pmf;
	delete psidata;

	return output;
}

mxArray *diagnostics ( const mxArray * inp ) {
	PsiData *psidata;
	PsiPsychometric *pmf;
	mxArray *output;
	unsigned int i,j;
	const char *fieldnames[25];
	std::vector<double> cuts ( get_cuts ( inp ) );
	unsigned int ncuts ( cuts.size() );

	make_dataset_and_pmf ( inp, &psidata, &pmf );

	unsigned int nblocks ( psidata->getNblocks() ), nparams ( pmf->getNparams() );
	std::vector<double> params ( get_params ( inp ) );
	std::vector<double> deviance_residuals ( pmf->getDevianceResiduals ( params, psidata ) );

	mxArray *predicted    = mxCreateNumericMatrix ( nblocks, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *deviance_res = mxCreateNumericMatrix ( nblocks, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *deviance     = mxCreateNumericMatrix ( 1, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *thres        = mxCreateNumericMatrix ( ncuts, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *slope        = mxCreateNumericMatrix ( ncuts, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *rpd          = mxCreateNumericMatrix ( nblocks, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *rkd          = mxCreateNumericMatrix ( nblocks, 1, mxDOUBLE_CLASS, mxREAL );

	double *cpredicted     = mxGetPr ( predicted );
	double *cdeviance_res  = mxGetPr ( deviance_res );
	double *cdeviance      = mxGetPr ( deviance );
	double *cthres         = mxGetPr ( thres );
	double *cslope         = mxGetPr ( slope );
	double *crpd           = mxGetPr ( rpd );
	double *crkd           = mxGetPr ( rkd );

	for ( i=0; i<nblocks; i++ ) {
		cpredicted[i]    = pmf->evaluate ( psidata->getIntensity ( i ), params );
		cdeviance_res[i] = deviance_residuals[i];
		crpd[i]          = pmf->getRpd ( deviance_residuals, params, psidata );
		crkd[i]          = pmf->getRkd ( deviance_residuals, psidata );
	}
	*cdeviance = pmf->deviance ( params, psidata );
	for ( i=0; i<ncuts; i++ ) {
		cthres[i] = pmf->getThres ( params, cuts[i] );
		cslope[i] = pmf->getSlope ( params, cthres[i] );
	}

	fieldnames[0] = "predicted";
	fieldnames[1] = "deviance_residuals";
	fieldnames[2] = "deviance";
	fieldnames[3] = "thres";
	fieldnames[4] = "slope";
	fieldnames[5] = "rpd";
	fieldnames[6] = "rkd";
	output = mxCreateStructMatrix ( 1, 1, 7, fieldnames );
	mxSetField ( output, 0, "predicted", predicted );
	mxSetField ( output, 0, "deviance_residuals", deviance_res );
	mxSetField ( output, 0, "deviance", deviance );
	mxSetField ( output, 0, "thres", thres );
	mxSetField ( output, 0, "slope", slope );
	mxSetField ( output, 0, "rpd", rpd );
	mxSetField ( output, 0, "rkd", rkd );

	delete pmf;
	delete psidata;

	return output;
}

mxArray *asir ( const mxArray * inp ) {
	PsiData *psidata;
	PsiPsychometric *pmf;
	mxArray *output;
	unsigned int i,j;
	const char *fieldnames[25];
	std::vector<double> cuts ( get_cuts ( inp ) );
	unsigned int ncuts ( cuts.size() );
	unsigned int nsamples ( get_nsamples ( inp ) );
	std::vector<double> grid1d, marginal1d;

	make_dataset_and_pmf ( inp, &psidata, &pmf );

	unsigned int nblocks ( psidata->getNblocks() ), nparams ( pmf->getNparams() );

	PsiIndependentPosterior marginals ( independent_marginals ( pmf, psidata ) );
	MCMCList samples ( sample_posterior ( pmf, psidata, marginals, get_nsamples ( inp ), 25 ) );
	sample_diagnostics ( pmf, psidata, &samples );

	logmsg ( "Finished analysis\n" );

	mxArray *mcestimates                   = mxCreateNumericMatrix ( nsamples, nparams, mxDOUBLE_CLASS, mxREAL );
	mxArray *mcdeviance                    = mxCreateNumericMatrix ( nsamples, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *mcRpd                         = mxCreateNumericMatrix ( nsamples, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *mcRkd                         = mxCreateNumericMatrix ( nsamples, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *posterior_predictive_data     = mxCreateNumericMatrix ( nsamples, nblocks, mxDOUBLE_CLASS, mxREAL );
	mxArray *posterior_predictive_deviance = mxCreateNumericMatrix ( nsamples, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *posterior_predictive_Rpd      = mxCreateNumericMatrix ( nsamples, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *posterior_predictive_Rkd      = mxCreateNumericMatrix ( nsamples, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *logposterior_ratios           = mxCreateNumericMatrix ( nsamples, nblocks, mxDOUBLE_CLASS, mxREAL );
	mxArray *duplicates                    = mxCreateNumericMatrix ( 1, 1, mxDOUBLE_CLASS, mxREAL );
	mxArray *posterior_grids               = mxCreateNumericMatrix ( nparams, 1000, mxDOUBLE_CLASS, mxREAL );
	mxArray *posterior_margin              = mxCreateNumericMatrix ( nparams, 1000, mxDOUBLE_CLASS, mxREAL );
	mxArray *resampling_entropy            = mxCreateNumericMatrix ( 1, 1, mxDOUBLE_CLASS, mxREAL );

	logmsg ( "Finished allocating mx objects\n" );

	double *cmcestimates = mxGetPr ( mcestimates );
	double *cmcdeviance  = mxGetPr ( mcdeviance );
	double *cmcRpd       = mxGetPr ( mcRpd );
	double *cmcRkd       = mxGetPr ( mcRkd );
	double *cposterior_predictive_data     = mxGetPr ( posterior_predictive_data );
	double *cposterior_predictive_deviance = mxGetPr ( posterior_predictive_deviance );
	double *cposterior_predictive_Rpd      = mxGetPr ( posterior_predictive_Rpd );
	double *cposterior_predictive_Rkd      = mxGetPr ( posterior_predictive_Rkd );
	double *clogposterior_ratios           = mxGetPr ( logposterior_ratios );
	double *cduplicates                    = mxGetPr ( duplicates );
	double *cposterior_grids               = mxGetPr ( posterior_grids );
	double *cposterior_margin              = mxGetPr ( posterior_margin );
	double *cresampling_entropy            = mxGetPr ( resampling_entropy );
	char cposterior_approximations_m[4][100];
	char cposterior_approximations_str[4][100];

	logmsg ( "Finished getting pointers\n" );

	for ( i=0; i<nsamples; i++ ) {
		for ( j=0; j<nparams; j++ )
			cmcestimates[j*nsamples+i] = samples.getEst ( i, j );
		cmcdeviance[i] = samples.getdeviance ( i );
		cmcRpd[i]      = samples.getRpd ( i );
		cmcRkd[i]      = samples.getRkd ( i );
		for ( j=0; j<nblocks; j++ ) {
			cposterior_predictive_data[j*nsamples+i] = samples.getppData ( i, j );
			clogposterior_ratios[j*nsamples+i]       = samples.getlogratio ( i, j );
		}
		cposterior_predictive_deviance[i] = samples.getppDeviance ( i );
		cposterior_predictive_Rpd[i]      = samples.getppRpd ( i );
		cposterior_predictive_Rkd[i]      = samples.getppRkd ( i );
	}
	*cduplicates                       = samples.get_accept_rate ();
	*cresampling_entropy               = samples.get_entropy ();


	for ( j=0; j<nparams; j++ ) {
		grid1d     = marginals.get_grid ( j );
		marginal1d = marginals.get_margin ( j );
		for ( i=0; i<1000; i++ ) {
			cposterior_grids[j+i*nparams]  = grid1d[i];
			cposterior_margin[j+i*nparams] = marginal1d[i];
		}
	}

	PsiPrior *posterior;
	double al, bt;
	posterior = marginals.get_posterior ( 0 );
	al = posterior->getprm ( 0 );
	bt = posterior->getprm ( 1 );
	sprintf ( cposterior_approximations_m[0], "@(x) normpdf ( x, %lf, %lf );", al, bt );
	sprintf ( cposterior_approximations_str[0], "N(%.2f,%.2f)", al, bt );
	logmsg ( cposterior_approximations_m[0] );
	logmsg ( cposterior_approximations_str[0] );
	logmsg ( "\n" );

	posterior = marginals.get_posterior ( 1 );
	al = posterior->getprm ( 0 );
	bt = posterior->getprm ( 1 );
	sprintf ( cposterior_approximations_m[1], "@(x) gampdf ( x, %lf, %lf );", al, bt );
	sprintf ( cposterior_approximations_str[1], "Gamma(%.2f,%.2f)", al, bt );
	logmsg ( cposterior_approximations_m[1] );
	logmsg ( cposterior_approximations_str[1] );
	logmsg ( "\n" );

	posterior = marginals.get_posterior ( 2 );
	al = posterior->getprm ( 0 );
	bt = posterior->getprm ( 1 );
	sprintf ( cposterior_approximations_m[2], "@(x) betapdf ( x, %lf, %lf );", al, bt );
	sprintf ( cposterior_approximations_str[2], "Beta(%.2f,%.2f)", al, bt );
	logmsg ( cposterior_approximations_m[2] );
	logmsg ( cposterior_approximations_str[2] );
	logmsg ( "\n" );

	if ( nparams==4 ) {
		posterior = marginals.get_posterior ( 3 );
		al = posterior->getprm(0);
		bt = posterior->getprm(1);
		sprintf ( cposterior_approximations_m[3], "@(x) betapdf ( x, %lf, %lf );", al, bt );
		sprintf ( cposterior_approximations_str[3], "Beta(%.2f,%.2f)", al, bt );
		logmsg ( cposterior_approximations_m[3] );
		logmsg ( cposterior_approximations_str[3] );
		logmsg ( "\n" );
	}

	mxArray *posterior_approximations_m = mxCreateCellMatrix ( nparams, 1 );
	mxArray *posterior_approximations_str = mxCreateCellMatrix ( nparams, 1 );
	mxArray *value;
	for ( i=0; i<nparams; i++ ) {
		value = mxCreateString ( cposterior_approximations_m[i] );
		mxSetCell ( posterior_approximations_m, i, value );
		value = mxCreateString ( cposterior_approximations_str[i] );
		mxSetCell ( posterior_approximations_str, i, value );
	}


	fieldnames[0] = "mcestimates";
	fieldnames[1] = "mcdeviance";
	fieldnames[2] = "mcRpd";
	fieldnames[3] = "mcRkd";
	fieldnames[4] = "posterior_predictive_data";
	fieldnames[5] = "posterior_predictive_deviance";
	fieldnames[6] = "posterior_predictive_Rpd";
	fieldnames[7] = "posterior_predictive_Rkd";
	fieldnames[8] = "logposterior_ratios";
	fieldnames[9] = "duplicates";
	fieldnames[10] = "posterior_grids";
	fieldnames[11] = "posterior_margin";
	fieldnames[12] = "resampling_entropy";
	fieldnames[13] = "posterior_approximations_m";
	fieldnames[14] = "posterior_approximations_str";
	output = mxCreateStructMatrix ( 1, 1, 15, fieldnames );

	mxSetField ( output, 0, "mcestimates", mcestimates );
	mxSetField ( output, 0, "mcdeviance", mcdeviance );
	mxSetField ( output, 0, "mcRpd", mcRpd );
	mxSetField ( output, 0, "mcRkd", mcRkd );
	mxSetField ( output, 0, "posterior_predictive_data", posterior_predictive_data );
	mxSetField ( output, 0, "posterior_predictive_deviance", posterior_predictive_deviance );
	mxSetField ( output, 0, "posterior_predictive_Rpd", posterior_predictive_Rpd );
	mxSetField ( output, 0, "posterior_predictive_Rkd", posterior_predictive_Rkd );
	mxSetField ( output, 0, "logposterior_ratios", logposterior_ratios );
	mxSetField ( output, 0, "duplicates", duplicates );
	mxSetField ( output, 0, "posterior_grids", posterior_grids );
	mxSetField ( output, 0, "posterior_margin", posterior_margin );
	mxSetField ( output, 0, "posterior_approximations_m", posterior_approximations_m );
	mxSetField ( output, 0, "posterior_approximations_str", posterior_approximations_str );
	mxSetField ( output, 0, "resampling_entropy", resampling_entropy );
	
	return output;
}

void mexFunction ( int nhls, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
	/* Allfunctions can be called from this one in the form:
	 *
	 * mex_psignifit ( command, parameters )
	 *
	 * which will return a struct with the results.
	 * In the above call, command is the type of inference you want to do -- right now, the following keywords
	 * are supported:
	 *
	 * "bootstrap"
	 * "mcmc"   (not implementation yet, though!)
	 * "asir"
	 * "diagnostic"
	 * "map"
	 *
	 * parameters is a struct with all the parameters that you need included. An example can be found in the
	 * mkconfig.m file
	 */
	char method[20], errormsg[40];

	logmsg ( "Checking input\n" );
	if ( nrhs != 2 )
		mexErrMsgTxt ( "Incorrect number of input arguments\n" );

	if ( nhls != 1 )
		mexErrMsgTxt ( "Exactly one output argument should be given" );

	mxGetString ( prhs[0], method, 20 );

	if ( !strcmp ( "bootstrap", method ) ) {
		/*********************** bootstrap ********************************/
		logmsg ( "Performing bootstrap\n" );
		plhs[0] = bootstrap ( prhs[1] );
	} else if ( !strcmp ( "mcmc", method ) ) {
		logmsg ( "Performing mcmc\n" );
	} else if ( !strcmp ( "asir", method ) ) {
		logmsg ( "Performing asir\n" );
		plhs[0] = asir ( prhs[1] );
	} else if ( !strcmp ( "diagnostic", method ) ) {
		logmsg ( "Performing diagnostic\n" );
		plhs[0] = diagnostics ( prhs[1] );
	} else if (!strcmp ( "map", method ) ) {
		logmsg ( "Performing map estimation\n" );
		plhs[0] = mapestimate ( prhs[1] );
	} else {
		sprintf ( errormsg, "Unknown method '%s' - aborting!", method );
		mexErrMsgTxt ( errormsg );
	}
}
