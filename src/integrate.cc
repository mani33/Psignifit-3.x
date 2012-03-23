#include "integrate.h"
#include "errors.h"
#include "linalg.h"

// #define DEBUG_INTEGRATE

PsiIndependentPosterior::PsiIndependentPosterior ( unsigned int nprm,
				std::vector<PsiPrior*> posteriors,
				std::vector< std::vector<double> > x,
				std::vector< std::vector<double> > fx
				) : nparams (nprm), fitted_posteriors ( posteriors ), grids ( x ), margins ( fx ) {
	unsigned int i,j,k;
	double p;
	std::vector<double> w;
	Matrix M ( grids[0].size(), 2 );

	for ( i=0; i<nparams; i++ ) {
		for ( j=0; j<grids[i].size(); j++ ) {
			M(j,0) = margins[i][j];
			p = posteriors[i]->pdf ( grids[i][j] );
			k = 1;
			while ( std::isinf ( p ) ) {
				p = posteriors[i]->pdf ( grids[i][j+k] );
				k++;
				if ( j+k>= grids[i].size() ) p = 1e40;
			}
			if ( p!=p ) p=0;
			M(j,1) = p;
			// fitted_posteriors[i] = posteriors[i]->clone();
		}
		w = leastsq ( &M ); // Fit the determined margin using the analytical posterior (is this needed?)
#ifdef DEBUG_INTEGRATE
		std::cerr << "w = " << w[0] << "\n";
#endif
		if ( w[0]==w[0] )
			for ( j=0; j<margins[i].size(); j++ )
				margins[i][j] *= w[0];
	}
}

std::vector<double> lingrid ( double xmin, double xmax, unsigned int gridsize ) {
	unsigned int i;
	double dx;
	std::vector<double> x (gridsize);

	if ( xmin>xmax ) {
		dx = xmin;
		xmin = xmax;
		xmax = dx;
	}

	// std::cerr << "xmax = " << xmax << " xmin = " << xmin << "\n";
	dx = (xmax-xmin)/(gridsize-1);
	// std::cerr << "dx = " << dx << "\n";
	for ( i=0; i<gridsize; i++ )
		x[i] = xmin + i*dx;

	return x;
}

double max_vector ( const std::vector<double>& fx ) {
	double mx(-1e5);
	unsigned int i;

	for (i=0; i<fx.size(); i++ ) {
		if ( fx[i] > mx )
			mx = fx[i];
	}
	return mx;
}

double min_vector ( const std::vector<double>& fx ) {
	double mn(1e5);
	unsigned int i;

	for (i=0; i<fx.size(); i++ ) {
		if ( fx[i] < mn )
			mn = fx[i];
	}
	return mn;
}

// Numerical integration ///////////////////////

void normalize_probability ( const std::vector<double>& x, std::vector<double>& fx ) {
	double Z(0);
	double last_f ( fx[0] ), last_x ( x[0] );
	double f, d;
	unsigned int i;

	for ( i=1; i<x.size(); i++ ) {
		if ( fx[i]!=fx[i] || std::isinf ( fx[i] ) ) {
			continue;
		}
		// Trapez Quadradure
		f = 0.5*(fx[i]+last_f);
		d = x[i]-last_x;
		Z += f*d;
		last_f = fx[i];
		last_x = x[i];
#ifdef DEBUG_INTEGRATE
		std::cerr << "fx["<<i<< "] = " << fx[i] << "\n";
#endif
	}

	for ( i=0; i<x.size(); i++ ) {
		fx[i] /= Z;
	}

#ifdef DEBUG_INTEGRATE
	std::cerr << "Z = " << Z << "\n";
#endif
}

double numerical_mean ( const std::vector<double>& x, const std::vector<double>& fx ) {
	double m (0.);
	double last_f(fx[0]),last_x(x[0]);
	double f,d;
	unsigned int i;

#ifdef DEBUG_INTEGRATE
	std::cerr << last_f << "," << last_x << "\n";
#endif

	for ( i=1; i<x.size(); i++ ) {
		// Trapez Quadrature
		if ( fx[i] != fx[i] || std::isinf ( fx[i] ) ) {
			continue;
		}
		// f = 0.25*(x[i]+last_x)*(fx[i]+last_f);
		f = 0.5*(fx[i]*x[i] + last_f);
		d = x[i]-last_x;
		m += f*d;
		last_x = x[i];
		last_f = fx[i]*x[i];
	}

#ifdef DEBUG_INTEGRATE
	std::cerr << "E(X) = " << m << "\n";
#endif

	return m;
}

double numerical_variance ( const std::vector<double>& x, const std::vector<double>& fx, double m ) {
	double v (0.);
	double last_f(fx[0]),last_x(x[0]);
	double mx,mf,d,f;
	unsigned int i;

	for ( i=0; i<x.size(); i++ ) {
		// Trapze Quadrature
		if ( fx[i] != fx[i] || std::isinf ( fx[i] ) )
			continue;
		// mx = 0.5*(x[i]+last_x);
		// mf = 0.5*(fx[i]+last_f);
		f = (x[i]-m);
		f *= f;
		f *= fx[i];
		d  = x[i] - last_x;
		v += f*d;
		last_f = f;
		last_x = x[i];
	}

#ifdef DEBUG_INTEGRATE
	std::cerr << "V(X) = " << v << "\n";
#endif

	return v;
}
// End numerical integration /////////////////

// Moment matching ///////////////////////////
std::vector<double> match_gauss ( const std::vector<double>& x, const std::vector<double>& fx ) {
	std::vector<double> out ( 3 );

	out[0] = numerical_mean ( x, fx );
	out[1] = sqrt ( numerical_variance ( x, fx, out[0] ) );

	return out;
}

std::vector<double> match_gamma ( const std::vector<double>& x, const std::vector<double>& fx ) {
	double m ( numerical_mean ( x, fx ) );
	double v ( numerical_variance ( x, fx, m ) );
	std::vector<double> out (3);
	double k, theta;

	theta = v/m;
	k = m/theta;

	out[0] = k;
	out[1] = theta;

	return out;
}

std::vector<double> match_beta ( const std::vector<double>& x, const std::vector<double>& fx ) {
	double m ( numerical_mean ( x, fx ) );
	double v ( numerical_variance ( x, fx, m ) );
	std::vector<double> out (3);
	double al, bt;

	al = m*m*(1-m)/v - m;
	bt = m*(1-m)*(1-m)/v - (1-m);
	// al = -m + m*m/v - m*m*m/v;
	// bt = (al-m*al)/m;

	// al = eta/(eta2*eta2*eta2*v) - 1./eta2;
	// al = m*((1-m)*m/v -1 );
	// bt = al/m - al;

	out[0] = al;
	out[1] = bt;
	
#ifdef DEBUG_INTEGRATE
	std::cerr << "Matched Beta ( " << al << ", " << bt << " ), m = " << m << ", v = " << v << "\n";
#endif

	return out;
}
// End Moment matching ///////////////////////////

PsiIndependentPosterior independent_marginals (
		const PsiPsychometric *pmf,
		const PsiData *data
		)
{
	unsigned int gridsize (100);

	unsigned int nprm ( pmf->getNparams() ), i, j;
	unsigned int maxntrials ( 0 );
	double minp,minm,maxm,maxw,s,Z,p,m, count;
	std::vector< std::vector<double> > grids ( nprm );
	std::vector< std::vector<double> > margin ( nprm, std::vector<double>(gridsize) );
	std::vector< std::vector<double> > distparams (nprm, std::vector<double>(3) );
	std::vector<PsiPrior*> fitted_posteriors (nprm);

	PsiOptimizer * opt = new PsiOptimizer ( pmf, data );
	std::vector<double> MAP ( opt->optimize ( pmf, data ) );
	std::vector<double> prm ( MAP );
	delete opt;
	for ( i=0; i<nprm; i++ ) {
		// std::cerr << "================= prm " << i << "================\n";
		// Determine parameter ranges using the same routine as for starting values
		parameter_range ( data, pmf, i, &minm, &maxm );
		Z = 0;
		// This routine is slightly more narrow than we would like for our purposes
		if ( i>1 ) { minm=0; maxm=1.; }
		if ( i==1 ) { minm=0; maxm*=2; }

		grids[i] = lingrid ( minm, maxm, gridsize );

		for ( j=0; j<nprm; j++ ) { prm[j] = MAP[j]; }

		count = 0;
		for ( j=0; j<gridsize; j++ ) {
			prm[i] = grids[i][j];
			p = pmf->neglpost ( prm, data );
#ifdef DEBUG_INTEGRATE
			std::cerr << "p(" << i << "," << j << ") = " << p << " " << prm[i];
#endif
			// Clip
			if ( std::isinf ( p ) || p!=p ) {
				if ( j>0 ) {
					grids[i][j] = grids[i][j-1];
					margin[i][j] = margin[i][j-1];
				} else {
					grids[i][j] = grids[i][j+1];
					prm[i] = grids[i][j+1];
					margin[i][j] = pmf->neglpost ( prm, data );
					if ( std::isinf ( margin[i][j] ) ) std::cerr << "WARNING: margin("<<i<<","<<j<<") is inf\n";
				}
			} else margin[i][j] = (p > -1e10 ? p : -1e10 );
			if ( std::isinf ( margin[i][j] ) ) margin[i][j] = 1e20;
#ifdef DEBUG_INTEGRATE
			std::cerr << " margin = " << margin[i][j] << "\n";
#endif

			// And compute average Z online
			if (p>-1e10 && p<1e10) {
				Z += ( margin[i][j] );
				count ++;
			}
		}
		Z /= count;

		// Include Z for numerical stability
		for ( j=0; j<gridsize; j++ ) {
			margin[i][j] = exp ( Z - margin[i][j] );
		}

	}

	for ( i=0; i<nprm; i++ ) {
		normalize_probability ( grids[i], margin[i] );
		switch (i) {
			case 0:
				distparams[i] = match_gauss ( grids[i], margin[i] );
				fitted_posteriors[i] = new GaussPrior ( distparams[i][0], distparams[i][1] );
				break;
			case 1:
				distparams[i] = match_gamma ( grids[i], margin[i] );
				fitted_posteriors[i] = new GammaPrior ( distparams[i][0], distparams[i][1] );
				break;
			case 2: case 3:
				distparams[i] = match_beta ( grids[i], margin[i] );
				fitted_posteriors[i] = new BetaPrior ( distparams[i][0], distparams[i][1] );
				break;
		}
	}

	PsiIndependentPosterior out ( nprm, fitted_posteriors, grids, margin );

	for ( i=0; i<nprm; i++ ) {
		// delete fitted_posteriors[i];
	}

	return out;
}

MCMCList sample_posterior (
		const PsiPsychometric *pmf,
		const PsiData *data,
		PsiIndependentPosterior& post,
		unsigned int nsamples,
		unsigned int propose
		)
{
	unsigned int nprm ( pmf->getNparams() ), i, j, k;
	unsigned int nproposals ( nsamples*propose );
	MCMCList finalsamples ( nsamples, nprm, data->getNblocks() );
	double q,p,q_raw,p_raw;
	double nduplicate ( 0 );
	PsiRandom rng;
	PsiPrior * posteri;
	std::vector < PsiPrior* > posteriors ( nprm );
	double H(0),N(0);

	std::vector< std::vector<double> > proposed ( nproposals, std::vector<double> (nprm) );
	std::vector<double> weights ( nproposals );
	std::vector<double> cum_probs ( nproposals );
	std::vector<double> rnumbers ( nsamples );

	for ( j=0; j<nprm; j++ ) {
		posteriors[j] = post.get_posterior (j);
	}

	for ( i=0; i<nproposals; i++ ) {
		// Propose
		for ( j=0; j<nprm; j++ )
			proposed[i][j] = posteriors[j]->rand();
		// determine weight
		q = 1.;
        for ( j=0; j<nprm; j++ ) {
			posteri = post.get_posterior(j);
			q_raw = posteri->pdf ( proposed[i][j] );
			if ( q_raw > 1e10 )
				q_raw = 1e10;
			if ( q_raw != q_raw )
				q_raw = 1e5;
			if ( q_raw<1e-5 )
				q_raw = 1e-5;
            q *= q_raw;
			delete posteri;
		}
        p = - pmf->neglpost ( proposed[i], data );
		if ( std::isinf ( p ) || p!=p )
			weights[i] = 0;
		else
			weights[i] = exp ( p - log (q) );
#ifdef DEBUG_INTEGRATE
		if ( weights[i] != weights[i] || weights[i] > 1e10) {
			std::cerr << "Index: " << i << ", weight: " << weights[i] << ", p: " << p << ", q: " << q << ", th: (" << proposed[i][0];
			for ( j=1; j<nprm; j++ ) std::cerr << ", " << proposed[i][j];
			std::cerr << ")\n";
			std::cerr.flush();
		}
#endif

		// make a cumulative distribution vector for the weights
		if (i>0)
			cum_probs[i] = cum_probs[i-1] + weights[i];
		else {
			cum_probs[0] = weights[0];
#ifdef DEBUG_INTEGRATE
			std::cerr << "w0 = " << weights[0] << "\n";
#endif
		}
	}

	for ( i=0; i<nsamples; i++ ) {
		// And generate random numbers
		rnumbers[i] = rng.rngcall();
	}

	// Normalize the cumulative distribution
	for ( i=0; i<nproposals; i++ )
		cum_probs[i] /= cum_probs[nproposals-1];

	// Calculate entropy for diagnosis
	if ( cum_probs[0] > 0 )
		H = - cum_probs[0] * log(cum_probs[0]);
	else
		H = 0;
#ifdef DEBUG_INTEGRATE
	std::cerr << "H=" << H << ", cum_probs[0] = " << cum_probs[0] << "\n";
#endif
	N = 1.;
	// Avoid zeros
	for ( i=0; i<nproposals-1; i++ ) {
		if ((cum_probs[i+1]-cum_probs[i])>0) {
			H -= (cum_probs[i+1]-cum_probs[i]) * log ( cum_probs[i+1]-cum_probs[i] );
			N += 1;
		}
#ifdef DEBUG_INTEGRATE
		if ( H!=H ) {
			std::cerr << "H became nan in the " << i << "th iteration, cum_probs["<<i+1<<"] = " << cum_probs[i+1] << ", cum_probs["<<i<<"] = " << cum_probs[i] <<"\n";
			break;
		}
#endif
	}
	H /= log(N);
#ifdef DEBUG_INTEGRATE
	std::cerr << "H = " << H << "\n";
#endif

	sort ( rnumbers.begin(), rnumbers.end() );

	// resampling
	i = j = 0;
    while (i<nsamples) {
        k = 0;
        while (rnumbers[i] <= cum_probs[j]) {
			finalsamples.setEst ( i, proposed[j], pmf->deviance ( proposed[j], data ) );
			nduplicate += k;
			k=1;
            i++;
            if (i>=nsamples)
                break;
		}
        j++;
		if (j>nproposals) {
#ifdef DEBUG_INTEGRATE
			std::cerr << "What's going on here? i=" <<  i << ", cum_probs.max() = " << cum_probs[nproposals-1] << "\n";
#endif
			break;
		}
	}

	finalsamples.set_accept_rate ( double(nduplicate)/nsamples );
	finalsamples.set_entropy ( H );

	for ( i=0; i<nprm; i++ )
		delete posteriors[i];

	return finalsamples;
}

void sample_diagnostics (
		const PsiPsychometric *pmf,
		const PsiData *data,
		MCMCList *samples
		)
{
	unsigned int i,j,k, nprm ( pmf->getNparams() ), nblocks ( data->getNblocks() );
	std::vector<double> probs ( nblocks );
	std::vector<double> est ( nprm );
	PsiData *localdata = new PsiData ( data->getIntensities(), data->getNtrials(), data->getNcorrect(), data->getNalternatives() );
	std::vector<int> posterior_predictive ( nblocks );

	std::vector<double> reducedx ( data->getNblocks()-1 );
	std::vector<int> reducedk ( data->getNblocks()-1 );
	std::vector<int> reducedn ( data->getNblocks()-1 );
	std::vector< PsiData* > reduceddata ( data->getNblocks() );

	for ( i=0; i<nblocks; i++ ) {
		j = 0;
		for (k=0; k<nblocks; k++ ) {
			if ( i!=k ) {
				reducedx[j] = data->getIntensity(k);
				reducedk[j] = data->getNcorrect(k);
				reducedn[j] = data->getNtrials(k);
				j++;
			}
		}
		reduceddata[i] = new PsiData ( reducedx, reducedn, reducedk, data->getNalternatives() );
	}

	for ( i=0; i<samples->getNsamples(); i++ ) {
		for ( j=0; j<nprm; j++ )
			est[j] = samples->getEst ( i, j );

		for ( j=0; j<nblocks; j++ )
			probs[j] = pmf->evaluate ( data->getIntensity(j), est );
		newsample ( localdata, probs, &posterior_predictive );
		localdata->setNcorrect ( posterior_predictive );
		samples->setppData ( i, posterior_predictive, pmf->deviance ( est, localdata ) );

		probs = pmf->getDevianceResiduals ( est, data );
		samples->setRpd ( i, pmf->getRpd ( probs, est, data ) );
		samples->setRkd ( i, pmf->getRkd ( probs, data ) );

		probs = pmf->getDevianceResiduals ( est, localdata );
		samples->setppRpd ( i, pmf->getRpd ( probs, est, localdata ) );
		samples->setppRkd ( i, pmf->getRkd ( probs, localdata ) );

		// Store log posterior ratios for reduced data sets
		for ( j=0; j<nblocks; j++ )
			samples->setlogratio ( i, j, pmf->neglpost(est,data) - pmf->neglpost(est,reduceddata[j]) );
	}

	for ( i=0; i<nblocks; i++ ) {
		delete reduceddata[i];
	}
	delete localdata;
}
