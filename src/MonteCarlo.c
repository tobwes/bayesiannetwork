/**
   Monte Carlo methods:
   - importance sampling
   - rejection sampling
   - metropolis sampling
   - (TODO: gibbs sampling) 
*/
#include"MonteCarlo.h"

MonteCarloEngine* initializeMonteCarloEngine(ProbabilitySpace* ps, Density* density, Density* sampling_prototype, int sampling_algorithm, int burnin, int gap, double rejection_sampling_factor){
  MonteCarloEngine* mce = ps->mce;

  if(mce == 0) mce = ps->mce = malloc(sizeof(MonteCarloEngine));

  mce->density = density;
  mce->sampling_prototype = sampling_prototype;
  mce->rejection_sampling_factor = rejection_sampling_factor;
  
  mce->burnin = burnin;
  mce->gap = gap;

  mce->sampling_algorithm = sampling_algorithm;

  mce->end = mce->variables = malloc(sizeof(double)*NPARALLEL*ps->dim*(burnin+gap*MAX_RUNS));

  return mce;
}

void destroyMonteCarloEngine(MonteCarloEngine* mce){
  if(mce->sampling_prototype) free(mce->sampling_prototype);
  if(mce->variables) free(mce->variables);
}

void clear(MonteCarloEngine* mce){
  mce->end = mce->variables;
}

const double* getRVs(ProbabilitySpace* ps, int n, bool generate_new){
  MonteCarloEngine* mce = ps->mce;
  int m;
  if(generate_new){
    m = n;
  }else if( (m = NPARALLEL*mce->burnin + n - (mce->end-mce->variables)/ps->dim) < 0 ){
    m = 0;
  }

  if(m>0){ drawRVs(ps, m); }

  return mce->end - n*ps->dim;
  
}

void new_burnin(ProbabilitySpace* ps, bool new_starting_point){
  MonteCarloEngine* mce = ps->mce;

  if(mce->variables == mce->end && mce->burnin != 0 && !new_starting_point){
    printf("*** Critical error: cannot burn in without a starting point!\n");
    exit(0);
  }else if(!new_starting_point){
    /** If not choosing a new starting point, make sure that there is enough memory for the burn in (by starting at the beginning of the array) **/
    memcpy(mce->variables,mce->end-ps->dim*NPARALLEL,sizeof(double)*ps->dim*NPARALLEL);
    mce->end = mce->variables + NPARALLEL*ps->dim;
  }
  
  for(int i=0; i<NPARALLEL; ++i){
    burnin(ps, mce->end+i*ps->dim,new_starting_point); // on first run: burn in
  }
  mce->end = mce->end + ps->dim * mce->burnin * NPARALLEL;
}

void drawRVs(ProbabilitySpace* ps, int n){
  MonteCarloEngine* mce = ps->mce;

  if(mce == 0){
    printf("Critical error: Monte Carlo Engine needs to be set up before drawing RVs\n");
    exit(0);
  }
  
  if(mce->end == mce->variables && mce->burnin!=0){
    new_burnin(ps,true);
  }

  int array_size = NPARALLEL*(MAX_RUNS*mce->gap+mce->burnin)*ps->dim;

  if( n*ps->dim+NPARALLEL*mce->burnin*ps->dim > array_size ){
    printf("*** Critical error: Monte Carlo Engine buffer overflow\n");
    exit(0);
  } else if( mce->variables + array_size < mce->end + n*ps->dim  ){
#if VERB == 3
    printf("MonteCarlo.c: reallocating RV cache...");
#endif
    memcpy(mce->variables,mce->end - NPARALLEL*ps->dim*100,sizeof(double)*NPARALLEL*ps->dim*100);
    mce->end = mce->variables+100*ps->dim*NPARALLEL;
  }

  for(int i=0; i<NPARALLEL; ++i){
    if(i<n%NPARALLEL){
      start_sampling(ps,(int)(n/NPARALLEL)+1,mce->end + i*ps->dim);
    }else{
      start_sampling(ps,(int)(n/NPARALLEL),mce->end + i*ps->dim);
    }
  }
  
  mce->end += ps->dim*n;
  
}

void burnin(ProbabilitySpace* ps, double* start, bool new_starting_point){
  MonteCarloEngine* mce = ps->mce;

  if(new_starting_point){
      /** choosing random starting point **/
    for(int i=0; i<ps->dim; i++){
      start[i] = ((double)rand()/(double)RAND_MAX)*STARTING_WINDOW_SIDE_LENGTH;
    }
    start += NPARALLEL*ps->dim;
  }
  
  /**starting algorithm**/
  start_sampling(ps,mce->burnin-1,start);
}

void start_sampling(ProbabilitySpace* ps, int m, double* start){
  MonteCarloEngine* mce = ps->mce;
  
  if(mce->sampling_algorithm == IMPORTANCE_SAMPLING)
    importance_sampling(ps,m,start);
  else if(mce->sampling_algorithm == REJECTION_SAMPLING)
    rejection_sampling(ps,m,start);
  else if(mce->sampling_algorithm == METROPOLIS_SAMPLING)
    metropolis_sampling(ps,m,start);
  //else if(mce->sampling_algorithm == GIBBS_SAMPLING)
  else
    printf("*** Critical error: unknown sampling technique");
}  

/** Importance sampling**/

/**
   Importance sampling
   Given:
   - density function p*, not necessarily normalized
   - probability density function p=cp*, c unknown
   Steps:
   - choose a sampling distribution Q with density q, "close" to p
   - use the following formula to approximate any integral, where \{\xi_i\} are drawn from distribution Q

   \int f(x)p(x)dx &= \int f(x)\frac{p(x)}{q(x)}q(x)dx\\
   &= \frac{ \int f(x)\frac{p*(x)}{q(x)}q(x)dx }{ \int p*(x) dx }
   &= \frac{ \int f(x)\frac{p*(x)}{q(x)}q(x)dx }{ \int \frac{p*(x)}{q(x)}q(x)dx }
   &\tilde = \frac{ \sum_{i=1}^n f(\xi_i)\frac{p*(\xi_i)}{q(\xi_i)} }{\sum_{i=1}^n \frac{p*(\xi_i)}{q(\xi_i)}}
**/

void importance_sampling(ProbabilitySpace* ps, int m, double* start){
  MonteCarloEngine* mce = ps->mce;

  while(m>0){
    mce->sampling_prototype->rnd(mce->sampling_prototype, start);
    start += NPARALLEL*ps->dim;
    m--;
  }
}

/** Rejection sampling **/

/**
   Rejection sampling:
   Given:
   - one dimensional density function p*, not necessarily normalized
   - probability density function p=cp*, c unknown
   Steps:
   - choose a sampling distribution Q with density q, such that for some known k, kq(x)>p*(x) for all x
   - draw random samples x from Q
   - draw random number in the interval [0,kq(x)]
   - reject x, if u>p*(x)
**/

void rejection_sampling(ProbabilitySpace* ps, int m, double* start){
  MonteCarloEngine* mce = ps->mce;

  if(mce->sampling_prototype == 0){
    printf("*** Critical error: need to setup sampling density first");
    exit(0);
  }
  
  int abort = MAX_REJECTION_TRIALS;
  while(m>0){
    mce->sampling_prototype->rnd(mce->sampling_prototype, start);
    double u = (double)rand()/(double)RAND_MAX;
    if( mce->density->evaluate(mce->density,start,start+ps->dim-1) / ( mce->rejection_sampling_factor * mce->sampling_prototype->evaluate(mce->sampling_prototype,start,start+ps->dim-1) ) > u ){
      start += NPARALLEL*ps->dim;
      m--;
    }else if( (abort--) == 0 ){
      printf("*** Critical error: too many trials were rejected.\n");
      exit(0);
    }
  }
}

/** Metropolis sampling **/

/**
   Given:
   - density function p*, not necessarily normalized
   - probability density function p=cp*, c unknown
   Steps:
   - generate a random starting position x_0 (this is done in burnin())
   - choose a set of proposal transition distributions \{Q(x,x')\} with densities \{q(x,x')\} to represent the probability to go from x' to x
   - for given x_t, generate a sample x_{t+1} according to the distribution Q(x,x_{t+1})
   - accept x_{t+1} if 
           a = \frac{ p*(x_{t+1}) q(x_t,x_{t+1}) }{ p*(x_{t}) q(x_{t+1},x_{t}) } >= 1, 
     otherwise accept x_{t+1} with probability a.
   Remarks:
   - let L be the length scale of p
   - let \epsilon be the length scale of q
   - to get independent samples, one should choose a gap of at least (L/\epsilon)^2 between the used sample variables (this is controversial)
**/

void metropolis_sampling(ProbabilitySpace* ps, int m, double* start){
  MonteCarloEngine* mce = ps->mce;

  while(m>0){
    double* previous = start - NPARALLEL*ps->dim; // exists because of burnin

    mce->sampling_prototype->setTransitionPoint(mce->sampling_prototype, previous, previous+ps->dim-1);
    mce->sampling_prototype->rnd(mce->sampling_prototype,start);

    double (*evaluateDensity)(Density*,const double*,const double*);
    double (*evaluatePrototype)(Density*,const double*,const double*);
    bool log = false;
    if(mce->density->evaluateLog && mce->sampling_prototype->evaluateLog){
      log = true;
      evaluateDensity = mce->density->evaluateLog;
      evaluatePrototype = mce->sampling_prototype->evaluateLog;
    }else{
      //printf("** WARNING: Rounding errors may be large...\n");
      evaluateDensity = mce->density->evaluate;
      evaluatePrototype = mce->sampling_prototype->evaluate;
    }

    double probability_previous_state,
      probability_new_state,
      transition_probability_walk_back,
      transition_probability_walk_forward;

    probability_new_state = evaluateDensity(mce->density,start,start+ps->dim-1);
    probability_previous_state = evaluateDensity(mce->density,previous,previous+ps->dim-1);

    transition_probability_walk_forward = evaluatePrototype(mce->sampling_prototype,start,start+ps->dim-1);
    mce->sampling_prototype->setTransitionPoint(mce->sampling_prototype, start, start+ps->dim-1);
    transition_probability_walk_back = evaluatePrototype(mce->sampling_prototype,previous,previous+ps->dim-1);
    
    double a;
    if(log){
      a = exp( probability_new_state - probability_previous_state + transition_probability_walk_back - transition_probability_walk_forward );
    } else {
      a = (probability_new_state / probability_previous_state) * (transition_probability_walk_back / transition_probability_walk_forward);
    }
    
    if( a >= 1 || ((double)rand()/(double)RAND_MAX < a) ){
      start += NPARALLEL*ps->dim;
      m--;
    }else{
      //printf("Proposal rejected...\n"); 
      memcpy(start,previous,sizeof(double)*ps->dim);
      start += NPARALLEL*ps->dim;
      m--;
    }

    //char weights[2048];
    //vectorToString(start-NPARALLEL*ps->dim,ps->dim,&weights);
    //printf("New sample: w=%s",weights);
  }
}

/** Gibbs sampling **/

/**
    Gibbs sampling:
    Given:
    - density function p*, not necessarily normalized, living on R^n, n>1
    - probability density function p=cp*, c unknown
    - For x\in R^n use the notation x=(x^{(1)},...,x^{(n)})
    Steps:
    - generate a random starting position x_0
    - choose a set of one dimensional proposal transition distributions \{Q^{(1)}(x^{(1)},x'^{(1)}),...,Q^{(n)}(x^{(n)},x'^{(n)})\} to represent the probabilities for the components of x to go from x^{(i)} to x'^{(i)}
    - do one Metropolis step for each component x^{(i)}, i=1..n
**/

/** adaptive Monte Carlo Method (AM) **/

/**
   AM:
   Given:
   - density function p*, not necessarily normalized, living on a bounded subset of R^n, n>1
   - probability density function p=cp*, c unknown
   - gaussian sampling_prototypes with mean at x_t
   Steps:
   - do a Metropolis step
   - after an initial warmup period, update the covariance of sampling_prototype's covariance matrix according to:
   - C_{t+1} = \frac{t-1}{t} C_t + \frac{s_d}{t} ( t m(X_{t-1})m(X_{t-1})^T - (t+1) m(X_t)m(X_t)^T + \epsilon I_d )
 **/
