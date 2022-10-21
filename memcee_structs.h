struct McObject{
  double *values;
  //int nparams;
  double logL;
  int unaccepted;

};
  
struct McParams{
  int samplesize; 
  int walkers;// preferably even integer
  //double kfactor; //scale the scaled cov with this 
  //int maximum_iterations;
  //double exit_factor;
  char *outfile_all;
  char *outfile_select;
  //char *livefile;
  //int checkpoint;
  int checkpoint_continue;
};

struct mc_thread_params {
  double (*loglike)(double *,void *, int); 
  void *loglike_params;
  double *values;
  double *logL;
  struct McObject **Samples;
  int *left; 
  int *right;
  int *nparams;
  int *tsamp;
  struct McObject *walker;
  struct McObject *newwalker;
  int mc_threadnum;
};



