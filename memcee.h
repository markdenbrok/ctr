int memcee(struct McParams, int, double (*)(double *,void *, int), void *, int);

void memcee_initialize_sample(struct McObject **,int, double (*)(double *,void *, int), void *, int, int);

void *thread_eval_ll(void *);

void memcee_update_samples(struct McObject **,int, int, double (*)(double *,void *,int), void *, int, int);

void memcee_update_subsample(int, int, struct McObject **,int,int,int,int, double (*)(double *,void *,int), void *, int, int);

double memcee_gz(double);

void memcee_generate_new_point(struct McObject *, struct McObject *, int , struct McObject **, int , int , double (*)(double *,void *,int), void *, int, int );

int num_lines_in_file(char *);

void memcee_initialize_sample_from_file(struct McObject **,int, int, int, char *);

int can_open_file(char *);

void write_output_to_file(struct McParams*,int ,struct McObject **, int, int);

