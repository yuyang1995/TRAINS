#ifndef _PROTO_H
#define _PROTO_H

void begin_run();
void end_run();

void read_parset();
void check_param_range(char *name, int flag, int min, int max);
void get_num_particles(char *fname);

void read_parset_pt();
void read_data_pt();
void alloc_data_pt();
void free_data_pt();

void init();
void init_num_params();
void init_func();
void allocate_memory();
void free_memory();
void set_par_range_model();
void set_par_fix_model();
void set_par_prior_model();
void set_par_fix_str(int *label, double *value, const char *label_str, const char *value_str, int num, const char *tag);
UPDATE_FLAG get_update_flag(int pos);

void set_source_range_model(double **range);
void output_source_range_model(double *mean, double *std);
void input_source_range_model(double **range, int size);

void get_omega(const double *psource);
void residual_cal(const void *pm, PSRmodel *p, double **t, double **res, double *phase, double *dis, int is, int jp);
void residual_sum(double **res_src, double **res, int jp);
void residual_shift(double **res, int jp);
double LLR_Mx_Av(const void *pm, PSRmodel *p, double **t, double **res_data, double *phase);
void LLR_initial(int Nt);
void LLR_end();

void sim();
void sim_init();
void sim_end();
void set_source_value_sim(double *pbh);

void from_prior_gen(void *model);
void print_particle_gen(FILE *fp, const void *model);
void read_particle_gen(FILE *fp, void *model);
double perturb_gen(void *model, int which);
double log_likelihoods_cal_exam(const void *model);
void print_par_names_gen(char *fname);
void restart_action_gen(int iflag);

void sample_stats(double *sample, double *mean, double *std, int num_ps);
void postprocess_gen();
void reconstruct_gen();

void output_rec_pt(void *posterior_sample, int size_of_modeltype, int num_ps);
void output_pre_pt(const void *model);
void reconstruct_pt_init();
void reconstruct_pt_end();
double prob_initial_pt(const void *model);
double prob_pt(const void *model);
int dnest_pt(int argc, char **argv);
void accept_action_pt();
void kill_action_pt(int i, int i_copy);

int command_line_options(int argc, char **argv);
void print_version();
void fprint_version(FILE *fp);
void print_help();

#endif