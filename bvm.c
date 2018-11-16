/******************************** Inclusions. ********************************/
#include "hypocrisy.h"

/**************************** Function prototypes. ***************************/
double xlogx (double x);
double predicted_mean_bvm (int n, int nr_init, double rc);

/*************** Utility functions for the mean consensus time. **************/
double xlogx (double x) {
  return (x * log(x));
}
double predicted_mean_bvm (int n, int nr_init, double rc)
{
  double rho_r = (double) nr_init / n;
  return (-n * (xlogx(rho_r) + xlogx(1-rho_r)) / rc);
}

/*********************************** Main. ***********************************/
int main (int argc, char* argv[])
{
  double rc,         /* Rate of copying a neighbour. */
    cons_time1 = 0,  /* Sum of consensus times. */
    cons_time2 = 0,  /* Sum of squared consensus times. */
    se_mean;         /* Standard error of the mean consensus time. */
  FILE *traj_file = fopen("dynamics.dat", "w");
  gsl_rng *rng;
  int n,           /* Number of agents. */
    nr_init,       /* Initial number of agents whose opinion is red. */
    n_run,         /* Number of runs. */
    n_red_wins;    /* Number of red victories. */
  long seed = 0;

  /************************** Parameter checks. ******************************/
  if (argc < 5 || argc > 6) {
    fprintf(stderr,"ERROR: Incorrect number of command-line arguments.\n");
    fprintf(stderr," Use as: ./bvm n_run rc n nr_init (seed)\n");
    exit(1);
  }
  n_run = atoi(argv[1]);    /* Number of independent runs. */
  rc = atof(argv[2]);       /* Rate of copying a neighbour. */
  n = atoi(argv[3]);        /* Number of agents. */
  nr_init = atoi(argv[4]);  /* Initial number of red opinions. */
  if (nr_init<0 || nr_init>n) {
    fprintf(stderr, "ERROR: Invalid initial conditions.\n");
    exit(1);
  }
  printf("n_run = %d, rc = %f, n = %d, nr_init = %d\n",
	 n_run, rc, n, nr_init);

  /******************** Start the random number generator. *******************/
  rng = gsl_rng_alloc(gsl_rng_mt19937);
  if (seed == 0) {
    seed = time(NULL);
  }
  printf("SEED OF RANDOM NUMBER GENERATOR: %li\n", seed);  
  gsl_rng_set(rng, seed);

  n_red_wins = 0;  /* Count how many times red wins. */

  /* During the first run, print the trajectory to a file. */
  fprintf(traj_file, "t,rho_R\n");
  fprintf(traj_file, "0.000000,%f\n", (double) nr_init / n);
  
  /******************* Outer loop: n_run independent runs. *******************/
  for (int run = 0; run < n_run; run++) {
    if (run % 250 == 0) {
      printf("working on run %d out of %d\n", run, n_run);
    }
    int nr = nr_init;  /* Number of agents whose opinion is red. */
    double t = 0.0;    /* Current time. */

    /************* Inner loop: repeat until we reach a consensus. ************/
    while (nr!=0 && nr!=n) {

      /* Rate of nr -> nr + 1. */
      double q_plus = nr * (n - nr) * rc / n;
      
      /* Rate of nr -> nr - 1. */
      double q_minus = nr * (n - nr) * rc / n;
      
      /* Jump rate. */
      double lambda = q_plus + q_minus;

      /* Generate the waiting time until the jump. It is exponentially       */
      /* distributed with rate lambda. Note that the log() in the following  */
      /* expression is negative.                                             */
      t += gsl_ran_exponential(rng, 1.0 / lambda);

      /* Choose next state. */
      double rnd = lambda * gsl_rng_uniform(rng);

      /* Change number of red agents. */
      nr += 2*(rnd < q_plus) - 1;

      /* During the first run, print the trajectory to a file. */
      if (run == 0)
	fprintf(traj_file, "%f,%f\n", t, (double) nr / n);
    }
    cons_time1 += t;
    cons_time2 += t * t;
    n_red_wins += (nr == n);
  }

  /************************** Fraction of red wins. **************************/
  printf("observed fraction of red wins = %f\n", (double) n_red_wins / n_run);
  printf("expected fraction of red wins = %f\n", (double) nr_init / n);

  /***************** Mean consensus time with error estimate. ****************/
  cons_time1 /= n_run;  /* First moment. */
  cons_time2 /= n_run;  /* Second moment. */
  se_mean =             /* Standard error of mean. */
    sqrt((cons_time2 - cons_time1 * cons_time1) / (n_run - 1));  
  printf("mean consensus time: %f +/- %f (95%% CI)\n",
	 cons_time1, 1.959964*se_mean);
  printf("predicted mean consensus time: %f\n",
	 predicted_mean_bvm(n, nr_init, rc));

  return 0;
}
