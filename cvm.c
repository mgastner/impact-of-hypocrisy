/******************************** Inclusions. ********************************/
#include "hypocrisy.h"

/**************************** Function prototypes. ***************************/
double xlogx (double x);
double predicted_mean_cvm (int n, int nr_ext_init, int nr_int_init,
		           double rc, double re, double ri);
  
/*************** Utility functions for the mean consensus time. **************/
double xlogx (double x) {
  return (x * log(x));
}
double predicted_mean_cvm (int n, int nr_ext_init, int nr_int_init,
		           double rc, double re, double ri)
{
  double rho_ext_init = (double) nr_ext_init / n;
  double rho_int_init = (double) nr_int_init / n;
  double m = (ri*rho_ext_init + re*rho_int_init) / (re + ri);
  double prefac = (rc+re+ri)*(re+ri)*(re+ri) * n /
    (rc * ri * ((re+ri)*(re+ri) + rc*ri));
  
  return (-prefac * (xlogx(m) + xlogx(1-m)));
}

/*********************************** Main. ***********************************/
int main (int argc, char* argv[])
{
  double rc,         /* Rate of copying a neighbour. */
    re,              /* Rate of externalizing. */
    ri,              /* Rate of internalizing. */
    cons_time1 = 0,  /* Sum of consensus times. */
    cons_time2 = 0,  /* Sum of squared consensus times. */
    se_mean;         /* Standard error of the mean consensus time. */
  FILE *traj_file = fopen("dynamics.dat", "w");
  gsl_rng *rng;
  int n,          /* Number of agents. */
    nr_ext_init,  /* Initial number of agents whose external opinion is red. */
    nr_int_init,  /* Initial number of agents whose internal opinion is red. */
    nrr_init,     /* Initial number of agents whose opinion is red, both */
                  /* externally and internally. */
    n_run,        /* Number of runs. */
    n_red_wins;   /* Number of red victories. */
  long seed = 0;
  
  /************************** Parameter checks. ******************************/
  if (argc<9 || argc>10) {
    fprintf(stderr,"ERROR: Incorrect number of command-line arguments.\n");
    fprintf(stderr," Use as: ./cvm n_run rc re ri n ");
    fprintf(stderr, "nr_ext_init nr_int_init nrr_init (seed)\n");
    exit(1);
  }
  n_run = atoi(argv[1]);  /* Number of independent runs. */
  rc = atof(argv[2]);     /* Rate of copying a neighbour. */
  re = atof(argv[3]);     /* Rate of externalizing. */
  ri = atof(argv[4]);     /* Rate of internalizing. */
  n = atoi(argv[5]);      /* Number of agents. */
  nr_ext_init = atoi(argv[6]);  /* Initial number of external red opinions. */
  nr_int_init = atoi(argv[7]);  /* Initial number of internal red opinions. */
  nrr_init = atoi(argv[8]);     /* Initial number of agents whose opinion is */
                                /* red, both externally and internally. */
  if (nr_ext_init<0 || nr_ext_init>n || nr_int_init<0 || nr_int_init>n ||
      nrr_init < MAX(0, nr_ext_init + nr_int_init - n) ||
      nrr_init > MIN(nr_ext_init, nr_int_init)) {
    fprintf(stderr, "ERROR: Invalid initial conditions.\n");
    exit(1);
  }
  printf("n_run = %i, rc = %f, re = %f, ri = %f, n = %i\n",
  	 n_run, rc, re, ri, n);
  printf("nr_ext_init = %i, nr_int_init = %i, nrr_init = %i\n",
   	 nr_ext_init, nr_int_init, nrr_init);

  /******************** Start the random number generator. *******************/
  rng = gsl_rng_alloc(gsl_rng_mt19937);
  if (seed == 0) {
    seed = time(NULL);
  }
  printf("SEED OF RANDOM NUMBER GENERATOR: %li\n", seed);  
  gsl_rng_set(rng, seed);

  n_red_wins = 0;  /* Count how many times red wins. */

  /* During the first run, print the trajectory to a file. */
  fprintf(traj_file, "t,rho_R,rho_r,rho_Rr\n");
  fprintf(traj_file, "0.000000,%f,%f,%f\n",
	  (double) nr_ext_init / n, (double) nr_int_init / n,
	  (double) nrr_init / n);
  
  /******************* Outer loop: n_run independent runs. *******************/
  for (int run = 0; run < n_run; run++) {
    if (run % 250 == 0) {
      printf("working on run %i out of %i\n", run, n_run);
    }

    /* Number of agents whose external opinion is red. */
    int nr_ext = nr_ext_init;

    /* Number of agents whose internal opinion is red. */
    int nr_int = nr_int_init;

    /* Number of agents whose opinion is red, both externally and internally. */
    int nrr = nrr_init;
    double t = 0.0;  /* Current time. */

    /************* Inner loop: repeat until we reach a consensus. ************/
    while ((nr_ext!=0 || nr_int!=0) && (nr_ext!=n || nr_int!=n)) {

      /* Rate of (nr_ext, nr_int, nrr) -> (nr_ext + 1, nr_int, nrr). */
      double q_plus_ext = nr_ext * (n - nr_ext - nr_int + nrr) * rc / n;
      
      /* Rate of (nr_ext, nr_int, nrr) -> (nr_ext + 1, nr_int, nrr + 1). */
      double q_plus_ext_rr =  (nr_int - nrr) * (nr_ext * rc / n + re);
      
      /* Rate of (nr_ext, nr_int, nrr) -> (nr_ext, nr_int + 1, nrr + 1). */
      double q_plus_int_rr = ri * (nr_ext - nrr);

      /* Rate of (nr_ext, nr_int, nrr) -> (nr_ext - 1, nr_int, nrr). */
      double q_minus_ext = (nr_ext - nrr) * ((n - nr_ext) * rc / n + re);

      /* Rate of (nr_ext, nr_int, nrr) -> (nr_ext, nr_int - 1, nrr). */
      double q_minus_int = ri * (nr_int - nrr);
      
      /* Rate of (nr_ext, nr_int, nrr) -> (nr_ext - 1, nr_int, nrr - 1). */
      double q_minus_ext_rr = nrr * (n - nr_ext) * rc / n;
      
      /* Jump rate. */
      double lambda = q_plus_ext + q_plus_ext_rr + q_plus_int_rr +
   	q_minus_ext + q_minus_int + q_minus_ext_rr;
      
      /* Generate the waiting time until the jump. It is exponentially       */
      /* distributed with rate lambda. Note that the log() in the following  */
      /* expression is negative.                                             */
      t += gsl_ran_exponential(rng, 1.0 / lambda);
      
      /* Choose next state. */
      double rnd = lambda * gsl_rng_uniform(rng);

      /* Change number of agents who are externally and internally red. */
      if (rnd < q_plus_ext) {
       	nr_ext++;
      } else if (rnd < q_plus_ext + q_plus_ext_rr) {
 	nr_ext++;
   	nrr++;
      } else if (rnd < q_plus_ext + q_plus_ext_rr + q_plus_int_rr) {
   	nr_int++;
   	nrr++;
      } else if (rnd < q_plus_ext + q_plus_ext_rr + q_plus_int_rr +
		 q_minus_ext) {
   	nr_ext--;
      } else if (rnd < q_plus_ext + q_plus_ext_rr + q_plus_int_rr +
		 q_minus_ext + q_minus_int) {
   	nr_int--;
      } else {
  	nr_ext--;
  	nrr--;
      }

      /* During the first run, print the trajectory to a file. */
      if (run == 0)
	fprintf(traj_file, "%f,%f,%f,%f\n",
		t, (double) nr_ext / n, (double) nr_int / n, (double) nrr / n);
    }
    cons_time1 += t;
    cons_time2 += t * t;
    n_red_wins += (nr_ext == n);
  }

  /************************** Fraction of red wins. **************************/
  printf("observed fraction of red wins = %f\n", (double) n_red_wins / n_run);
  printf("expected fraction of red wins = %f\n", (double)
   	 (ri*nr_ext_init + re*nr_int_init) / (n * (re + ri)));

  /***************** Mean consensus time with error estimate. ****************/
  cons_time1 /= n_run;  /* First moment. */
  cons_time2 /= n_run;  /* Second moment. */
  se_mean =             /* Standard error of mean. */
    sqrt((cons_time2 - cons_time1 * cons_time1) / (n_run - 1));  
  printf("mean consensus time: %f +/- %f (95%% CI)\n",
	 cons_time1, 1.959964*se_mean);
  printf("predicted mean consensus time: %f\n",
	 predicted_mean_cvm(n, nr_ext_init, nr_int_init, rc, re, ri));
  
  return 0;
}
