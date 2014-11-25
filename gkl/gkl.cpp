/*
GKL S4 & S6
Credits:
Rolf Simoes
Eduardo Dias
Based on lectures of Prof. J. Ricardo Mendonça (Programa em Modelagem de Sistemas Complexos/EACH/USP) and
  on paper Gach, Kurdyumov, Levin (1978).
*/

#include <iostream>
#include <random>
#include <iomanip>

#define __GKL_MAX_STATES 6
#define __GKL_MAX_LEN 5000
#define __GKL_MAIN_STATES 2

/* data type definitions */
enum __GKL_CONVERGENCE {
	__gkl_convergence_none = 0,
	__gkl_convergence_right = 1,
	__gkl_convergence_wrong = -1
};

enum __GKL_PRINT_DETAIL {
	__gkl_print_detail_x0 = 1 << 1,
	__gkl_print_detail_x1 = 1 << 2,
	__gkl_print_detail_x2 = 1 << 3,
	__gkl_print_detail_z0 = 1 << 4,
	__gkl_print_detail_z1 = 1 << 5,
	__gkl_print_detail_convergence = 1 << 6
};

#define __GKL_PAST_MEMORY 10
int _past_rho[__GKL_PAST_MEMORY][__GKL_MAIN_STATES];

/* begin global variables definitions */
bool stop_on_convergence = false;
bool x0_informed = false;
bool x1_informed = false;
bool x2_informed = false;
bool z0_informed = false;
bool z1_informed = false;
bool seed_informed = false;
bool max_steps_informed = false;
bool print_out_rho = false;
bool print_out_ca = false;
bool print_out_stats = false;
bool print_ca_mapped = true;
int num_states = 4;
int ca_len = 400;
int max_steps = 4000;
int simulations = 1;
int stats_from = 0;
int stop_noise_on_step = 0;
double x0 = 0.0;
double x1 = 0.0;
double x2 = 0.0;
double z0 = 0.0;
double z1 = 0.0;
double noise = 0.0;

int simulation = 0;
int step = 0;
double entropy_before_noise = 0.0;
double entropy_after_noise = 0.0;
__GKL_CONVERGENCE convergence = __gkl_convergence_none;
int ca[2][__GKL_MAX_LEN + 2] = { 0 };
int psi[__GKL_MAX_STATES][__GKL_MAX_STATES][__GKL_MAX_STATES] = { 0 };
int rho[__GKL_MAX_STATES] = { 0 };
int s_count[__GKL_MAX_STATES] = { 0 };
int zeta[__GKL_MAX_STATES][__GKL_MAX_LEN + 1] = { 0 };
double sum_rho[__GKL_MAX_STATES] = { 0 };
double sum_rho_sqr[__GKL_MAX_STATES] = { 0 };
/* end global variables definitions */

using namespace std;

/* begin defining random functions */
random_device __rd;
//mt19937 __gen(__rd());
//uniform_real_distribution<> random(0, 1);
//void __init_gen(){ ; }
/* 
Public domain code for JKISS RNG - The period of JKISS is aproximately 2**127 (MT's period is much larger: 2**19937-1) 
but KISS is about 40% more fast than MT PRNG. (Source: http://www0.cs.ucl.ac.uk/staff/d.jones/GoodPracticeRNG.pdf) 
*/
// KISS Seed variables;
unsigned int __KISS_x = 123456789, __KISS_y = 987654321, __KISS_z = 43219876, __KISS_c = 6543217;
// initialyzing KISS generator seed;
void __init_rand_gen(){
	if (!seed_informed){
		__KISS_x = __rd();
		while (!(__KISS_y = __rd())); /* y must not be zero! */
		__KISS_z = __rd();
		__KISS_c = __rd() % 698769068 + 1; /* Should be less than 698769069 */
	}
	cout << "#Seed: x=" << __KISS_x << "; y=" << __KISS_y << "; z=" << __KISS_z << "; c=" << __KISS_c << endl;
}
// main function KISS generator proposed by G. Marsaglia (see source above to references);
unsigned int __gen(){
	unsigned long long t;
	__KISS_x = 314527869 * __KISS_x + 1234567;
	__KISS_y ^= __KISS_y << 5; __KISS_y ^= __KISS_y >> 7; __KISS_y ^= __KISS_y << 22;
	t = 4294584393ULL * __KISS_z + __KISS_c;
	__KISS_c = static_cast<unsigned int>(t >> 32);
	__KISS_z = static_cast<unsigned int>(t);
	return __KISS_x + __KISS_y + __KISS_z;
}
// returns a random number between 0.0 (inclusive) and 1.0 (exclusive);
double random(unsigned int (*x)(void)){
	return (*x)() / 4294967296.0;
}
/* end defining random functions */

/* begin auxiliary functions */
void __usage(char *program_name){
	cout << "Usage: " << program_name << " <gkl_type> [-L <int>][-n <float>][-runs <int>]" << endl;
	cout << "  [-T <int>][-c][-rho][-ca][-nomap][-stats <int>][-seed <int>]" << endl;
	cout << "  [-x0 <float>][-x1 <float>][-x2 <float>][-z <float>]" << endl;
	cout << endl;
	cout << "<gkl_type>:    Type of GKL. Choose one of (S4, S6). Default S4." << endl;
	cout << "-L <int>:      CA size up to " << __GKL_MAX_LEN << ". Default 400;" << endl;
	cout << "-n <float>:    PCA noise level between 0.0 and 1.0. Default 0.0;" << endl;
	cout << "-rep <int>:    Repeat simulation. Default 1;" << endl;
	cout << "-T <int>:      Max CA steps. Default 4000;" << endl;
	cout << "-c:            Stop on convergence is reached;" << endl;
	cout << "-rho:          Prints rho (state densities) aint steps;" << endl;
	cout << "-ca:           Prints CA configuration over steps;" << endl;
	cout << "-nomap:        Prints CA configuration with states' numbers;" << endl;
	cout << "-stats <int>:  Supress -c, -rho and -ca options and prints out rho stats;" << endl;
	cout << "               You must inform -T greater than <int> to get stats;" << endl;
	cout << "For GKL S4:" << endl;
	cout << "  -x0 <float>: Initial s0/s1 proportion between 0.0 and 1.0. Default random;" << endl;
	cout << "  -x1 <float>: Initial s2/s3 proportion between 0.0 and 1.0. Default random;" << endl;
	cout << "  -z <float>:  Initial x0/x1 proportion between 0.0 and 1.0. Default random;" << endl;
	cout << "For GKL S6:" << endl;
	cout << "  -x0 <float>: Initial s0/s1 proportion between 0.0 and 1.0. Default random;" << endl;
	cout << "  -x1 <float>: Initial s2/s3 proportion between 0.0 and 1.0. Default random;" << endl;
	cout << "  -x2 <float>: Initial s4/s5 proportion between 0.0 and 1.0. Default random;" << endl;
	cout << "  -z <float>:  Initial x0/(x1 + x2) proportion between 0.0 and 1.0. Default random." << endl;
	cout << endl;
}
void __handle_options(int argc, char *argv[]){
	cout << "#Command: " << argv[0] << " " << argv[1] << " ";
	// obligatory paramters;
	if ((strcmp(argv[1], "s4") == 0) || (strcmp(argv[1], "S4") == 0)){ num_states = 4; }
	else if ((strcmp(argv[1], "s6") == 0) || (strcmp(argv[1], "S6") == 0)){ num_states = 6; }
	else {
		__usage(argv[0]);
		exit(1);
	}
	// optional parameters;
	for (int i = 2; i < argc; ++i){ 
		if ((strcmp(argv[i], "-L") == 0) && (argc - 1 > i)){
			cout << argv[i] << " " << argv[i + 1] << " ";
			ca_len = atol(argv[++i]);
			if (ca_len > __GKL_MAX_LEN){
				__usage(argv[0]);
				exit(1);
			}
		} 
		else if ((strcmp(argv[i], "-x0") == 0) && (argc - 1 > i)){
			cout << argv[i] << " " << argv[i + 1] << " ";
			x0 = atof(argv[++i]);
			x0_informed = true;
			if (x0 < 0 || x0 > 1){ 
				__usage(argv[0]);
				exit(1);
			}
		}
		else if ((strcmp(argv[i], "-x1") == 0) && (argc - 1 > i)){
			cout << argv[i] << " " << argv[i + 1] << " ";
			x1 = atof(argv[++i]);
			x1_informed = true;
			if (x1 < 0 || x1 > 1){
				__usage(argv[0]);
				exit(1);
			}
		}
		else if ((strcmp(argv[i], "-x2") == 0) && (argc - 1 > i)){
			cout << argv[i] << " " << argv[i + 1] << " ";
			x2 = atof(argv[++i]);
			x2_informed = true;
			if (x2 < 0 || x2 > 1){
				__usage(argv[0]);
				exit(1);
			}
		}
		else if ((strcmp(argv[i], "-z0") == 0) && (argc - 1 > i)){
			cout << argv[i] << " " << argv[i + 1] << " ";
			z0 = atof(argv[++i]);
			z0_informed = true;
			if (z0 < 0 || z0 > 1){
				__usage(argv[0]);
				exit(1);
			}
		}
		else if ((strcmp(argv[i], "-z1") == 0) && (argc - 1 > i)){
			cout << argv[i] << " " << argv[i + 1] << " ";
			z1 = atof(argv[++i]);
			z1_informed = true;
			if (z1 < 0 || z1 > 1){
				__usage(argv[0]);
				exit(1);
			}
		}
		else if ((strcmp(argv[i], "-n") == 0) && (argc - 1 > i)){
			cout << argv[i] << " " << argv[i + 1] << " ";
			noise = atof(argv[++i]);
			if (noise < 0 || noise > 1){ 
				__usage(argv[0]);
				exit(1);
			}
		}
		else if ((strcmp(argv[i], "-rep") == 0) && (argc - 1 > i)){
			cout << argv[i] << " " << argv[i + 1] << " ";
			simulations = atol(argv[++i]);
			if (simulations < 1){
				__usage(argv[0]);
				exit(1);
			}
		}
		else if ((strcmp(argv[i], "-T") == 0) && (argc - 1 > i)){
			cout << argv[i] << " " << argv[i + 1] << " ";
			max_steps_informed = true;
			max_steps = atol(argv[++i]);
			if (max_steps < 0 || (print_out_stats && max_steps < stats_from)){
				__usage(argv[0]);
				exit(1);
			}
		}
		/*
		else if (strcmp(argv[i], "-c") == 0){
			cout << argv[i] << " ";
			stop_on_convergence = true;
		}
		else if (strcmp(argv[i], "-rho") == 0){
			cout << argv[i] << " ";
			print_out_rho = true;
		}
		else if (strcmp(argv[i], "-ca") == 0){
			cout << argv[i] << " ";
			print_out_ca = true;
		}
		else if (strcmp(argv[i], "-nomap") == 0){
			cout << argv[i] << " ";
			print_ca_mapped = false;
		}
		*/
		else if ((strcmp(argv[i], "-stats") == 0) && (argc - 1 > i)){
			cout << argv[i] << " " << argv[i + 1] << " ";
			print_out_stats = true;
			stats_from = atol(argv[++i]);
			if ((max_steps < stats_from) && max_steps_informed){
				__usage(argv[0]);
				exit(1);
			}
		}
		else if ((strcmp(argv[i], "-seed") == 0) && (argc - 4 > i)){
			cout << argv[i] << " " << argv[i + 1] << " " << argv[i + 2] << " " << argv[i + 3] << " " << argv[i + 4] << " ";
			seed_informed = true;
			__KISS_x = atol(argv[++i]);
			__KISS_y = atol(argv[++i]);
			__KISS_z = atol(argv[++i]);
			__KISS_c = atol(argv[++i]) % 698769068 + 1;
			if (!(__KISS_y > 0)){
				__usage(argv[0]);
				exit(1);
			}
		}
		else if (strcmp(argv[i], "-h") == 0){
			cout << argv[i] << " ";
			__usage(argv[0]);
			exit(0);
		}
		else {
			__usage(argv[0]);
			exit(1);
		}
	}
	cout << endl;
}
void __init_psi(){
	if (num_states == 4){
		/* begin initialize psi for gkl s4 */
		psi[0][0][0] = 0;    psi[0][0][1] = 0;    psi[0][0][2] = 0;    psi[0][0][3] = 0;
		psi[0][1][0] = 3;    psi[0][1][1] = 1;    psi[0][1][2] = 3;    psi[0][1][3] = 1;
		psi[0][2][0] = 0;    psi[0][2][1] = 2;    psi[0][2][2] = 0;    psi[0][2][3] = 0;
		psi[0][3][0] = 0;    psi[0][3][1] = 2;    psi[0][3][2] = 0;    psi[0][3][3] = 0;
		psi[1][0][0] = 3;    psi[1][0][1] = 3;    psi[1][0][2] = 3;    psi[1][0][3] = 3;
		psi[1][1][0] = 3;    psi[1][1][1] = 1;    psi[1][1][2] = 3;    psi[1][1][3] = 1;
		psi[1][2][0] = 2;    psi[1][2][1] = 1;    psi[1][2][2] = 2;    psi[1][2][3] = 2;
		psi[1][3][0] = 2;    psi[1][3][1] = 1;    psi[1][3][2] = 2;    psi[1][3][3] = 2;
		psi[2][0][0] = 3;    psi[2][0][1] = 3;    psi[2][0][2] = 3;    psi[2][0][3] = 3;
		psi[2][1][0] = 3;    psi[2][1][1] = 1;    psi[2][1][2] = 3;    psi[2][1][3] = 1;
		psi[2][2][0] = 2;    psi[2][2][1] = 1;    psi[2][2][2] = 2;    psi[2][2][3] = 2;
		psi[2][3][0] = 2;    psi[2][3][1] = 1;    psi[2][3][2] = 2;    psi[2][3][3] = 2;
		psi[3][0][0] = 0;    psi[3][0][1] = 0;    psi[3][0][2] = 0;    psi[3][0][3] = 0;
		psi[3][1][0] = 3;    psi[3][1][1] = 1;    psi[3][1][2] = 3;    psi[3][1][3] = 1;
		psi[3][2][0] = 2;    psi[3][2][1] = 1;    psi[3][2][2] = 2;    psi[3][2][3] = 2;
		psi[3][3][0] = 2;    psi[3][3][1] = 1;    psi[3][3][2] = 2;    psi[3][3][3] = 2;
		/* end initialize psi for gkl s4 */
	}
	else if (num_states == 6){
		/* begin initialize psi for gkl s6 */
		psi[0][0][0] = 0;    psi[0][0][1] = 5;    psi[0][0][2] = 0;    psi[0][0][3] = 5;
		psi[0][0][4] = 0;    psi[0][0][5] = 0;    psi[0][1][0] = 0;    psi[0][1][1] = 5;
		psi[0][1][2] = 4;    psi[0][1][3] = 1;    psi[0][1][4] = 1;    psi[0][1][5] = 1;
		psi[0][2][0] = 0;    psi[0][2][1] = 5;    psi[0][2][2] = 0;    psi[0][2][3] = 0;
		psi[0][2][4] = 0;    psi[0][2][5] = 0;    psi[0][3][0] = 0;    psi[0][3][1] = 5;
		psi[0][3][2] = 3;    psi[0][3][3] = 3;    psi[0][3][4] = 3;    psi[0][3][5] = 3;
		psi[0][4][0] = 0;    psi[0][4][1] = 5;    psi[0][4][2] = 2;    psi[0][4][3] = 2;
		psi[0][4][4] = 2;    psi[0][4][5] = 2;    psi[0][5][0] = 0;    psi[0][5][1] = 5;
		psi[0][5][2] = 3;    psi[0][5][3] = 3;    psi[0][5][4] = 3;    psi[0][5][5] = 3;
		psi[1][0][0] = 4;    psi[1][0][1] = 1;    psi[1][0][2] = 0;    psi[1][0][3] = 5;
		psi[1][0][4] = 0;    psi[1][0][5] = 0;    psi[1][1][0] = 4;    psi[1][1][1] = 1;
		psi[1][1][2] = 4;    psi[1][1][3] = 1;    psi[1][1][4] = 1;    psi[1][1][5] = 1;
		psi[1][2][0] = 4;    psi[1][2][1] = 1;    psi[1][2][2] = 2;    psi[1][2][3] = 2;
		psi[1][2][4] = 2;    psi[1][2][5] = 2;    psi[1][3][0] = 4;    psi[1][3][1] = 1;
		psi[1][3][2] = 1;    psi[1][3][3] = 1;    psi[1][3][4] = 1;    psi[1][3][5] = 1;
		psi[1][4][0] = 4;    psi[1][4][1] = 1;    psi[1][4][2] = 2;    psi[1][4][3] = 2;
		psi[1][4][4] = 2;    psi[1][4][5] = 2;    psi[1][5][0] = 4;    psi[1][5][1] = 1;
		psi[1][5][2] = 3;    psi[1][5][3] = 3;    psi[1][5][4] = 3;    psi[1][5][5] = 3;
		psi[2][0][0] = 4;    psi[2][0][1] = 4;    psi[2][0][2] = 2;    psi[2][0][3] = 1;
		psi[2][0][4] = 4;    psi[2][0][5] = 4;    psi[2][1][0] = 1;    psi[2][1][1] = 1;
		psi[2][1][2] = 2;    psi[2][1][3] = 1;    psi[2][1][4] = 1;    psi[2][1][5] = 1;
		psi[2][2][0] = 2;    psi[2][2][1] = 1;    psi[2][2][2] = 2;    psi[2][2][3] = 1;
		psi[2][2][4] = 2;    psi[2][2][5] = 2;    psi[2][3][0] = 0;    psi[2][3][1] = 3;
		psi[2][3][2] = 2;    psi[2][3][3] = 1;    psi[2][3][4] = 3;    psi[2][3][5] = 3;
		psi[2][4][0] = 2;    psi[2][4][1] = 2;    psi[2][4][2] = 2;    psi[2][4][3] = 2;
		psi[2][4][4] = 2;    psi[2][4][5] = 2;    psi[2][5][0] = 3;    psi[2][5][1] = 3;
		psi[2][5][2] = 2;    psi[2][5][3] = 3;    psi[2][5][4] = 3;    psi[2][5][5] = 3;
		psi[3][0][0] = 0;    psi[3][0][1] = 0;    psi[3][0][2] = 0;    psi[3][0][3] = 3;
		psi[3][0][4] = 0;    psi[3][0][5] = 0;    psi[3][1][0] = 5;    psi[3][1][1] = 5;
		psi[3][1][2] = 0;    psi[3][1][3] = 3;    psi[3][1][4] = 5;    psi[3][1][5] = 5;
		psi[3][2][0] = 2;    psi[3][2][1] = 1;    psi[3][2][2] = 0;    psi[3][2][3] = 3;
		psi[3][2][4] = 2;    psi[3][2][5] = 2;    psi[3][3][0] = 0;    psi[3][3][1] = 3;
		psi[3][3][2] = 0;    psi[3][3][3] = 3;    psi[3][3][4] = 3;    psi[3][3][5] = 3;
		psi[3][4][0] = 2;    psi[3][4][1] = 2;    psi[3][4][2] = 2;    psi[3][4][3] = 3;
		psi[3][4][4] = 2;    psi[3][4][5] = 2;    psi[3][5][0] = 3;    psi[3][5][1] = 3;
		psi[3][5][2] = 3;    psi[3][5][3] = 3;    psi[3][5][4] = 3;    psi[3][5][5] = 3;
		psi[4][0][0] = 0;    psi[4][0][1] = 0;    psi[4][0][2] = 0;    psi[4][0][3] = 5;
		psi[4][0][4] = 4;    psi[4][0][5] = 0;    psi[4][1][0] = 1;    psi[4][1][1] = 1;
		psi[4][1][2] = 4;    psi[4][1][3] = 1;    psi[4][1][4] = 4;    psi[4][1][5] = 1;
		psi[4][2][0] = 2;    psi[4][2][1] = 1;    psi[4][2][2] = 2;    psi[4][2][3] = 2;
		psi[4][2][4] = 4;    psi[4][2][5] = 2;    psi[4][3][0] = 0;    psi[4][3][1] = 3;
		psi[4][3][2] = 3;    psi[4][3][3] = 3;    psi[4][3][4] = 4;    psi[4][3][5] = 3;
		psi[4][4][0] = 2;    psi[4][4][1] = 2;    psi[4][4][2] = 2;    psi[4][4][3] = 2;
		psi[4][4][4] = 4;    psi[4][4][5] = 2;    psi[4][5][0] = 3;    psi[4][5][1] = 3;
		psi[4][5][2] = 3;    psi[4][5][3] = 3;    psi[4][5][4] = 4;    psi[4][5][5] = 3;
		psi[5][0][0] = 0;    psi[5][0][1] = 0;    psi[5][0][2] = 0;    psi[5][0][3] = 5;
		psi[5][0][4] = 0;    psi[5][0][5] = 5;    psi[5][1][0] = 1;    psi[5][1][1] = 1;
		psi[5][1][2] = 4;    psi[5][1][3] = 1;    psi[5][1][4] = 1;    psi[5][1][5] = 5;
		psi[5][2][0] = 2;    psi[5][2][1] = 1;    psi[5][2][2] = 2;    psi[5][2][3] = 2;
		psi[5][2][4] = 2;    psi[5][2][5] = 5;    psi[5][3][0] = 0;    psi[5][3][1] = 3;
		psi[5][3][2] = 3;    psi[5][3][3] = 3;    psi[5][3][4] = 3;    psi[5][3][5] = 5;
		psi[5][4][0] = 2;    psi[5][4][1] = 2;    psi[5][4][2] = 2;    psi[5][4][3] = 2;
		psi[5][4][4] = 2;    psi[5][4][5] = 5;    psi[5][5][0] = 3;    psi[5][5][1] = 3;
		psi[5][5][2] = 3;    psi[5][5][3] = 3;    psi[5][5][4] = 3;    psi[5][5][5] = 5;
		/* end initialize psi for gkl s6 */
	}
}
void __init_main(int argc, char *argv[]){
	// config output stream;
	cout << setiosflags(ios::fixed) << setprecision(6);
	// handling command-line options;
	__handle_options(argc, argv);
	// initialyzing random function (after __handle_options);
	__init_rand_gen();
	// initialyzing psi (after __handle_options);
	__init_psi();
}
void __init_simulation_parameters(){
	if (!x0_informed) { x0 = random(__gen); }
	if (!x1_informed) { x1 = random(__gen); }
	if (!x2_informed) { x2 = random(__gen); }
	if (!z0_informed) { z0 = random(__gen); }
	if (!z1_informed) { z1 = random(__gen); }
	/* begin reseting rho & rho stats */
	// used to calc rho stats if -stats option is informed;
	if (print_out_stats){
		for (int state = 0; state < __GKL_MAIN_STATES; ++state){
			sum_rho[state] = 0.0;
			sum_rho_sqr[state] = 0.0;
		}
	}
	/* end reseting rho & rho stats */
	/* begin generating initial ca configuration */
	if (num_states == 4){
		// calculating cumulated states' quantity for gkl s4 based on simulation parameters;
		s_count[0] = static_cast<int>(round(round(ca_len * z0) * x0));
		s_count[1] = static_cast<int>(round(round(ca_len * z0) * (1.0 - x0)) + s_count[0]);
		s_count[2] = static_cast<int>(round(round(ca_len * (1.0 - z0)) * x1) + s_count[1]);
		s_count[3] = static_cast<int>(round(round(ca_len * (1.0 - z0)) * (1.0 - x1)) + s_count[2]);
	}
	else if (num_states == 6){
		// calculating cumulated states' quantity for gkl s6 based on simulation parameters;
		s_count[0] = static_cast<int>(round(round(ca_len * z0) * x0));
		s_count[1] = static_cast<int>(round(round(ca_len * z0) * (1.0 - x0)) + s_count[0]);
		s_count[2] = static_cast<int>(round(round(round(ca_len * (1.0 - z0)) * x1) * z1) + s_count[1]);
		s_count[3] = static_cast<int>(round(round(round(ca_len * (1.0 - z0)) * (1.0 - x1)) * z1) + s_count[2]);
		s_count[4] = static_cast<int>(round(round(round(ca_len * (1.0 - z0)) * x2) * (1.0 - z1)) + s_count[3]);
		s_count[5] = static_cast<int>(round(round(round(ca_len * (1.0 - z0)) * (1.0 - x2)) * (1.0 - z1)) + s_count[4]);
	}
}

void __print_out_header_total(){
	if (!print_out_stats){
		if (!print_out_ca){
			cout << "# GKL initial configuration and final convergence:" << endl;
			// printing for gkl s4;
			if (num_states == 4){
				cout << "x0:      \tx1:      \tz:      \tconv.:\tsteps:" << endl;
			}
			// printing for gkl s6;
			else if (num_states == 6){
				cout << "x0:      \tx1:      \tx2:     \tz:      \tconv.:\tsteps:" << endl;
			}
		}
	}
	else {
		// printing rho stats header;
		cout << "# GKL Rho statistics:" << endl;
		cout << "s0 average:\ts0 std.dev.:\ts1 average:\ts1 std.dev.:\tconv.:";
		cout << endl;
	}
}
void __print_out_data_total(){
	if (!print_out_stats){
		if (print_out_ca){
			cout << endl;
			// printing for -ca option in gkl s4;
			if (num_states == 4){
				cout << "x0=" << x0 << "\tx1=" << x1 << "\tz=" << z0 << "\tconv.=" << convergence << endl;
			}
			// printing for -ca option in gkl s6;
			else if (num_states == 6){
				cout << "x0=" << x0 << "\tx1=" << x1 << "\tx2=" << x2 << "\tz=" << z0 << "\tconv.=" << convergence << endl;
			}
		}
		else {
			// printing for gkl s4;
			if (num_states == 4){
				cout << x0 << "\t" << x1 << "\t" << z0 << "\t" << convergence << "\t" << step << endl;
			}
			// printing for gkl s6;
			else if (num_states == 6){
				cout << x0 << "\t" << x1 << "\t" << x2 << "\t" << z0 << "\t" << convergence << "\t" << step << endl;
			}
		}
	}
	else {
		// printing rho stats;
		for (int state = 0; state < __GKL_MAIN_STATES; ++state){
			cout << sum_rho[state] / (step - stats_from + 1) / ca_len << "\t" << sqrt(sum_rho_sqr[state] / (step - stats_from + 1) - pow(sum_rho[state] / (step - stats_from + 1), 2)) / ca_len << "\t";
		}
		cout << convergence << "\t" << endl;
	}
}
void __print_out_header_step(){
	if (!print_out_stats){
		if (print_out_ca && print_out_rho){
			cout << endl;
			cout << "# GKL CA configuration dynamics with densities:" << endl;
			cout << "step:\tconfiguration (states' densities & convergence):" << endl;
		}
		else if (print_out_ca){
			cout << endl;
			cout << "# GKL CA configuration dynamics:" << endl;
			cout << "step:\tconfiguration (convergence):" << endl;
		}
		else if (print_out_rho){
			cout << endl;
			cout << "# GKL Rho (states' densities) dynamics:" << endl;
			if (num_states == 4) {
				cout << "step:\ts0:   \ts1:   \ts2:   \ts3:   \tconv.:" << endl;
			}
			else if (num_states == 6) {
				cout << "step:\ts0:   \ts1:   \ts2:   \ts3:   \ts4:   \ts5:   \tconv.:" << endl;
			}
		}
	}
}
void __print_out_data_step(){
	unsigned char __KISS_charmap[__GKL_MAX_STATES];
	if (print_ca_mapped){
		if (num_states == 4) {
			__KISS_charmap[0] = 0xb0;
			__KISS_charmap[1] = 0xdb;
			__KISS_charmap[2] = 0xb1;
			__KISS_charmap[3] = 0xb2;
		}
		else if (num_states == 6){
			__KISS_charmap[0] = 0xb0;
			__KISS_charmap[1] = 0xdb;
			__KISS_charmap[2] = 0xb1;
			__KISS_charmap[3] = 0xb2;
			__KISS_charmap[4] = 0x00;
			__KISS_charmap[5] = 0xFE;
		}
	}
	else {
		if (num_states == 4) {
			__KISS_charmap[0] = 'r';
			__KISS_charmap[1] = 'l';
			__KISS_charmap[2] = 'u';
			__KISS_charmap[3] = 'd';
		}
		else if (num_states == 6){
			__KISS_charmap[0] = '0';
			__KISS_charmap[1] = '1';
			__KISS_charmap[2] = '2';
			__KISS_charmap[3] = '3';
			__KISS_charmap[4] = '4';
			__KISS_charmap[5] = '5';
		}
	}
	if (!print_out_stats){
		if (print_out_ca && print_out_rho){
			// printing ca;
			cout << step << "\t";
			for (int cell = 1; cell < ca_len + 1; ++cell){ cout << __KISS_charmap[ca[step % 2][cell]]; }
			// printing rho;
			for (int state = 0; state < num_states; ++state){ cout << " " << __KISS_charmap[state] << "=" << rho[state] << ";"; }
			// printing convergence;
			cout << "\t" << convergence;
			cout << endl;
		}
		else if (print_out_ca){
			// printing ca;
			cout << step << "\t";
			for (int cell = 1; cell < ca_len + 1; ++cell){ cout << __KISS_charmap[ca[step % 2][cell]]; }
			// printing convergence;
			cout << "\t" << convergence;
			cout << endl;
		}
		else if (print_out_rho){
			// printing rho;
			cout << step;
			for (int state = 0; state < num_states; ++state){ cout << "\t" << rho[state]; }
			// printing convergence;
			cout << "\t" << convergence;
			cout << "\t" << entropy_before_noise;
			cout << "\t" << entropy_after_noise;
			cout << endl;
		}
	}
}

void __prt_simulation_header(){
	cout << "zeta_index";
	for (int state = 0; state < num_states; ++state){
		cout << '\t' << "s" << state;
	}
	cout << endl;
}
void __prt_simulation_footer(){
	for (int zeta_index = 0; zeta_index <= ca_len; ++zeta_index){
		cout << static_cast<double>(zeta_index) / ca_len;
		for (int state = 0; state < num_states; ++state){
			cout << '\t' << static_cast<double>(zeta[state][zeta_index]) / simulations / max_steps;
		}
		cout << endl;
	}
}
void __prt_steps_header(){

}
void __prt_steps_footer(){

}
void __prt_detail(){

}
/* end auxiliary functions */

/* begin gkl functions dynamics */
__GKL_CONVERGENCE verify_convergence(){
	double _tolerance = pow(noise, 2.0 / 3.0) * ca_len;
	if (__GKL_MAIN_STATES == 2){
		if ((rho[0] <= _tolerance && rho[1] >= ca_len - _tolerance) || (rho[1] <= _tolerance && rho[0] >= ca_len - _tolerance)) {
			if (((s_count[0] > s_count[1] - s_count[0]) && (rho[0] > rho[1])) || 
				((s_count[0] < s_count[1] - s_count[0]) && (rho[0] < rho[1]))){
				return __gkl_convergence_right;
			}
			else { return __gkl_convergence_wrong; }
		}
		else { return __gkl_convergence_none; }
	}
}
double __calc_entropy(){
	double entropy = 0.0;
	for (int state = 0; state < num_states; ++state){
		if (rho[state] > 0){
			entropy += (static_cast<double>(rho[state]) / ca_len) * log2(static_cast<double>(rho[state]) / ca_len) / log2(num_states);
		}
	}
	return -1.0 * entropy;
}
void calc_stats(){
	if (step >= stats_from){
		for (int state = 0; state < num_states; ++state){
			sum_rho[state] += rho[state];
			sum_rho_sqr[state] += pow(rho[state], 2);
		}
	}
}
void generate_new_ca(){
	step = 0;
	// reseting rho;
	for (int state = 0; state < num_states; ++state){ rho[state] = 0; }
	// generating linear array with states' quantity;
	for (int cell = 1; cell < ca_len + 1; ++cell){
		for (int state = 0; state < num_states; ++state){
			if (cell <= s_count[state]){ ca[step][cell] = state; ++rho[state];break; }
		}
	}
	for (int state = 0; state < num_states; ++state){ ++zeta[state][rho[state]]; }
	// shuffle array (Fisher-Yates shuffle);
	for (int cell = ca_len; cell > 1; --cell){
		int _swap_i = static_cast<int>(random(__gen) * (cell + 1));
		int cell_value = ca[step][cell];
		ca[step][cell] = ca[step][_swap_i];
		ca[step][_swap_i] = cell_value;
	}
	//entropy_before_noise = __calc_entropy();
	//entropy_after_noise = __calc_entropy();
	// periodic border condition;
	ca[step % 2][0] = ca[step % 2][ca_len];  ca[step % 2][ca_len + 1] = ca[step % 2][1];
	// printing detail;
	__prt_detail();
	// defining first convergence;
	convergence = __gkl_convergence_none;
}
/* end gkl functions dynamics */
/*
int rho_stats_count = 0;
int rho_stats_sum[__GKL_MAX_STATES] = { 0 };
int rho_stats_sum2[__GKL_MAX_STATES] = { 0 };
void calc_rho_stats(){
	for (int state = 0; state < num_states; ++state){
		rho_stats_sum[state] += rho[state];
		rho_stats_sum2[state] += pow(rho[state], 2);
		++rho_stats_count;
	}
}

void print_rho_stats(char separator = '\t', char end = '\n'){
	for (int state = 0; state < num_states; ++state){
		cout << sum_rho[state] / rho_stats_count / ca_len << separator << sqrt(rho_stats_sum2[state] / rho_stats_count - pow(rho_stats_sum[state] / rho_stats_count, 2)) / ca_len << (state + 1 == num_states ? end : separator);
	}
}
*/

// program entry point;
int main(int argc, char *argv[]){
	__init_main(argc, argv);
	// printing simulation header;
	__prt_simulation_header();
	/* begin iterating simulations */
	for (simulation = 0; simulation < simulations; ++simulation){
		// initializing simulation parameters;
		__init_simulation_parameters();
		// generate new ca configuration as step = 0;
		generate_new_ca();
		//printing step header;
		__prt_steps_header();
		/* begin iterating steps */
		for (step = 1; step < max_steps; ++step){
			// reseting rho;
			for (int state = 0; state < num_states; ++state){ rho[state] = 0; }
			/* begin applying psi */
			for (int cell = 1; cell < ca_len + 1; ++cell){
				ca[step % 2][cell] = psi[ca[(step + 1) % 2][cell - 1]][ca[(step + 1) % 2][cell]][ca[(step + 1) % 2][cell + 1]];
				++rho[ca[step % 2][cell]];
			}
			/* end applying psi */
			//entropy_before_noise = __calc_entropy();
			/* begin applying noise error */
			if (!stop_noise_on_step || (stop_noise_on_step && (step < stop_noise_on_step))){
				for (int cell = 1; cell < ca_len + 1; ++cell){
					if (noise && (random(__gen) < noise)){
						--rho[ca[step % 2][cell]];
						ca[step % 2][cell] = static_cast<int>(random(__gen) * num_states);
						++rho[ca[step % 2][cell]];
					}
				}
			}
			for (int state = 0; state < num_states; ++state){ ++zeta[state][rho[state]]; }
			/* end applying noise error */
			//entropy_after_noise = __calc_entropy();
			// periodic border condition;
			ca[step % 2][0] = ca[step % 2][ca_len];  ca[step % 2][ca_len + 1] = ca[step % 2][1];
			// verifying convergence;
			//convergence = verify_convergence();
			// printing detail;
			__prt_detail();
			// stop if ca converged if -c option informed
			if ((convergence != __gkl_convergence_none) && stop_on_convergence){ break; }
		}
		/* end iterating steps */
		// printing step footer;
		__prt_steps_footer();
		//// printing final results;
		//if (simulation == 0){ __print_out_header_total(); }
		//__print_out_data_total();
	}
	__prt_simulation_footer();
	/* end iterating simulations */
}