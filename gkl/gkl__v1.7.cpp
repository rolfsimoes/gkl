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
enum __GKL_CONVERGENCE
{
	__gkl_no_convergence = 0,
	__gkl_right_convergence = 1,
	__gkl_wrong_convergence = -1
};

#define __GKL_PAST_MEMORY 10
int _past_rho[__GKL_PAST_MEMORY][__GKL_MAIN_STATES];

/* begin global variables definitions */
bool stop_on_convergence = false;
bool x0_informed = false;
bool x1_informed = false;
bool x2_informed = false;
bool z_informed = false;
bool print_out_rho = false;
bool print_out_ca = false;
bool print_out_stats = false;
bool print_ca_mapped = true;
bool seed_informed = false;
bool max_steps_informed = false;
int num_states = 4;
int ca_len = 400;
int max_steps = 4000;
int stats_from = 0;
int simulations = 1;
int step = 0;
unsigned int seed = 0;
__GKL_CONVERGENCE convergence = __gkl_no_convergence;
double x0 = 0.0;
double x1 = 0.0;
double x2 = 0.0;
double z = 0.0;
double noise = 0.0;
int ca[2][__GKL_MAX_LEN + 2] = { 0 };
int psi[__GKL_MAX_STATES][__GKL_MAX_STATES][__GKL_MAX_STATES] = { 0 };
int rho[__GKL_MAX_STATES] = { 0 };
int s_count[__GKL_MAX_STATES] = { 0 };
double sum_rho[__GKL_MAIN_STATES] = { 0 };
double sum_rho_sqr[__GKL_MAIN_STATES] = { 0 };
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
unsigned int __x = 123456789, __y = 987654321, __z = 43219876, __c = 6543217;
// initialyzing KISS generator seed;
void __init_rand_gen(){
	if (seed_informed){
		__x = seed;
	}
	else {
		__x = __rd();
		while (!(__y = __rd())); /* y must not be zero! */
		__z = __rd();
		__c = __rd() % 698769068 + 1; /* Should be less than 698769069 */
	}
}
// main function KISS generator proposed by G. Marsaglia (see source above to references);
unsigned int __gen(){
	unsigned long long t;
	__x = 314527869 * __x + 1234567;
	__y ^= __y << 5; __y ^= __y >> 7; __y ^= __y << 22;
	t = 4294584393ULL * __z + __c;
	__c = static_cast<unsigned int>(t >> 32);
	__z = static_cast<unsigned int>(t);
	return __x + __y + __z;
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
	cout << "# Command: " << argv[0] << " " << argv[1] << " ";
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
		else if ((strcmp(argv[i], "-z") == 0) && (argc - 1 > i)){
			cout << argv[i] << " " << argv[i + 1] << " ";
			z = atof(argv[++i]);
			z_informed = true;
			if (z < 0 || z > 1){
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
		else if ((strcmp(argv[i], "-stats") == 0) && (argc - 1 > i)){
			cout << argv[i] << " " << argv[i + 1] << " ";
			print_out_stats = true;
			stats_from = atol(argv[++i]);
			if ((max_steps < stats_from) && max_steps_informed){
				__usage(argv[0]);
				exit(1);
			}
		}
		else if ((strcmp(argv[i], "-seed") == 0) && (argc - 1 > i)){
			cout << argv[i] << " " << argv[i + 1] << " ";
			seed_informed = true;
			seed = atol(argv[++i]);
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
				cout << "x0=" << x0 << "\tx1=" << x1 << "\tz=" << z << "\tconv.=" << convergence << endl;
			}
			// printing for -ca option in gkl s6;
			else if (num_states == 6){
				cout << "x0=" << x0 << "\tx1=" << x1 << "\tx2=" << x2 << "\tz=" << z << "\tconv.=" << convergence << endl;
			}
		}
		else {
			// printing for gkl s4;
			if (num_states == 4){
				cout << x0 << "\t" << x1 << "\t" << z << "\t" << convergence << "\t" << step << endl;
			}
			// printing for gkl s6;
			else if (num_states == 6){
				cout << x0 << "\t" << x1 << "\t" << x2 << "\t" << z << "\t" << convergence << "\t" << step << endl;
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
	unsigned char __charmap[__GKL_MAX_STATES];
	if (print_ca_mapped){
		if (num_states == 4) {
			__charmap[0] = 0xb0;
			__charmap[1] = 0xdb;
			__charmap[2] = 0xb1;
			__charmap[3] = 0xb2;
		}
		else if (num_states == 6){
			__charmap[0] = 0xb0;
			__charmap[1] = 0xdb;
			__charmap[2] = 0xb1;
			__charmap[3] = 0xb2;
			__charmap[4] = 0x00;
			__charmap[5] = 0xFE;
		}
	}
	else {
		if (num_states == 4) {
			__charmap[0] = 'r';
			__charmap[1] = 'l';
			__charmap[2] = 'u';
			__charmap[3] = 'd';
		}
		else if (num_states == 6){
			__charmap[0] = '0';
			__charmap[1] = '1';
			__charmap[2] = '2';
			__charmap[3] = '3';
			__charmap[4] = '4';
			__charmap[5] = '5';
		}
	}
	if (!print_out_stats){
		if (print_out_ca && print_out_rho){
			// printing ca;
			cout << step << "\t";
			for (int cell = 1; cell < ca_len + 1; ++cell){ cout << __charmap[ca[step % 2][cell]]; }
			// printing rho;
			for (int state = 0; state < num_states; ++state){ cout << " " << __charmap[state] << "=" << rho[state] << ";"; }
			// printing convergence;
			cout << "\t" << convergence;
			cout << endl;
		}
		else if (print_out_ca){
			// printing ca;
			cout << step << "\t";
			for (int cell = 1; cell < ca_len + 1; ++cell){ cout << __charmap[ca[step % 2][cell]]; }
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
			cout << endl;
		}
	}
}
/* end auxiliary functions */

__GKL_CONVERGENCE verify_convergence(){
	double _tolerance = pow(noise, 2.0 / 3.0) * ca_len;
	if (__GKL_MAIN_STATES == 2){
		if ((rho[0] <= _tolerance && rho[1] >= ca_len - _tolerance) || (rho[1] <= _tolerance && rho[0] >= ca_len - _tolerance)) {
			if (((s_count[0] > s_count[1] - s_count[0]) && (rho[0] > rho[1])) || 
				((s_count[0] < s_count[1] - s_count[0]) && (rho[0] < rho[1]))){
				return __gkl_right_convergence;
			}
			else { return __gkl_wrong_convergence; }
		}
		else { return __gkl_no_convergence; }
	}
}
__GKL_CONVERGENCE verify_convergence_eduardo(){
	for (int state = 0; state < __GKL_MAIN_STATES; ++state){
		double _verifica_converge = 1.0;
		_past_rho[step % __GKL_PAST_MEMORY][state] = rho[state];
		if (step > 500){
			for (int delta_step = 0; delta_step < __GKL_PAST_MEMORY; ++delta_step){
				_verifica_converge *= 1.0 + (static_cast<double>(_past_rho[(step - delta_step) % __GKL_PAST_MEMORY][state]) -
					static_cast<double>(_past_rho[(step - delta_step - 1) % __GKL_PAST_MEMORY][state])) /
					static_cast<double>(_past_rho[(step - delta_step - 1) % __GKL_PAST_MEMORY][state]);
			}
		}
		if (_verifica_converge < 1.1 && _verifica_converge > 0.9 && rho[state] > ca_len / 2 && step > 500){
			if ((x0 >= 0.5 && rho[0] > rho[1]) || (x0 <= 0.5 && rho[0] < rho[1])){ return __gkl_right_convergence; }
			else { return __gkl_wrong_convergence; }
		}
	}
	return __gkl_no_convergence;
}
void calc_stats(){
	if (step >= stats_from){
		for (int state = 0; state < __GKL_MAIN_STATES; ++state){
			sum_rho[state] += rho[state];
			sum_rho_sqr[state] += pow(rho[state], 2);
		}
	}
}

/* begin gkl functions dynamics */
void run_ca(){
	do {
		++step;
		// reseting rho;
		for (int state= 0; state < num_states; ++state){ rho[state] = 0; }
		/* begin applying psi */
		for (int cell = 1; cell < ca_len + 1; ++cell){
			ca[step % 2][cell] = psi[ca[(step + 1) % 2][cell - 1]][ca[(step + 1) % 2][cell]][ca[(step + 1) % 2][cell + 1]];
			/* begin applying noise error */
			if ((noise > 0) && (random(__gen) < noise)){ ca[step % 2][cell] = static_cast<int>(random(__gen) * num_states); }
			/* end applying noise error */
			++rho[ca[step % 2][cell]];
		}
		/* end applying psi */
		// periodic border condition;
		ca[step % 2][0] = ca[step % 2][ca_len];
		ca[step % 2][ca_len + 1] = ca[step % 2][1];
		// verifying convergence;
		convergence = verify_convergence();
		// calc rho stats;
		if (print_out_stats){ calc_stats(); }
		// printing out data;
		__print_out_data_step();
	} while (step < max_steps && !((convergence != __gkl_no_convergence) && stop_on_convergence && !print_out_stats));
}
// program entry point;
int main(int argc, char *argv[]){
	// config output stream;
	cout << setiosflags(ios::fixed) << setprecision(6);
    // handling command-line options;
	__handle_options(argc, argv);
	// initilyzing random function (after __handle_options);
	__init_rand_gen();
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
	/* begin iterating simulations */
	for (int simulation = 0; simulation < simulations; ++simulation){
		/* begin initializing simulation parameters */
		step = 0;
		convergence = __gkl_no_convergence;
		if (!x0_informed) { x0 = random(__gen); }
		if (!x1_informed) { x1 = random(__gen); }
		if (!x2_informed) { x2 = random(__gen); }
		if (!z_informed) { z = random(__gen); }
		/* end initializing simulation parameters */
		/* begin reseting rho & rho stats */
		for (int state = 0; state < num_states; ++state){ rho[state] = 0; }
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
			s_count[0] = static_cast<int>(round(round(ca_len * z) * x0));
			s_count[1] = static_cast<int>(round(round(ca_len * z) * (1.0 - x0)) + s_count[0]);
			s_count[2] = static_cast<int>(round(round(ca_len * (1.0 - z)) * x1) + s_count[1]);
			s_count[3] = static_cast<int>(round(round(ca_len * (1.0 - z)) * (1.0 - x1)) + s_count[2]);
		}
		else if (num_states == 6){
			// calculating cumulated states' quantity for gkl s6 based on simulation parameters;
			s_count[0] = static_cast<int>(round(round(ca_len * z) * x0));
			s_count[1] = static_cast<int>(round(round(ca_len * z) * (1.0 - x0)) + s_count[0]);
			s_count[2] = static_cast<int>(round(round(ca_len * (1.0 - z) / 2) * x1) + s_count[1]);
			s_count[3] = static_cast<int>(round(round(ca_len * (1.0 - z) / 2) * (1.0 - x1)) + s_count[2]);
			s_count[4] = static_cast<int>(round(round(ca_len * (1.0 - z) / 2) * x2) + s_count[3]);
			s_count[5] = static_cast<int>(round(round(ca_len * (1.0 - z) / 2) * (1.0 - x2)) + s_count[4]);
		}
		// generating linear array with states' quantity;
		for (int cell = 1; cell < ca_len + 1; ++cell){
			for (int state = 0; state < num_states; ++state){
				if (cell < s_count[state]){ ca[step][cell] = state; ++rho[state]; break; }
			}
		}
		// shuffle array (Fisher-Yates shuffle);
		for (int cell = ca_len; cell > 1; --cell){
			int _swap_i = static_cast<int>(random(__gen) * (cell + 1));
			int cell_value = ca[step][cell];
			ca[step][cell] = ca[step][_swap_i];
			ca[step][_swap_i] = cell_value;
		}
		/* end generating initial ca configuration */
		// periodic border condition;
		ca[step % 2][0] = ca[step % 2][ca_len];
		ca[step % 2][ca_len + 1] = ca[step % 2][1];
		// printing out data for step = 0;
		__print_out_header_step();
		__print_out_data_step();
		// calc rho stats for step = 0 if -stats option is informed;
		if (print_out_stats){ calc_stats(); }
		// run ca for steps > 0 until max steps or convergence is reached if -c is informed;
		run_ca();
		// printing final results;
		if (simulation == 0){ __print_out_header_total(); }
		__print_out_data_total();
	}
	/* end iterating simulations */
}
/* end gkl functions dynamics */