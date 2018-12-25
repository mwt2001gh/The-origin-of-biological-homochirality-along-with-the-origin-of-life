// "The origin of biological homochirality along with the origin of life"
//by Yong Chen and Wentao Ma*.
//C source codes for the simulation program
//The case is correponding to Fig.3a in the article

#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cstring>

// Random number generator
#define A1 471
#define B1 1586
#define C1 6988
#define D1 9689
#define M 16383
#define RIMAX 2147483648.0        // = 2^31 
#define RandomInteger (++nd, ra[nd & M] = ra[(nd-A1) & M] ^ ra[(nd-B1) & M] ^ ra[(nd-C1) & M] ^ ra[(nd-D1) & M])
#define RANDD ((double)RandomInteger / RIMAX)

//Nucleotides
#define C 2
#define G 3
#define A 1
#define U 4

#define NSRSEQ         A,C,U,G,G,C // The characteristic (catalytic domain) sequence of a nucleotide synthetase ribozyme (NSR)
#define NSRCOMSEQ      G,C,C,A,G,U

#define REPSEQ	      G,U,U,C,A,G  // The characteristic (catalytic domain) sequence of an RNA replicase (REP)
#define REPCOMSEQ     C,U,G,A,A,C

#define RECTXT "data.txt"   //Text file for output
#define STEPNUM 20000000    // Total time steps of Monte Carlo simulation
#define STAREC 0            // The step to start record
#define RECINT 10000        // The interval steps of recording
#define MAX_RNA_LENGTH 100   // Maximum RNA length allowed in the simulation
#define MAXRNASINGRID 10000  // Maximum number of RNA allowed in a single grid room
#define LEN sizeof(struct rna)
#define ROOMNUM (N*N)       // The room number in the grid
#define N 20                // The system surface is defined as an N¡ÁN grid
#define RMRW  (pow(p->length1+p->length2,1/2.0))  	// The relationship between the movement of RNAs and their molecular weight

struct rna								// Nucleotide or RNA
{
	char information[2][MAX_RNA_LENGTH];
	char chirality[2][MAX_RNA_LENGTH];
	int length1;
	int length2;
	int nick[MAX_RNA_LENGTH];
	struct rna *next;
	struct rna *prior;
};

int    SD = 5;          // The random seed
/************* Parameters *************/
int    TNPB = 50000;    // Total nucleotide precursors (quotients in measurement of nucleotides) introduced in the beginning
double FCSS = 0.5;      // The factor for chiral selection in surface-mediated synthesis
double FCST = 0.5;      // Factor of chirality selectivity in the template-directed replication
double PAT = 0.5;	    // Probability of an RNA template attracting a substrate (nucleotide or oligomer)
double PBB = 0.00001;   // Probability of a phosphodiester bond breaking within an RNA chain
double PCIC = 0.5;      // Probability of the chirality inter-conversion of nucleotide precursors 
double PFP = 0.001;     // Probability of the false base-pairing (relating to RNA's replicating fidelity)
double PMN = 0.0001;    // Probability of the movement of nucleotides
double PMPN = 0.002;    // Probability of the movement of nucleotide precursors
double PND = 0.01;      // Probability of a nucleotide decaying into its precursor
double PNDE = 0.0001;   // Probability of a nucleotide residue decaying at RNA¡¯s chain end
double PNF = 0.001;     // Probability of a nucleotide forming from its precursor (non-enzymatic)
double PNFR = 0.2;	    // Probability ofA nucleotide forming from its precursor catalyzed by NSR
double PRL = 0.000002;  // Probability of the random ligation of nucleotides and RNAs (surface-mediated) 
double PSP = 0.3;	    // Probability of the separation of a base pair
double PTL = 0.002;	    // Probability of the template-directed ligation (non-enzymatic)
double PTLR = 0.0;      // Probability of the template-directed ligation catalyzed by Rep
int    CT = 8;          // Collision times of nucleotides and RNAs within a grid room in a Monte-Carlo step

static long ra[M + 1], nd;  // For the random number generator
struct rna *room_head[2][N][N]; // Head nodes for chain-tables to record nucleotides and RNAs in grid rooms
struct rna *p, *p1, *p2, *p3;
struct rna *rna_arr[MAXRNASINGRID];

char nsrseq[50] = { NSRSEQ };
char nsrcomseq[50] = { NSRCOMSEQ };
char repseq[50] = { REPSEQ };
char repcomseq[50] = { REPCOMSEQ };

int lpn_arr[2][N][N];	// Precursors of L-nucleotides in grid rooms
int dpn_arr[2][N][N];	// Precursors of D-nucleotides in grid rooms

int x, y;				// The coordinate of rooms in the grid 
long i;					// Cycle variable for Monte Carlo steps
int g = 0;				// Recording times
int over_max_len = 0;	 

int k, pn_bef, randcase;
int replength, repcomlength, nsrcomlength, nsrlength;
int flag, flag1, flagntsyn, flagtdlig;

// For data-recording
int nsr_num[(STEPNUM - STAREC) / RECINT + 1];
int d_nsr_num[(STEPNUM - STAREC) / RECINT + 1];
int l_nsr_num[(STEPNUM - STAREC) / RECINT + 1];

int rep_num[(STEPNUM - STAREC) / RECINT + 1];
int d_rep_num[(STEPNUM - STAREC) / RECINT + 1];
int l_rep_num[(STEPNUM - STAREC) / RECINT + 1];

int total_nt_mat[(STEPNUM - STAREC) / RECINT + 1];				
int rna_num[(STEPNUM - STAREC) / RECINT + 1];					
int pn_num[(STEPNUM - STAREC) / RECINT + 1];					

int chnum_d[(STEPNUM - STAREC) / RECINT + 1][MAX_RNA_LENGTH]; 
int	chnum_l[(STEPNUM - STAREC) / RECINT + 1][MAX_RNA_LENGTH];

int total_d;
int total_l;

/////////////////////////////////////////////////////////////////////////
void seed(long seed)   // Random number initialization    
{
	int a;
	if (seed<0) { printf("SEED error."); exit(0); }
	ra[0] = (long)fmod(16807.0*(double)seed, 2147483647.0);
	for (a = 1; a <= M; a++)
	{
		ra[a] = (long)fmod(16807.0 * (double)ra[a - 1], 2147483647.0);
	}
}


/////////////////////////////////////////////////////////////////////////
void fresh_rna(int h)     // Updating a nucleotide or RNA for the next time step
{
	p1 = p->prior;
	p2 = p->next;
	p3 = room_head[!h][y][x]->next;
	room_head[!h][y][x]->next = p;
	p->next = p3;
	p->prior = room_head[!h][y][x];
	if (p3 != room_head[!h][y][x])p3->prior = p;
	p1->next = p2;
	if (p2 != room_head[h][y][x])p2->prior = p1;
	p = p1;
}

/////////////////////////////////////////////////////////////////////////
void rna_shuffle(int h) //Random order-arrangement of nodes in the chain-table for recording nucleotides and RNA in a grid room 
{
	int a, b, len;
	for (a = 0, p = room_head[h][y][x]->next; p != room_head[h][y][x]; p = p->next)
	{
		if (a == MAXRNASINGRID) {printf("Too many RNA molecules in a single grid room"); exit(0);}
		rna_arr[a] = p;
		a++;
	}
	len = a;
	rna_arr[a] = NULL;
	for (a = len - 1; a >= 1; a--)
	{
		p = rna_arr[a];
		b = RandomInteger % (a);
		rna_arr[a] = rna_arr[b];
		rna_arr[b] = p;
	}
	p = room_head[h][y][x];
	for (a = 0; a < len; a++)
	{
		p->next = rna_arr[a];
		rna_arr[a]->prior = p;
		p = p->next;
	}
	p->next = room_head[h][y][x];

}

/////////////////////////////////////////////////////////////////////////
int last_nick(struct rna *p)    //Finding the last unlinked site on the synthesizing chain on the template
{
	for (int z = p->length2 - 1; z>0; z--)
	{
		if (p->nick[z] != 0) return z;
	}
	return 0; // no nick
}

/////////////////////////////////////////////////////////////////////////
void save_parameters(const char * FILENAME) // Recording the parameters used in the case
{
	FILE *fp;

	printf("----------------------------------------------------------------------------------------------------------------------------------------\n");
	printf("||   SD:%-7d         N:%-7d      TNPB:%-7d      FCST:%-7f      FCSS:%-7f  ||\n", SD, N, TNPB, FCST, FCSS);
	printf("||  PAT:%-7f       PFP:%-7f       PSP:%-7f       PND:%-7f      PNDE:%-7f  ||\n", PAT, PFP, PSP, PND, PNDE);
	printf("||  PNF:%-7f      PNFR:%-7f       PRL:%-7f       PCIC:%-7f       PBB:%-7f  ||\n", PNF, PNFR, PRL, PCIC, PBB);
	printf("||  PTL:%-7f      PTLR:%-7f       PMN:%-7f      PMPN:%-7f                 ||\n", PTL, PTLR, PMN, PMPN);
	printf("-------------------------------------------------------------------------------------------------------------------------------------------------------\n");

	if ((fp = fopen(FILENAME, "at")) == NULL) { printf("cannot open file"); exit(-1); }
	fprintf(fp, "----------------------------------------------------------------------------------------------------------------------------------------\n");
	fprintf(fp, "||   SD:%-7d         N:%-7d      TNPB:%-7d      FCST:%-7f      FCSS:%-7f  ||\n", SD, N, TNPB, FCST, FCSS);
	fprintf(fp, "||  PAT:%-7f       PFP:%-7f       PSP:%-7f       PND:%-7f      PNDE:%-7f  ||\n", PAT, PFP, PSP, PND, PNDE);
	fprintf(fp, "||  PNF:%-7f      PNFR:%-7f       PRL:%-7f       PCIC:%-7f       PBB:%-7f  ||\n", PNF, PNFR, PRL, PCIC, PBB);
	fprintf(fp, "||  PTL:%-7f      PTLR:%-7f       PMN:%-7f      PMPN:%-7f                 ||\n", PTL, PTLR, PMN, PMPN);
	fprintf(fp, "-------------------------------------------------------------------------------------------------------------------------------------------------------\n");
	fclose(fp);
}

int is_nsr_chirality(rna *p)  // Return the chirality type if find NSR's characteristic sequence,else return 0
{
	int flag2, a, b;

	if (p->length1 >= nsrlength)
	{
		for (b = 0; p->length1 - nsrlength - b >= 0; b++)
		{
			flag2 = 1;
			for (a = 0; a < nsrlength; a++)
			{
				if (p->information[0][b + a] == nsrseq[a] && p->chirality[0][b + a] == p->chirality[0][b]);
				else { flag2 = 0; break; }
			}
			if (flag2) return p->chirality[0][b];//found
		}

	}
	return 0;//Not found 
}

int is_rep_chirality(rna *p)  // Return the chirality type if find REP's characteristic sequence,else return 0
{
	int flag2, a, b;

	if (p->length1 >= replength)
	{
		for (b = 0; p->length1 - replength - b >= 0; b++)
		{
			flag2 = 1;
			for (a = 0; a < replength; a++)
			{
				if (p->information[0][b + a] == repseq[a] && p->chirality[0][b + a] == p->chirality[0][b]);
				else { flag2 = 0; break; }
			}
			if (flag2) return p->chirality[0][b];//found
		}

	}
	return 0;//Not found 
}

/////////////////////////////////////////////////////////////////////////

bool check_rna_chirality(rna *p, int *chirality, int *ulength)   // For recording chain_length of D or L type RNAs
{
	if (p->chirality[0][0] == p->chirality[0][p->length1 - 1])  // For nucleotides or handedness-uniform RNAs 
	{
		*chirality = p->chirality[0][0];     // its chirality is recorded
		*ulength = p->length1;               // its length is recorded as 'ulength'
		return true;
	}
	else // For handedness-mixed RNAs, the chirality of the longer domain is recorded, and the length of the longer domain is recorded as 'ulength'. 
	{
		int l_length = 0, r_length = 0;
		for (int i = 0; i < p->length1; i++) {
			if (p->chirality[0][i] != p->chirality[0][0])
				break;
			else
				l_length++;
		}
		r_length = p->length1 - l_length;
		if (l_length > r_length) { *chirality = p->chirality[0][0]; *ulength = l_length; }
		else { *chirality = p->chirality[0][p->length1 - 1]; *ulength = r_length; }

		if (*ulength > 1) return true;
		else return false; // For ulength=1, which means handedness-mixed dimer (i.e., D-L dimer is not counted)
	}
}

/////////////////////////////////////////////////////////////////////////
void inits(void)         // Initialization of the system
{
	int j, m, k;
	nd = 0;
	seed(SD);

	nsrlength = 0;         // Caculate the length of the active domain of the Nsr
	for (j = 0; nsrseq[j] != 0; j++)
		nsrlength++;

	replength = 0;         // Caculate the length of the active domain of Rep
	for (j = 0; repseq[j] != 0; j++)
		replength++;

	for (m = 0; m<2; m++)  // Initiate the data chain array of nucleotides and RNA
	{
		for (y = 0; y<N; y++)
		{
			for (x = 0; x<N; x++)
			{
				p1 = (struct rna *)malloc(LEN);
				if (!p1) { printf("\tinit1--memeout\n"); exit(-1); }
				room_head[m][y][x] = p1;
				p1->next = room_head[m][y][x];
			}
		}
	}

	for (k = 0; k<TNPB / 2; k++)  // Initial distribution of L-nucleotide precursors
	{
		x = RandomInteger % (N);
		y = RandomInteger % (N);
		lpn_arr[0][y][x]++;
	}
	for (k = 0; k<TNPB / 2; k++)  // Initial distribution of D-nucleotide precursors
	{
		x = RandomInteger % (N);
		y = RandomInteger % (N);
		dpn_arr[0][y][x]++;
	}
}

/////////////////////////////////////////////////////////////////////////
void unit_action(void)      // Action (movement and events) of units (molecules) in the system
{
	int a, b, j, randnt;
	double f, f1;

	//=================================== Movement of molecules

	//----------------- Nucleotides and RNA moving
	for (y = 0; y<N; y++)
	{
		for (x = 0; x<N; x++)
		{
			for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
			{
				if (RANDD*RMRW<PMN)
				{
					randcase = RandomInteger % (4);   // Four possible directions
					switch (randcase)
					{
						case 0:
							p1 = p->prior;
							p2 = p->next;
							p3 = room_head[1][y][(N + x - 1) % N]->next;
							room_head[1][y][(N + x - 1) % N]->next = p;
							p->next = p3;
							p->prior = room_head[1][y][(N + x - 1) % N];
							if (p3 != room_head[1][y][(N + x - 1) % N])p3->prior = p;
							p1->next = p2;
							if (p2 != room_head[0][y][x])p2->prior = p1;
							p = p1;
							break;
						case 1:
							p1 = p->prior;
							p2 = p->next;
							p3 = room_head[1][y][(x + 1) % N]->next;
							room_head[1][y][(x + 1) % N]->next = p;
							p->next = p3;
							p->prior = room_head[1][y][(x + 1) % N];
							if (p3 != room_head[1][y][(x + 1) % N])p3->prior = p;
							p1->next = p2;
							if (p2 != room_head[0][y][x])p2->prior = p1;
							p = p1;
							break;
						case 2:
							p1 = p->prior;
							p2 = p->next;
							p3 = room_head[1][(N + y - 1) % N][x]->next;
							room_head[1][(N + y - 1) % N][x]->next = p;
							p->next = p3;
							p->prior = room_head[1][(N + y - 1) % N][x];
							if (p3 != room_head[1][(N + y - 1) % N][x])p3->prior = p;
							p1->next = p2;
							if (p2 != room_head[0][y][x])p2->prior = p1;
							p = p1;
							break;
						case 3:
							p1 = p->prior;
							p2 = p->next;
							p3 = room_head[1][(y + 1) % N][x]->next;
							room_head[1][(y + 1) % N][x]->next = p;
							p->next = p3;
							p->prior = room_head[1][(y + 1) % N][x];
							if (p3 != room_head[1][(y + 1) % N][x])p3->prior = p;
							p1->next = p2;
							if (p2 != room_head[0][y][x])p2->prior = p1;
							p = p1;
							break;
						default: printf("RNA moving error");
					}
				}
				else fresh_rna(0);
			}
		}
	}

	//----------------- Nucleotide precursors moving
	for (y = 0; y<N; y++)	// L-nucleotide precursors
	{
		for (x = 0; x<N; x++)
		{
			pn_bef = lpn_arr[0][y][x];
			lpn_arr[0][y][x] = 0;
			for (k = 0; k<pn_bef; k++)
			{
				if (RANDD<PMPN)   // Moving
				{
					randcase = RandomInteger % (4);   // Four possible directions
					switch (randcase)
					{
						case 0:
							lpn_arr[1][y][(N + x - 1) % N]++; 
							break;
						case 1:
							lpn_arr[1][y][(x + 1) % N]++;
							break;
						case 2:
							lpn_arr[1][(N + y - 1) % N][x]++;
							break;
						case 3:
							lpn_arr[1][(y + 1) % N][x]++;
							break;
						default: printf("pn moving error");
					}
				}
				else lpn_arr[1][y][x]++; // No moving
			}
		}
	}
	
	for (y = 0; y<N; y++)      // D-nucleotide precursors
	{
		for (x = 0; x<N; x++)
		{
			pn_bef = dpn_arr[0][y][x];
			dpn_arr[0][y][x] = 0;
			for (k = 0; k<pn_bef; k++)
			{
				if (RANDD<PMPN)   // Moving
				{
					randcase = RandomInteger % (4);   // Four possible directions
					switch (randcase)
					{
						case 0:
							dpn_arr[1][y][(N + x - 1) % N]++; 
							break;
						case 1:
							dpn_arr[1][y][(x + 1) % N]++;
							break;
						case 2:
							dpn_arr[1][(N + y - 1) % N][x]++;
							break;
						case 3:
							dpn_arr[1][(y + 1) % N][x]++;
							break;
						default: printf("pn moving error");
					}
				}
				else dpn_arr[1][y][x]++; // No moving
			}
		}
	}

	//========== Events of molecules
	for (y = 0; y<N; y++)
	{
		for (x = 0; x<N; x++)
		{
			//===================Nucleotides and RNA's events

			// crash test
			for (int collision = 0; collision < CT; collision++)
			{
				rna_shuffle(1); // Random order-arrangement of nodes in the chain-table for recording nucleotides and RNA in a grid room 
				for (p = room_head[1][y][x]->next; p != room_head[1][y][x];)
				{
					if (p->next == room_head[1][y][x])
						break;
					p3 = p->next;

					switch (RandomInteger % 2)
					{
						case 0: // An RNA template attracting a substrate (note the three points below associated with the consideration of chirality)
							// 1.Cross inhibition: if an handedness-opposite residue is at the end of the synthesizing strand, the elongation will be terminated.
							if (p->length2>0 && p->chirality[0][p->length2 - 1] != p->chirality[1][p->length2 - 1]) break;

							if (p3->length2 != 0) break;
							if (p3->length1 > p->length1 - p->length2) break; // p3 is too long to combine with p

							flag = 0;  
							for (int b = 0; b < p3->length1; b++)
							{
								// 2. The handedness-opposite residue in the template strand would block the attraction.
								if (p->chirality[0][p->length2 + b] != p->chirality[0][0]) { flag = 1; break; }
								// 3. If the sustrate is not a mononucleotide, the handedness of all its residues must be the same with the template.
								else if (p3->length1 != 1 && p3->chirality[0][p3->length1 - 1 - b] != p->chirality[0][p->length2 + b]) { flag = 1; break; }

								// base-pairing check
								if ((p3->information[0][p3->length1 - 1 - b] + p->information[0][p->length2 + b]) == 5)	continue; // correct base-pairing
								else if (RANDD < PFP)	continue; // tolerable mispairing 
								else {flag = 1;	break;} // intolerable mispairing
							}

							if (flag == 0)
							{
								if (RANDD < PAT)
								{
									for (int a = 0; a < p3->length1; a++)
									{
										p->information[1][p->length2 + a] = p3->information[0][p3->length1 - 1 - a];
										p->chirality[1][p->length2 + a] = p3->chirality[0][p3->length1 - 1 - a];
									}
									p->information[1][p->length2 + p3->length1] = 0;
									p->nick[p->length2] = 1;
									for (int a = 1; a < p3->length1; a++) p->nick[p->length2 + a] = 0;

									p->length2 = p->length2 + p3->length1;
									(p3->prior)->next = p3->next;
									if (p3->next != room_head[1][y][x])(p3->next)->prior = p3->prior;
									free(p3);
									p3 = p;
								}
							}
							break;
						case 1: // Surface-mediated synthesis
							if (p3->length2 == 0)
							{
								// About the primer effect in the surface-mediated synthesis
								if (p->length1 == 1) f = 1;
								else if (p->length1 == 2) f = 10;
								else f = 20;
								//For the cases that do not consider the primer effect of the surface-mediated synthesis, f always equals to 1.

								if (p->chirality[0][0] != p->chirality[0][p->length1 - 1]) break;// uniform strand (only D or only L)
								if (p3->chirality[0][0] != p3->chirality[0][p3->length1 - 1]) break;// uniform strand (only D or only L)
								if (p3->chirality[0][0] != p->chirality[0][p->length1 - 1])
								{
									if (p3->length1 != 1) break; // Only mononucleotide can be incorporated when the handedness of the two molecules is different
									else f *= FCSS; //  In practice, FCSS<1, that is, the incorporation prefers the same handedness sustrates
								}

								if (RANDD < f*PRL / p3->length1)
								{
									if (p->length1 + p3->length1 > MAX_RNA_LENGTH - 1)
									{
										over_max_len++;
										break;
									}
									for (int a = 0; a < p3->length1; a++)
									{
										p->information[0][a + p->length1] = p3->information[0][a];
										p->chirality[0][a + p->length1] = p3->chirality[0][a];
									}
									p->length1 = p->length1 + p3->length1;
									p->information[0][p->length1] = 0;
									(p3->prior)->next = p3->next;
									if (p3->next != room_head[1][y][x]) (p3->next)->prior = p3->prior;
									free(p3);
									p3 = p; // p3 point to p
								}
							}
							break;
						default: printf("collision error");
					}
					p = p3->next;
				}
			}

			int l_rep_turn = 0;
			int d_rep_turn = 0;

			// Finding REP molecules in the grid room
			for (p = room_head[1][y][x]->next; p != room_head[1][y][x]; p = p->next)
			{	
				
				if (p->length1 > 2 * replength || p->length2 != 0 ) continue;  
								// Only a single chain containing the REP domain but no longer than twice its length may act as a REP
				flag1 = is_rep_chirality(p); // The value of flag1 represents the chirality type of the REP.
				if (flag1 != 0) {
					if (flag1 == -1) l_rep_turn++;  // Counting possible catalysis-turns in the room 
					else d_rep_turn++;
				}
			}

			// Template-directed ligating
			for (p = room_head[1][y][x]->next; p != room_head[1][y][x]; p = p->next)
			{
				for (a = 1; a < p->length2; a++)
				{
					if (p->nick[a] == 0)
						continue;
					//About the primer effect in template-directed synthesis
					if (a == 1)f = 1;
					else if (a == 2)f = 5;
					else f = 10;
					//For the cases that do not consider the primer effect of the template-directed synthesis, f always equals to 1.

					flagtdlig = 1;
					if (p->chirality[1][a - 1] == p->chirality[1][a])
					{
						if (p->chirality[1][a] == 1)
						{
							if (d_rep_turn>0)
							{
								d_rep_turn--;
								if (RANDD<PTLR*f) flagtdlig = 0; // Ligation catalyzed by REP
							}
						}
						else if (p->chirality[1][a] == -1)
						{
							if (l_rep_turn>0)
							{
								l_rep_turn--;
								if (RANDD<PTLR*f) flagtdlig = 0; // Ligation catalyzed by REP
							}
						}
						else {printf("unexpected error on the chirality of the substrate"); exit(0);}
					}

					if (p->chirality[1][a - 1] != p->chirality[1][a]) f *= FCST; // Template-directed ligation prefers substrates with the same handedness
					if (RANDD < PTL*f) flagtdlig = 0; // Uncatalyzed ligation

					if (flagtdlig == 0) p->nick[a] = 0;
				}
			}

			// Double-chain separating
			for (p = room_head[1][y][x]->next; p != room_head[1][y][x]; p = p->next)
			{
				if (p->length2 != 0)
				{
					a = last_nick(p);

					if (RANDD< pow(PSP, (p->length2 - a + 1) / 2.0)) //With consideration on the influence of single chain folding, which favors the double-chain separation    
					{
						p3 = (struct rna *)malloc(LEN);
						if (!p3) { printf("\t%ldsep--memeout\n", i); exit(-1); }

						for (b = 0; b < p->length2 - a; b++)
						{
							p3->information[0][b] = p->information[1][p->length2 - 1 - b];
							p3->chirality[0][b] = p->chirality[1][p->length2 - 1 - b];
						}
						p->information[1][a] = 0;

						p3->information[0][p->length2 - a] = 0;
						p3->information[1][0] = 0;

						p3->length1 = p->length2 - a;
						p->length2 = a;
						p3->length2 = 0;

						p3->prior = room_head[1][y][x];
						p3->next = room_head[1][y][x]->next;
						if (p3->next != room_head[1][y][x])(p3->next)->prior = p3;
						room_head[1][y][x]->next = p3;
					}
				}
			}

			// Nucleotide decaying and RNA degradating 
			for (p = room_head[1][y][x]->next; p != room_head[1][y][x]; p = p->next)
			{
				if (p->length1 == 1)	  // Nucleotide decaying
				{
					if (p->length2 == 0)
					{
						if (RANDD< PND)
						{
							if (p->chirality[0][0] == -1) lpn_arr[0][y][x]++; // Nucleotides decay into precursors with the same handedness
							else dpn_arr[0][y][x]++;

							(p->prior)->next = p->next;
							if (p->next != room_head[1][y][x])(p->next)->prior = p->prior;
							p3 = p;
							p = p->prior;
							free(p3);
						}
					}
					else if (p->length2 == 1)
					{
						if (RANDD< PND*sqrt(PND)) // The decay of the paired nucleotides at the same time
						{
							if (p->chirality[0][0] == -1) lpn_arr[0][y][x]++; // Nucleotides decay into precursors with the same handedness
							else dpn_arr[0][y][x]++;
							if (p->chirality[1][0] == -1) lpn_arr[0][y][x]++;
							else dpn_arr[0][y][x]++;

							(p->prior)->next = p->next;
							if (p->next != room_head[1][y][x])(p->next)->prior = p->prior;
							p3 = p;
							p = p->prior;
							free(p3);
						}
					}
					else {printf("unexpected error on chain-length"); exit(0);} 
				}
				else  // RNA-end decaying and RNA Degradating
				{
					// Nucleotide residue decaying at the end of RNA
					if (p->length1 == p->length2) {
						if (RANDD < PNDE*sqrt(PNDE)) { // The decay of the paired residues at the same time
							if (p->chirality[0][p->length1 - 1] == -1) lpn_arr[0][y][x]++;// nucleotide decay into the same handedness precursor
							else dpn_arr[0][y][x]++;
							p->information[0][p->length1 - 1] = 0;
							p->length1--;

							if (p->chirality[1][p->length2 - 1] == -1) lpn_arr[0][y][x]++;
							else dpn_arr[0][y][x]++;
							p->information[1][p->length2 - 1] = 0;
							p->length2--;
						}
					}
					else if (p->length1 > p->length2) { // Only the end of first chain may decay 
						if (RANDD < PNDE) {
							if (p->chirality[0][p->length1 - 1] == -1) lpn_arr[0][y][x]++;// nucleotide decay into the same handedness precursor
							else dpn_arr[0][y][x]++;
							p->information[0][p->length1 - 1] = 0;
							p->length1--;
						}
					}
					else {printf("unexpected error on chain-length"); exit(0);} 

					// Nucleotide residue decaying at the start of RNA
					if (p->length2 > 0) {
						if (RANDD < PNDE*sqrt(PNDE)) { // The decay of the paired residues at the same time
							if (p->chirality[0][0] == -1) lpn_arr[0][y][x]++;
							else dpn_arr[0][y][x]++;
							for (b = 1; b < p->length1; b++) {
								p->information[0][b - 1] = p->information[0][b];
								p->chirality[0][b - 1] = p->chirality[0][b];
							}
							p->information[0][p->length1 - 1] = 0;
							p->length1--;

							if (p->chirality[1][0] == -1) lpn_arr[0][y][x]++;
							else dpn_arr[0][y][x]++;
							for (b = 1; b < p->length2; b++) {
								p->information[1][b - 1] = p->information[1][b];
								p->chirality[1][b - 1] = p->chirality[1][b];
							}
							p->information[1][p->length2 - 1] = 0;
							p->length2--;
						}
					}
					else if (p->length2 == 0) {
						if (RANDD < PNDE) {  // The decay of the start residue on this single chain 
							if (p->chirality[0][0] == -1) lpn_arr[0][y][x]++;
							else dpn_arr[0][y][x]++;

							for (b = 1; b < p->length1; b++) {
								p->information[0][b - 1] = p->information[0][b];
								p->chirality[0][b - 1] = p->chirality[0][b];
							}
							p->information[0][p->length1 - 1] = 0;
							p->length1--;
						}
					}
					else {printf("unexpected error on chain-length"); exit(0);} 

					// RNA length check
					if (p->length1 == 0) {
						if (p->length2 != 0) {printf("unexpected error on chain-length"); exit(0);} 
						(p->prior)->next = p->next;
						if (p->next != room_head[1][y][x]) (p->next)->prior = p->prior;
						p3 = p;
						p = p->prior;
						free(p3);
						continue;
					}
					else if (p->length1 < p->length2) {printf("unexpected error on chain-length"); exit(0);} 

					// RNA degradation
					while (1)
					{
						for (j = p->length1; j>1; j--)  // Phosphodiester bonds are checked one by one.
						{
							f1 = PBB;					

							if (j <= p->length2 && p->nick[j - 1] == 0) f1 = f1*sqrt(f1);  // Falling into double chain region 

							if (RANDD < f1)
							{
								p3 = (struct rna *)malloc(LEN);
								if (!p3) { printf("\t%lddeg--memeout\n", i); exit(-1); }

								for (b = 0; b < p->length1 - j + 1; b++) {
									p3->information[0][b] = p->information[0][b + j - 1];
									p3->chirality[0][b] = p->chirality[0][b + j - 1];
								}
								p3->information[0][p->length1 - j + 1] = 0;
								p->information[0][j - 1] = 0;
								p3->length1 = p->length1 - j + 1;
								p->length1 = j - 1;

								if (p->length2>j - 1)
								{
									for (b = 0; b < p->length2 - j + 1; b++) {
										p3->information[1][b] = p->information[1][b + j - 1];
										p3->chirality[1][b] = p->chirality[1][b + j - 1];
									}
									p3->information[1][p->length2 - j + 1] = 0;
									p->information[1][j - 1] = 0;
									p3->length2 = p->length2 - j + 1;
									p->length2 = j - 1;
								}
								else
								{
									p3->information[1][0] = 0;
									p3->length2 = 0;
								}

								if (p3->length2 != 0)
								{
									p3->nick[0] = 1;
									for (a = 1; a<p3->length2; a++)
										p3->nick[a] = p->nick[j + a - 1];
								}

								p3->prior = room_head[1][y][x];
								p3->next = room_head[1][y][x]->next;
								if (p3->next != room_head[1][y][x])(p3->next)->prior = p3;
								room_head[1][y][x]->next = p3;
								break;
							}
						}
						if (j == 1) break;
					}

				}
			}

			for (p = room_head[1][y][x]->next; p != room_head[1][y][x]; p = p->next) fresh_rna(1); // Updating nucleotides and RNA for the next time step's checking

			//--------------------------------  Nucleotide precursors' events
			// Nucleotide precursors to Nucleotides
			int lnt_turn = 0;
			int dnt_turn = 0;

			for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
			{	// Finding NSR in the room
				// The value of flag1 represents the chirality of the NSR
				if (p->length1 > 2 * nsrlength || p->length2 != 0) continue;
									// Only a single chain containing the NSR domain but no longer than twice its length may act as an NSR
				flag1 = is_nsr_chirality(p); // The value of flag1 represents the chirality type of the NSR.
				if (flag1 != 0) {
					if (flag1 == -1) lnt_turn++;  // Counting possible catalysis-turns in the room 
					else dnt_turn++;
				}
			}

			//L-type
			pn_bef = lpn_arr[1][y][x];
			lpn_arr[1][y][x] = 0;
			for (k = 0; k<pn_bef; k++)
			{
				flagntsyn = 1;
				if (lnt_turn>0)
				{
					lnt_turn--;
					if (RANDD<PNFR) flagntsyn = 0; // Catalyzed synthesis
				}
				if (RANDD<PNF) flagntsyn = 0; // Uncatalyzed synthesis

				if (flagntsyn == 0)
				{				      // Forming a nucleotide 
					p3 = (struct rna *)malloc(LEN);
					if (!p3) { printf("\t%d form_monomer--memeout\n", k + 1); exit(-1); }
					randnt = RandomInteger % (4) + 1;
					switch (randnt)
					{
						case 1:  p3->information[0][0] = A; break;
						case 2:  p3->information[0][0] = C; break;
						case 3:  p3->information[0][0] = G; break;
						case 4:  p3->information[0][0] = U; break;
						default: printf("form randnt error");
					}
					p3->information[0][1] = 0;
					p3->information[1][0] = 0;
					p3->length1 = 1;
					p3->length2 = 0;
					p3->chirality[0][0] = -1;//L

					p3->prior = room_head[0][y][x];
					p3->next = room_head[0][y][x]->next;
					if (p3->next != room_head[0][y][x])(p3->next)->prior = p3;
					room_head[0][y][x]->next = p3;
				}
				else lpn_arr[0][y][x]++;
			}
			//D-type
			pn_bef = dpn_arr[1][y][x];
			dpn_arr[1][y][x] = 0;
			for (k = 0; k<pn_bef; k++)
			{
				flagntsyn = 1;
				if (dnt_turn>0)
				{
					dnt_turn--;
					if (RANDD<PNFR) flagntsyn = 0; // Catalyzed synthesis
				}
				if (RANDD<PNF) flagntsyn = 0; // Uncatalyzed synthesis

				if (flagntsyn == 0)
				{				      // Forming a nucleotide 
					p3 = (struct rna *)malloc(LEN);
					if (!p3) { printf("\t%d form_monomer--memeout\n", k + 1); exit(-1); }
					randnt = RandomInteger % (4) + 1;
					switch (randnt)
					{
						case 1:  p3->information[0][0] = A; break;
						case 2:  p3->information[0][0] = C; break;
						case 3:  p3->information[0][0] = G; break;
						case 4:  p3->information[0][0] = U; break;
						default: printf("form randnt error");
					}
					p3->information[0][1] = 0;
					p3->information[1][0] = 0;
					p3->length1 = 1;
					p3->length2 = 0;
					p3->chirality[0][0] = 1;//D

					p3->prior = room_head[0][y][x];
					p3->next = room_head[0][y][x]->next;
					if (p3->next != room_head[0][y][x])(p3->next)->prior = p3;
					room_head[0][y][x]->next = p3;
				}
				else dpn_arr[0][y][x]++;
			}

			// L-nucleotide precursors <--> D-nucleotide precursors
			int lpn_bef = lpn_arr[0][y][x];
			int dpn_bef = dpn_arr[0][y][x];
			dpn_arr[0][y][x] = 0;
			lpn_arr[0][y][x] = 0;
			for (k = 0; k < lpn_bef; k++) {
				if (RANDD < PCIC) dpn_arr[0][y][x]++;//Chirality changed
				else lpn_arr[0][y][x]++;
			}
			for (k = 0; k < dpn_bef; k++) {
				if (RANDD < PCIC) lpn_arr[0][y][x]++;//Chirality changed
				else dpn_arr[0][y][x]++;
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////
void freepool(void)        // Memory releasing  
{
	int m;
	for (m = 0; m<2; m++)
	{
		for (y = 0; y<N; y++)
		{
			for (x = 0; x<N; x++)
			{
				while (1)
				{
					if (room_head[m][y][x]->next != room_head[m][y][x])
					{
						p = room_head[m][y][x]->next;
						room_head[m][y][x]->next = p->next;
						free(p);
					}
					else break;
				}
				free(room_head[m][y][x]);
			}
		}
	}

}

/////////////////////////////////////////////////////////////////////////
void record(void)        // Data recording at every interval step (RECINT) 
{
	FILE* fptxt;
	if ((fptxt = fopen(RECTXT, "at")) == NULL) { printf("cannot open file"); exit(-1); }

	total_nt_mat[g] = 0;   // Total materials in quotient of nucleotides (including nucleotide precursors, nucleotides and nucleotide residues within RNA)
	pn_num[g] = 0;         // Number of nucleotide precursors
	rna_num[g] = 0;        // Number of nucleotides and RNAs

	nsr_num[g] = 0;       // Number of RNA molecules containing the NSR domain
	d_nsr_num[g] = 0;	  // Number of RNA molecules containing the D-NSR domain
	l_nsr_num[g] = 0;	  // Number of RNA molecules containing the D-NSR domain

	rep_num[g] = 0;      // Number of RNA molecules containing the REP domain
	d_rep_num[g] = 0;    // Number of RNA molecules containing the D-REP domain
	l_rep_num[g] = 0;    // Number of RNA molecules containing the L-REP domain

	int lpn_num = 0;     // Number of L-nucleotide precursors
	int dpn_num = 0;	 // Number of D-nucleotide precursors
	total_l = 0;         // Total material of L-type in quotient of nucleotides (including nucleotide precursors, nucleotides and nucleotide residues within RNA)
	total_d = 0;         // Total material of D-type in quotient of nucleotides (including nucleotide precursors, nucleotides and nucleotide residues within RNA)
	int mixed_rna = 0;   // Number of RNA molecules containing both D- and L-type nucleotide residues

	for (y = 0; y<N; y++)
	{
		for (x = 0; x<N; x++)
		{
			lpn_num += lpn_arr[0][y][x];
			dpn_num += dpn_arr[0][y][x];

			for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
			{
				rna_num[g]++;
				if (p->chirality[0][0] != p->chirality[0][p->length1 - 1]) mixed_rna++;
				total_nt_mat[g] += p->length1 + p->length2;

				for (k = 0; k < p->length1; k++) {
					if (p->chirality[0][k] == -1) total_l++;
					else total_d++;
				}
				for (k = 0; k < p->length2; k++) {
					if (p->chirality[1][k] == -1) total_l++;
					else total_d++;
				}
				int a = is_nsr_chirality(p);
				if (a != 0)
				{
					nsr_num[g]++;
					if (a == 1)d_nsr_num[g]++;
					else if (a == -1)l_nsr_num[g]++;
					else {printf("unexpected error on the chirality of Nsr"); exit(0);} 
				}

				a = is_rep_chirality(p);
				if (a != 0)
				{
					rep_num[g]++;
					if (a == 1)d_rep_num[g]++;
					else if (a == -1)l_rep_num[g]++;
					else {printf("unexpected error on the chirality of Rep"); exit(0);}
				}
			}
		}
	}
	pn_num[g] = lpn_num + dpn_num;
	total_l += lpn_num;
	total_d += dpn_num;
	total_nt_mat[g] += pn_num[g];

	if (i == 0) save_parameters(RECTXT);	
	printf("Step=>%d\n", i);
	printf("lpn:%-7d    dpn:%-7d    total_l:%-7d    total_d:%-7d\n", lpn_num, dpn_num, total_l, total_d);
	printf("Tn:%-7d     RNA:%-7d (mixed: %-7d)      Pn:%-7d \n", total_nt_mat[g], rna_num[g], mixed_rna, pn_num[g]);
	printf("Nsr:%-7d   D_Nsr:%-7d   L_Nsr:%-7d \n", nsr_num[g], d_nsr_num[g], l_nsr_num[g]);
	printf("Rep:%-7d   D_Rep:%-7d   L_Rep:%-7d \n", rep_num[g], d_rep_num[g], l_rep_num[g]);

	fprintf(fptxt, "Step=>%d\n", i);
	fprintf(fptxt, "lpn:%-7d    dpn:%-7d    total_l:%-7d    total_d:%-7d \n", lpn_num, dpn_num, total_l, total_d);
	fprintf(fptxt, "Tn:%-7d     RNA:%-7d (mixed: %-7d)      Pn:%-7d \n", total_nt_mat[g], rna_num[g], mixed_rna, pn_num[g]);
	fprintf(fptxt, "Nsr:%-7d   D_Nsr:%-7d   L_Nsr:%-7d \n", nsr_num[g], d_nsr_num[g], l_nsr_num[g]);
	fprintf(fptxt, "Rep:%-7d   D_Rep:%-7d   L_Rep:%-7d \n", rep_num[g], d_rep_num[g], l_rep_num[g]);

	// length distribution
	printf("length    :");
	fprintf(fptxt, "length    :");
	for (int sc = 0; sc < 50; sc++)  // Only RNA molecules shorter than 50 nt are recorded here.
	{
		printf("%14d,", sc + 1);
		fprintf(fptxt, "%14d,", sc + 1);
	}
	printf("\n");
	fprintf(fptxt, "\n");

	int chirality, ulength;
	for (int si = 0; si<MAX_RNA_LENGTH; si++)
	{
		chnum_d[g][si] = 0;
		chnum_l[g][si] = 0;
	}
	for (y = 0; y<N; y++)
	{
		for (x = 0; x<N; x++)
		{
			for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
			{
				if (check_rna_chirality(p, &chirality, &ulength))    // For D-L dimers, the value is false, i.e., they are not counted.
				{
					if (chirality == 1) chnum_d[g][ulength - 1]++; 
					// D-type RNA: valid length (that is, there might actually be one L-type residue at one end of the RNA, and if so, the chain-length recorded here would be the length of the D-domain)
					else if (chirality == -1) chnum_l[g][ulength - 1]++; 
					// L-type RNA: valid length (that is, there might actually be one D-type residue at one end of the RNA, and if so, the chain-length recorded here would be the length of the L-domain)
				}
			}
		}
	}
	printf("val_len   :");
	fprintf(fptxt, "val_len   :");
	for (int sc = 0; sc < 50; sc++)   // Only RNA molecules shorter than 50 nt are recorded here.
	{
		printf("%6d(%6d),", chnum_l[g][sc], chnum_d[g][sc]);    // L-type (D-type)
		fprintf(fptxt, "%6d(%6d),", chnum_l[g][sc], chnum_d[g][sc]);
	}
	printf("\n\n");
	fprintf(fptxt, "\n\n");

	g++;
	fclose(fptxt);

}

/////////////////////////////////////////////////////////////////////////

int main()
{
	inits();	// Initialization of the system
	for (i = 0; i <= STEPNUM; i++)	// The Monte-Carlo cycle
	{
		if (i >= STAREC&&i%RECINT == 0) record(); // Data recording at every interval step (RECINT) 
		unit_action();	// Action of units (molecules) in the system
	}
	freepool(); // Memory releasing
}
