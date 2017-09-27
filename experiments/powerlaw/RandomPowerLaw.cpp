#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#ifndef MAXSIZE
#define MAXSIZE 1000
#endif

char AName[64] = "A0.smat";
char GName[64] = "A0-orig.smat";
char BName[64] = "B0.smat";
char LName[64] = "L0.data";
char SName[64] = "S0.smat";

double p[MAXSIZE], alpha = 0.00, beta = 0.4, Gamma = 1.8;
int numreps=10;
double flipA = 0, flipL = 0;
int n = 200, r, deg[MAXSIZE], deg_copy[MAXSIZE], a[4*MAXSIZE], tot, edges;
int g[MAXSIZE][MAXSIZE], g1[MAXSIZE][MAXSIZE], g2[MAXSIZE][MAXSIZE], t1[MAXSIZE*MAXSIZE], t2[MAXSIZE*MAXSIZE];

int findsum(int k, int l, int r)
{
	if (l == r) return a[k] = deg[l];
	int m = (l + r) / 2;
	return a[k] = findsum(2 * k + 1, l, m) + findsum(2 *  k + 2, m + 1, r);
}

int findindex(int v)
{
	int k = 0, l = 0, r = n - 1, m;
	while (l != r){
		m = (l + r) / 2;
		k = 2 * k + 1;
		if (a[k] >= v) r = m;
		else v -= a[k], k++, l = m + 1;
	}
	return l;
}

void update(int ind, int delta)
{
	int k = 0, l = 0, r = n - 1, m;
	while (l != r){
		a[k] += delta;
		m = (l + r) / 2;
		k = 2 * k + 1;
		if (ind <= m) r = m;
		else k++, l = m + 1;
	}
	a[k] += delta;
}

void OutputLinkedList()
{
	int i, j;
	for (i = 0; i < n; i++){
		printf("%d:",i);
		for (j = 0; j < n; j++){
			if (g[i][j]) printf(" %d", j);
		}
		printf("\n");
	}
}

int RandGen()
{
	int k;
	if (RAND_MAX == 32768) {
		k = (rand() * 32768 + rand()) % tot + 1;
	}
	else {
		k = rand() % tot + 1;
	}
	return findindex(k);
}

double coinFlip()
{
	return 1.0 * rand() / RAND_MAX;
}

void initPowerLog(int k, int n, double r) {
	int i;
	double sum = 0;
	for (i = 1; i <= n; i++) {
		p[i] = 1.0 * k / pow(i, r);
		sum += p[i];
	}
	for (i = 1; i <= n; i++) {
		p[i] /= sum;
	}
}

int genPowerLog(int n) {
	int i;
	double k = 1.0 * rand() / RAND_MAX, sum = 0;
	for (i = 1; i <= n; i++) {
		sum += p[i];
		if (sum > k) return i;
	}
	return n;
}

void genDeg() {
	int i, tot = 0;
	initPowerLog(n, 30, Gamma);
	for (i = 0; i < n; i++) {
		deg[i] = genPowerLog(50);
		tot += deg[i];
	}
	if (tot & 1) tot++, deg[0]++;
	edges = tot / 2;
	memcpy(deg_copy, deg, sizeof(deg));
}

void write_graph(FILE* f, int g[][1000]) {
    int i, j, k = edges;
    fprintf(f, "%d %d %d\n", n, n, k);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (g[i][j]) {
				fprintf(f, "%d %d 1\n", i, j);
			}
		}
	}
    fclose(f);
}
    

void generate(FILE *f, int g[][1000]) {
	int i, j, k = edges;

	for (i = 0; i < n; i++) {
	  for (j = i + 1; j < n; j++) {
	    if (g[i][j] == 1 && coinFlip() < flipA) {
	      g[i][j] = g[j][i] = 0;
	      k--;
	    }
	  }
	}
	for (i = 0; i < n; i++) {
		for (j = i + 1; j < n; j++) {
			if (coinFlip() < alpha) {
				g[i][j] = g[j][i] = 1;
				k++;
			}
		}
	}
	fprintf(f, "%d %d %d\n", n, n, k);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (g[i][j]) {
				fprintf(f, "%d %d 1\n", i, j);
			}
		}
	}
}

void constructL() {
	int i, j;
	memset(g, 0, sizeof(g));
	for (i = 0; i < n; i++) {
		if (coinFlip() < (1 - flipL)) g[i][i] = 1;
		for (j = 0; j < n; j++) {
			if (i != j && coinFlip() < beta) g[i][j] = 1;
		}
	}
	FILE *f = fopen(LName,"w");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (g[i][j]) {
				fprintf(f, "%d %d 1\n", i, j);
				t1[r] = i;
				t2[r] = j;
				r++;
			}
		}
	}
	fclose(f);
}

void findSquares() {
	int i, j, k = 0;
	for (i = 0; i < r; i++) {
		for (j = 0; j < r; j++) {
			if (g1[t1[i]][t1[j]] && g2[t2[i]][t2[j]]) k++;
		}
	}
	FILE *f = fopen(SName, "w");
	fprintf(f, "%d %d %d\n", r, r, k);
	for (i = 0; i < r; i++) {
		for (j = 0; j < r; j++) {
			if (g1[t1[i]][t1[j]] && g2[t2[i]][t2[j]]) fprintf(f, "%d %d 1\n", i, j);
		}
	}
	fclose(f);
}

void read_params(int argc, char** argv) {
        // n p 
        // alpha = 0.00, beta = 0.4, Gamma = 1.8;
        // double flipA = 0.1, flipL = 0.5;
        // int n = 200, 
        if (argc > 1) {
                n = atoi(argv[1]);
        }
        if (argc > 2) {
                Gamma = atof(argv[2]);
        }
        if (argc > 3) {
                alpha = atof(argv[3]);
        }
        if (argc > 4) {
                beta = atof(argv[4]);
        }
        if (argc > 5) {
                numreps =  atoi(argv[5]);
        }
        if (argc > 6) {
                flipA = atof(argv[6]);
        }
        if (argc > 7) {
                flipL =  atof(argv[7]);
        }
        printf("%7s = %i\n", "n", n);
        printf("%7s = %f\n", "gamma", Gamma);
        printf("%7s = %f\n", "alpha", alpha);
        printf("%7s = %f\n", "beta", beta);
        printf("%7s = %f\n", "flipA", flipA);
        printf("%7s = %f\n", "flipL", flipL);
        printf("%7s = %i\n", "numreps", numreps);
        
        if (n>=1000) {
            fprintf(stderr,"invalid n, n must be smaller than %i\n",
                    MAXSIZE);
            exit(-1);
        }
}

int main(int argc, char** argv) {
    read_params(argc, argv);
	int i, j, k, t, z;
	srand(time(0));
	FILE *f1, *f2, *gfile;
	for (z = 0; z < numreps; z++) {
        snprintf(AName, 63, "A%i.smat", z);
        snprintf(BName, 63, "B%i.smat", z);
        snprintf(GName, 63, "A%i-orig.smat", z);
        snprintf(LName, 63, "L%i.data", z);
        snprintf(SName, 63, "S%i.smat", z);
		f1 = fopen(AName, "w");
		f2 = fopen(BName, "w");
        gfile = fopen(GName, "w");

		genDeg();

/*		
		for (i = 0; i < n; i++) {
			printf("%d\n", deg[i]);
		}
		return 0;
*/		

		for (t = 0;; t++){
			tot = 0;
			memcpy(deg, deg_copy, sizeof(deg));
			tot = edges * 2;
		
			findsum(0, 0, n-1);
			memset(g, 0, sizeof(g));

			k = 0;
			int succ = 1;
			while(tot){
				i = RandGen();
				j = RandGen();
				if (i!=j && !g[i][j] && coinFlip() < (1 - .025 * deg_copy[i] * deg_copy[j] / edges)){
					k = 0;
					tot -= 2;
					deg[i]--, deg[j]--;
					g[i][j] = g[j][i] = 1;
					update(i, -1);
					update(j, -1);
				}else{
					k++;
				}
				if (k == 100){
					for(i = 0; i < n; i++){
						for(j = i + 1; j < n; j++){
							if (deg[i] && deg[j] && !g[i][j]) break;
						}
						if (j < n) break;
					}
					if (i == n){
						succ = 0;
						break;
					}
					else {
						k = 0;
					}
				}
			}
			if (succ){
				break;
			}
		}
        
        write_graph(gfile, g);

		r = 0;
		memcpy(g1, g, sizeof(g));
		memcpy(g2, g, sizeof(g));

		generate(f1, g1);
		generate(f2, g2);
		constructL();
		findSquares();
		fclose(f1);
		fclose(f2);
	}
	return 0;
}
