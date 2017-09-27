#include <cstdio>
#include <cstring>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>

using namespace std;

FILE *f1, *f2, *f3;

map<pair<int,int>,int> A, B;

int t1[5000], t2[5000];
int n1, n2, m1, m2;

int main(int argc, char **argv) {
  int i, j, k;
  f1 = fopen ("graph.EN.smat", "r");
  f2 = fopen ("graph.ES.smat", "r");
  fscanf(f1, "%d%*d%d", &n1, &m1);
  while (m1--) {
    fscanf(f1,"%d%d%*d",&i,&j);
    A[make_pair<int,int>(i,j)] = 1;
  }
  fscanf(f2, "%d%*d%d", &n2, &m2);
  while (m2--) {
    fscanf(f2, "%d%d%*d", &i, &j);
    B[make_pair<int,int>(i,j)] = 1;
  }
  f3 = fopen (argv[1], "r");
  k = 0;
  while(fscanf(f3, "%d%d", &i, &j) != EOF) {
    t1[k] = i;
    t2[k++] = j;
  }
  int ct = 0, mark;
  for (i = 0; i < k; i++) {
    mark = 0;
    for (j = 0; j < k; j++) {
      if (A.find(make_pair<int,int>(t1[i]-1,t1[j]-1)) != A.end()) {
	if (B.find(make_pair<int,int>(t2[i]-1,t2[j]-1)) != B.end()) {
	  ct++;
	  mark = 1;
	}
      }
    }
    if (mark) printf("%d %d\n",t1[i],t2[i]);
  }
  fprintf(stderr, "%d\n",ct);
  return 0;
}
