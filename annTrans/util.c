#include <R_ext/Utils.h>

void cOverlaps(int* length, int* qret, int* sret, int* qstart, int* qwidth, int* qnum, int* sstart, int* swidth, int* snum, double* overlap) {
	int k = 0;
	int qoverlap, soverlap;
	double send, qend;
	for (int i = 0;i < *qnum;i++) {
		qoverlap = qwidth[i] * *overlap;
		qend = qstart[i] + qwidth[i] - 1;
		for (int j = 0;j < *snum;j++) {
			send = sstart[j] + swidth[j] - 1;
			soverlap = swidth[j] * *overlap;
			if (sstart[j] <= qstart[i]) {
                                if (send >= qend || (send - qstart[i]) >= qoverlap || (send - qstart[i]) >= soverlap) {
                                        qret[k] = i + 1;
					sret[k] = j + 1;
                                        k++;
                                }
                        } else if (send <= qend || (qend - sstart[j]) >= qoverlap || (qend - sstart[j]) >= soverlap) {
                                        qret[k] = i + 1;
					sret[k] = j + 1;
                                        k++;
                        }
	
		}
	}
	*length = k;
}
