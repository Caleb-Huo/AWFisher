#include <iostream>
#include <R.h>
#include <Rmath.h>
#include <random>
#include <vector>
#include <algorithm>

using namespace std;

template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

  return idx;
}
	
void importanceSampling(double a, int n, int k, double *bTarget, int J, double *awStats){
		
	std::vector<double> pTargets(bTarget, bTarget+J);	
	for(int j=0;j<J;j++){
		pTargets[j] = -log(pTargets[j]);
	}
	
	
	std::vector<double> betaVec(k,0.0);
	std::vector<double> awStat_vec(n,0.0);
	std::vector<double> weight_vec(n,0.0);
	std::vector<double> awStat_sort(n,0.0);
	std::vector<double> weight_sort(n,0.0);
	std::vector<double> logProbNew(n,0.0);

	double awStat;
	double pCumLog;
	double alog;
	double pImp;
	double cumSumExp=0.0;
	double apTarget;
	
	for(int i=0;i<n;i++){
		awStat = 0.0;
		pCumLog = 0.0;

		for(int ak=0;ak<k;ak++){
			betaVec[ak] = qbeta(((double) rand() / (RAND_MAX)), a, 1.0, 1, 0);			
			//betaVec[ak] = rbeta(a,1.0);
		}

			

		sort(betaVec.begin(),betaVec.end()); 
	
		for(int ak=0;ak<k;ak++){
			alog = log(betaVec[ak]);
			pCumLog += alog;
			awStat = max(awStat, -pchisq(-2*pCumLog, 2*(ak+1), 0, 1));

			//cout <<  "awStat: " << awStat << endl;
			//cout <<  "-pchisq(-2*pCumLog, 2*(ak+1), 0, 1): " << -pchisq(-2*pCumLog, 2*(ak+1), 0, 1) << endl;
		}


		awStat_vec[i] = awStat;
	    weight_vec[i] = - k*log(a) - (a-1) * pCumLog;
	    
	}
	
	int index = 0;
	for (auto i: sort_indexes(awStat_vec)) {
		awStat_sort[index] = awStat_vec[i];
		weight_sort[index++] = weight_vec[i];
	}
	
	for(int i=0;i<n;i++){
		cumSumExp += exp(weight_sort[i])/(n+1);
		logProbNew[i] = -log(cumSumExp);		
	}	
	
	for(int j=0;j<J;j++){
		apTarget = pTargets[j];
		for(int i=0;i<n;i++){
			if(apTarget>logProbNew[i]){
				awStats[j] = awStat_sort[i];
				break;
			}
		}
	}    				
	
}



extern "C" {
  void importanceSampling_R(double *a, int *n, int *k, double *pTarget, int *J, double *awStat)
  {
    importanceSampling(*a, *n, *k, pTarget, *J, awStat);
  }
}

