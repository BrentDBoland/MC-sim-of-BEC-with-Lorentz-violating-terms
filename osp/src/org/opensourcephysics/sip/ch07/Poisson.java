package org.opensourcephysics.sip.ch07;

public class Poisson {
	double nAve;
	
	public Poisson(double nAve){
		this.nAve = nAve;
	}
	
	public double probability(double n){
		double numerator = Math.pow(nAve,n)*Math.exp(-nAve);
		double nFactor = 1;
		for(int i = 1; i <= n; i++ ){
			nFactor*=i;
		}
		return numerator/nFactor;
	}
	
	
	/* f is the function you are finding the expectation of.
	 * N is the number of terms of the sum you want to include.
	 * The higher N is the more accurate the result.
	 */
	public double expectation(Function f, int N){   
		double value = 0;
		for(int n=0; n< N; n++){
			value += probability(n)*f.evaluate(n);
		}
		return value;
	}

}
