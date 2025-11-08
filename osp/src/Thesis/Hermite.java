package Thesis;

public class Hermite {
	public double[] coeff;
	int n;
	double m=1;
	double omega=1;
	boolean even;
	
	
	public Hermite(int n, double m, double omega){
		this.n=n;
		this.m=m;
		this.omega=omega;
		
		if(n%2==0){
			even=true;
			coeff = new double[n/2+1];
		}else{
			even=false;
			coeff = new double[(n-1)/2+1];
		}
		coeff[0]=1;
		for(int i=1;i<coeff.length;i++){
			int j;
			if(even){
				j=2*(i-1);
			}else{
				j=2*(i-1)+1;
			}
			coeff[i]=nextCoeff(coeff[i-1],j,n);
		}
		scaleCoeff();
	}
	
	
	private void scaleCoeff(){
		double scale=Math.pow(2,n)/(coeff[coeff.length-1]);
		for(int i=0; i<coeff.length; i++){
			coeff[i]*=scale;
		}
	}
	
	

	private double nextCoeff(double a,int j,int n){
		double nextA=-2.0*(n-j)/(double)((j+1)*(j+2))*a;
		return nextA;
	}
	
	
	public double evaluate(double x){
		double y=0;
		for(int i=0;i<coeff.length;i++){
			y+=coeff[i]*Math.pow(x,2*i);
		}
		if(!even){
			y*=x;
		}
		return y;
	}
	
	public double d1Evaluate(double x){
		double y=0;
		if(even){
			for(int i=0;i<coeff.length;i++){
				y+=(2*i)*coeff[i]*Math.pow(x,2*i-1);
			}
		}else{
			for(int i=0;i<coeff.length;i++){
				y+=(2*i+1)*coeff[i]*Math.pow(x,2*i);
			}		
		}
		
		return y;
	}
	
	public double d2Evaluate(double x){
		double y=0;
		if(even){
			for(int i=0;i<coeff.length;i++){
				y+=(2*i-1)*(2*i)*coeff[i]*Math.pow(x,2*i-2);
			}
		}else{
			for(int i=0;i<coeff.length;i++){
				y+=(2*i)*(2*i+1)*coeff[i]*Math.pow(x,2*i-1);
			}		
		}
		
		return y;
	}
}
