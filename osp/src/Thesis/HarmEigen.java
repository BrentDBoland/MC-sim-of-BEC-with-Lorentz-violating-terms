package Thesis;

public class HarmEigen {
	public double[] coeff;
	Hermite[] hermite;
	double hBar = GPEigenDisplayApp.hBar;
	double m;
	double omega;
	
	public HarmEigen(int N, double hBar, double m, double omega){
		this.hBar=hBar;
		this.m=m;
		this.omega=omega;
		hermite = new Hermite[N];
		for(int i=0; i<hermite.length; i++){
			hermite[i]= new Hermite(i,m,omega);
		}
		coeff=new double[N];

		coeff[0]=1;
		/*for(int i=0; i<coeff.length; i++){
			coeff[i]=Math.exp(-i);
		}*/
		
		double scale=0;
		for(int i=0; i<coeff.length; i++){
			scale+=coeff[i]*coeff[i];
		}
		for(int i=0; i<coeff.length; i++){
			coeff[i]=coeff[i]/Math.sqrt(scale);
		}
	}
	
	public void coeffStep(double stepMax){
		int randCoeff1 ,randCoeff2;
    	do{
		randCoeff1 = (int) (Math.random()*(coeff.length));
		}while((coeff[randCoeff1]+.5)<Math.random());
    	
    	do{
    		randCoeff2 = (int) (Math.random()*(coeff.length));
    		}while((randCoeff1==randCoeff2));
    	
    	
    		double deltaC1 = Math.random()*stepMax*coeff[randCoeff1];   // random change in x1 is greater than -x_1  so deltaCxi are real
        	
        	double b=2*coeff[randCoeff2];
        	double c=deltaC1*deltaC1-2*deltaC1*coeff[randCoeff1];
        	double deltaC2 =(-b+Math.sqrt(b*b-4*c))/2;
        	coeff[randCoeff1]-=deltaC1;
        	coeff[randCoeff2]+=deltaC2;
        	coeff[randCoeff2]*=2*Math.round(Math.random())-1;
    	
    	
    	
    	double scale=0;
		for(int i=0; i<coeff.length; i++){
			scale+=coeff[i]*coeff[i];
		}
		for(int i=0; i<coeff.length; i++){
			coeff[i]=coeff[i]/Math.sqrt(scale);
		}
	}
	
	public double[] coeffDeepCopy(){
		double[] copy = new double[coeff.length];
		for(int i=0; i<copy.length; i++){
			copy[i]=coeff[i];
		}
		return copy;
	}
	
	public double evaluate(double x){
		double xi=Math.sqrt(m*omega/hBar)*x;
		double sum=0;
		for(int i=0; i<coeff.length; i++){
			double normconst=Math.pow(m*omega/Math.PI/GPEigenDisplayApp.hBar,.25)/Math.sqrt(Math.pow(2,i))*sqrtFactInv(i);
			sum+=normconst*coeff[i]*hermite[i].evaluate(xi)*Math.exp(-xi*xi/2);
		}
		return sum;
	}
	
	public double d2Evaluate(double x){
		double xi=Math.sqrt(m*omega/hBar)*x;
		double sum=0;
		for(int i=0; i<coeff.length; i++){
			double normconst=Math.sqrt(Math.sqrt(m*omega/Math.PI/GPEigenDisplayApp.hBar))/Math.sqrt(Math.pow(2,i))*sqrtFactInv(i);
			sum+=normconst*coeff[i]*m*omega/hBar*(
					hermite[i].d2Evaluate(xi)*Math.exp(-xi*xi/2)
					-2*hermite[i].d1Evaluate(xi)*xi*Math.exp(-xi*xi/2)
					+hermite[i].evaluate(xi)*(xi*xi-1)*Math.exp(-xi*xi/2));
		}
		return sum;
	}
	
	private double sqrtFactInv(int m){   //(1/n!)^1/2
		if(m==0){
			return 1;
		}
		else{
			return 1/Math.sqrt(m)*sqrtFactInv(m-1);
		}
	}
}
