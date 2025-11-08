package Thesis;

public class HarmonicOscillator implements SimpleQauntumSystemSolver {

	double m; // Mass of the particles in eV/c^2
	double omega; // strength of the harmonic potential
	int N; // number of points on our function considered
	public PointArray funct; // The function that will be modified funct[n][0] = xn  funct[n][1] = yn
	double width; // width of the domain
	int accepted; // number of accepted changes
	int rejected; // number of rejected changes
	double norm;// the standard normalization. Choosing a number larger than 1 helps negate effects caused by limited memory of the double type.
	double stepMax; // size is a fraction of the norm. Ex if norm was 100 and stepMax was .01 then the maximum change of any step two the whole system would be 1
	double Eexp;
	
	public HarmonicOscillator(int N, double m, double omega, double width, double norm, double stepMax){
		this.N=N;
		this.m=m;
		this.omega=omega;
		this.width=width;
		this.norm=norm;
		this.stepMax=stepMax;
		accepted = 0;
		rejected = 0;
		funct = new PointArray(N);
		
		double xStep=width/(double)(N-1);
		for(int i=0; i<N; i++){
			double[] point= new double[2];
			
			point[0]=xStep*i-width*.5;
			point[1]=1;//2*Math.random()-1;
			funct.set(i,point);
		}

		renormalize();
		PointArray newFunct = new PointArray(funct.getArrayDeepCopy());
	      
		PointArray tempArray = hamiltonian(newFunct);
	      
    	Eexp =tempArray.multiply(newFunct).definiteIntegral(0,tempArray.length()-1);
	      //System.out.println(chiSqr);
	}
	
	public double[][] eigenVal(PointArray psiInitial, PointArray psiFinal){
		double[][] eigenArray = new double[N][2];
		
		for(int i = 1; i<psiInitial.length()-1; i++){
				
				double eigenVal = psiFinal.get(i)[1]/psiInitial.get(i)[1];
				
				eigenArray[i][0]=psiInitial.get(i)[0];
				eigenArray[i][1]=eigenVal;
			}
		
			return eigenArray;
			
	}
	
	
	
	
	
	
	public PointArray hamiltonian(PointArray array){
		double dx = width/(double)(N-1);
		PointArray Vpsi = new PointArray(N);
		for(int i = 0; i<array.length(); i++){
			double x_minus = array.get(i)[0] -.5*dx;
			double x_plus = array.get(i)[0]+.5*dx;
			double V = omega*omega*m*((x_plus*x_plus*x_plus)-(x_minus*x_minus*x_minus))/(6*dx); // the average value of 1/2kx^2 about an interval of dx
			
			Vpsi.get(i)[0]=array.get(i)[0];
			Vpsi.get(i)[1]=V*array.get(i)[1];
		}
		
		double hBar = 1;//.000000000000000658211928; // units of eV
		
		PointArray tempArray = array.secondDerivative().multiply((-1)*hBar*hBar/(2*m)).add(Vpsi);
		
		return tempArray;
	}
	
	
	
	
	
		
	public void step() {
	    for(int i = 0;i<N; ++i){
	    	int randPoint;
	    	randPoint = (int) (Math.random()*N);
	    	double dy = (2.0*Math.random()-1.0)*stepMax;   // random change in y
	      
	    	PointArray newFunct = new PointArray(funct.getArrayDeepCopy());
	    	newFunct.get(randPoint)[1]+=dy;
	      
	    	PointArray tempArray = hamiltonian(newFunct);
	      
	    	double Eexp =tempArray.multiply(newFunct).definiteIntegral(0,tempArray.length()-1);

	      	//System.out.println(chiSqr+"  "+dy+"  "+randPoint);
	      
	      	if(this.Eexp > Eexp) {
	      		this.Eexp = Eexp;
	      		funct=newFunct;
	      		renormalize();
	      		accepted++;
	    	  
	      	}else{
	      		rejected++;
	      	}
	    }

	  }
	
	
	
	
	
	
	public void renormalize(){ //assumes even spacing between points and renormalizes the function so that the area underneath is equal to the norm variable
		funct.renormalize(width, norm);
	}
	
	
	public double[][] realValues(){
		return funct.getArray();
	}
	
	public double[][] realValuesDeep(){
		return funct.getArrayDeepCopy();
	}
	
	public double getAcceptanceRate(){
		return (double)accepted;//(double)(accepted+rejected);
	}
	
	public double getEexp(){
		return Eexp;
	}

	public double getChiSqr(){
		return 1.0/0.0;
	}
	

}
