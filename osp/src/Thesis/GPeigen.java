package Thesis;

public class GPeigen implements SimpleQauntumSystemSolver {

	double hBar = GPEigenDisplayApp.hBar;//.000000000000000658211928; // units of eV
	double m; // Mass of the particles in amu
	double omega; // strength of the harmonic potential
	double a; // effective scattering length
	double A;

	double F;
	int N; // number of the highest level eigen function
	int points; // number of data points used for the simulation.
	HarmEigen eigens;
	public PointArray funct; // The function that will be modified funct[n][0] = xn  funct[n][1] = yn
	PointArray V;
	double width; // width of the domain
	int accepted; // number of accepted changes
	int rejected; // number of rejected changes
	double norm;// the standard normalization. Choosing a number larger than 1 helps negate effects caused by limited memory of the double type.
	double cstepMax;
	double xstepMax;
	double Eexp;
	double E2exp;
	double error;
	static boolean metro=false;
	
	public GPeigen(int N, int points, double m, double omega, double width, double a, double[] AF, double norm, double cstepMax, double xstepMax, double error){
		this.N=N;
		this.points=points;
		this.m=m;
		this.omega=omega;
		this.width=width;
		this.a=a;
		A = AF[0];

		F = AF[1];
		this.norm=norm;
		this.cstepMax=cstepMax;
		this.xstepMax=xstepMax;
		this.error=error;
		accepted = 0;
		rejected = 0;
		eigens = new HarmEigen(N, hBar, m, omega);
		funct = new PointArray(points);
		V = new PointArray(points);
		if(metro){
			update();
		}else{
			double xStep=Math.pow(10,-width)/(double)(points-1);

			for(int i=0; i<points; i++){
				double x=xStep*i-Math.pow(10,-width)*.5;
				funct.get(i)[0]=x;
				double y =eigens.evaluate(x);
				funct.get(i)[1]=y;
			}
			


		    PointArray tempArray = hamiltonian(new PointArray(funct.getArrayDeepCopy())); 
	    	Eexp =tempArray.multiply(funct).definiteIntegral(0,funct.length()-1);
	    	E2exp=Eexp*Eexp;

		}
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

		PointArray Vpsi = new PointArray(points);
		
		double xStep=Math.pow(10,-width)/(double)(points-1);
		
		for(int i = 0; i<array.length(); i++){
			double Vx;
			double x = array.get(i)[0];
			if(metro){
				 Vx= omega*omega*m*(x*x)/2.0;
			}else{
				
				
				Vx = omega*omega*m*((x+.5*xStep)*(x+.5*xStep)*(x+.5*xStep)-(x-.5*xStep)*(x-.5*xStep)*(x-.5*xStep))/(6.0*xStep);
			}
			
			
			V.get(i)[0]=x;
			V.get(i)[1]=Vx;
		}
		final PointArray psi = new PointArray(array.getArrayDeepCopy());
		PointArray Gp = psi.multiply(psi).multiply(4*Math.PI*hBar*hBar*a/m);
		V=V.add(Gp);
		Vpsi=psi.multiply(V);
		
		double[][] d2psi=array.getArrayDeepCopy();
		double[][] d2psiCopy=array.getArrayDeepCopy();
		for(int i=0; i<d2psi.length;i++){
			d2psi[i][1]=eigens.d2Evaluate(d2psi[i][0]);
			d2psiCopy[i][1]=eigens.d2Evaluate(d2psi[i][0]);
		}
		PointArray Fpsi = (new PointArray(d2psi)).multiply((-1)*hBar*hBar/(2*m)*F);
		
		PointArray lorVio = Fpsi.add(psi.multiply(A));
		
		
		PointArray tempArray = (new PointArray(d2psiCopy)).multiply((-1)*hBar*hBar/(2*m)).add(Vpsi).add(lorVio);
		
		return tempArray;
	}
	
	
	
	
	
		
	public void step() {
	    for(int i = 0;i<N; ++i){
	    	double[] coeffsave = eigens.coeffDeepCopy();
	    	PointArray tempArray = new PointArray(funct.getArrayDeepCopy());
	    	eigens.coeffStep(cstepMax);
	    	
	    	double[] EandE2 =energy(tempArray);
	    	
	    	double Eexp =EandE2[0];
	    	double E2exp =EandE2[1];
	    	double std = Math.sqrt(E2exp-Eexp*Eexp);
	    	double std2 = Math.sqrt(this.E2exp-this.Eexp*this.Eexp);

	      if(metro){
	    	  if((this.Eexp-Math.max(std,std2)> Eexp)&&(std<Eexp*error*100)) {
		      		this.Eexp = Eexp;
		      		this.E2exp = E2exp;
		      		funct=tempArray;
		      		accepted++;
		    	  
		      	}else{
		      		eigens.coeff=coeffsave;
		      		rejected++;
		      	}
	    	  
	      }else{
	    	  if(this.Eexp> Eexp) {
		      		this.Eexp = Eexp;
		      		this.E2exp = E2exp;
		      		funct=tempArray;
		      		accepted++;
		      		for(int j=0; j<eigens.coeff.length; j++){
		    			System.out.println(eigens.coeff[j]);
		    		}
		      		System.out.println();
		    	  
		      	}else{
		      		eigens.coeff=coeffsave;
		      		rejected++;
		      	}
	      }
	      	
	    }

	  }
	
	private double[] energy(PointArray array){
		double runs=0;
    	double EexpSum=0;
    	double E2expSum=0;
    	if(metro){
    		do{
        		double x;

        		do{
        			
        			
        			
        			x=0+xstepMax*(2*Math.random()-1);
        		}while(eigens.evaluate(x)*eigens.evaluate(x)/eigens.evaluate(0)/eigens.evaluate(0)<Math.random());
        		array.get(0)[0]=x;
        		array.get(0)[1]=eigens.evaluate(x);
        		
        		for(int j=1; j<(int)points; j++){

        			
        			x =array.get(j-1)[0];
    				double xp;
    				do{
    					
    					xp=x+xstepMax*(2*Math.random()-1);
    				}while(eigens.evaluate(xp)*eigens.evaluate(xp)/eigens.evaluate(x)/eigens.evaluate(x)<Math.random());
    				array.get(j)[0]=xp;
    				array.get(j)[1]=eigens.evaluate(xp);
    			}
    	    	PointArray hamTemp = hamiltonian(array);
    	    	double E_Lsum=0;
    	    	for(int j=0;j<points;j++){
    	    		double E_L=hamTemp.get(j)[1]/array.get(j)[1];
    	    		E_Lsum+=E_L;
    	    	}
    	    	double E=E_Lsum/(double)points;
    	    	EexpSum+=E;
    	    	E2expSum+=E*E;
    	    	runs++;

        	}while(runs<10);
    	}else{
    		for(int j=0; j<(int)points; j++){

    			
    			double x =array.get(j)[0];
				array.get(j)[1]=eigens.evaluate(x);
			}
    		runs=1;
    		PointArray tempArray = hamiltonian(new PointArray(array.getArrayDeepCopy())); 
	    	EexpSum =tempArray.multiply(array).definiteIntegral(0,array.length()-1);
	    	E2expSum=EexpSum*EexpSum;
    		
    	}


    	
    	double[] EandE2 = new double[2];
    	EandE2[0]=EexpSum/runs;
    	EandE2[1]=E2expSum/runs;
    	
    	if(EandE2[1]<EandE2[0]*EandE2[0]){
    		EandE2[1]=EandE2[0]*EandE2[0];
    	}
    	
    	return EandE2;
	}
	
	
	
	
	
	public void renormalize(){ //assumes even spacing between points and renormalizes the function so that the area underneath is equal to the norm variable
		funct.renormalize(10, norm);
	}
	
	
	public double[][] realValues(){
		return funct.getArray();
	}
	
	public double[][] realValuesDeep(){
		return funct.getArrayDeepCopy();
	}
	
	public double[][] getVdeep(){
		return V.getArrayDeepCopy();
	}
	
	public double getAcceptanceRate(){
		return 100*(double)accepted/(double)(accepted+rejected);
	}
	
	public double getEexp(){
		return Eexp;
	}

	public double getChiSqr(){
		return 1.0/0.0;
	}
	
	public void save(){
		try {
			DataWriter.addToFile(eigens);
			DataWriter.addToFile(this);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	
	public void load(){
		double[] temp=DataReader.HarmEigenCoeff();
		if(temp != null){
			eigens.coeff=temp;
			N = temp.length;
			eigens.hermite = new Hermite[N];
			for(int i=0; i<eigens.hermite.length; i++){
				eigens.hermite[i]= new Hermite(i,m,omega);
			}
		}
		
		double accepted=DataReader.GPeigenConst("accepted");
		if(accepted != Double.NaN){
			this.accepted=(int) accepted;
		}
		
		double rejected=DataReader.GPeigenConst("rejected");
		if(rejected != Double.NaN){
			this.rejected=(int) rejected;
		}
		
		double Eexp=DataReader.GPeigenConst("Eexp");
		if(Eexp != Double.NaN){
			this.Eexp=Eexp;
		}
		
		update();
	}
	
	
	
	
	private void update(){

		double[] EandE2 =energy(funct);
		Eexp =EandE2[0];
		E2exp =EandE2[1];
	}
	
	public void setAF(double[] AF){
		A = AF[0];
		F = AF[1];
	}
	
	public void resetAJ(){
		accepted=0;
		rejected=0;
		update();
	}
	
	public void demo(int n_c){
		for(int i=0; (i<eigens.coeff.length)&&(n_c<eigens.coeff.length);i++){
			  if(i==n_c){
				  eigens.coeff[i]=1;
			  }else{
				  eigens.coeff[i]=0;
			  }
		  }
		update();
	}
	

}
