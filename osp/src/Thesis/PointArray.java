package Thesis;

public class PointArray {
	private double[][] array;

	
	
	
	public PointArray(int N){
		array = new double[N][2];
	}
	
	
	
	public PointArray(double[][] array){
		this.array= new double[array.length][2];
		this.array=array;
	}
	
	
	
	public double[] get(int i){
		return array[i];
	}
	
	
	
	
	public void set(int i, double[] point){
		array[i]=point;		
	}
	
	public int length(){
		return array.length;
	}
	
	
	
	public double[][] getArray(){
		return array;
	}
	
	public double[][] getArrayDeepCopy(){
		double[][] copy = new double[array.length][2];
		for(int i=0; i<copy.length; i++){
			for(int j=0;j<copy[i].length;j++){
				copy[i][j]=array[i][j];
			}
		}
		return copy;
	}
	
	
	public PointArray add(double input){
		PointArray tempArray = new PointArray(this.getArrayDeepCopy());
		for(int i=0; i < array.length; i++ ){
			tempArray.get(i)[1]+=input;
		}
		
		return tempArray;
	}
	
	public PointArray add(PointArray input){
		PointArray tempArray = new PointArray(this.getArrayDeepCopy());
		for(int i=0; i < array.length; i++ ){
			tempArray.get(i)[1]+=input.get(i)[1];
		}
		
		return tempArray;
	}
	
	public PointArray multiply(double m){
		PointArray tempArray = new PointArray(this.getArrayDeepCopy());
		for(int i=0; i < array.length; i++ ){
			tempArray.get(i)[1]*=m;
		}
		
		return tempArray;
	}
	
	public PointArray multiply(PointArray input){
		PointArray tempArray = new PointArray(this.getArrayDeepCopy());
		PointArray tempInput = new PointArray(input.getArrayDeepCopy());
		for(int i=0; i < array.length; i++ ){
			tempArray.get(i)[1]*=tempInput.get(i)[1];
		}
		
		return tempArray;
	}
	
	
	
	
	
	
	
	
	public int[] xSearch( double x){  // returns the indices of the two point in the searchable array that come before and after the input on the x axis.
		int i = 0;
		
		while(x >= array[i][0] && i<(array.length-1)){
			i++;
			if(i==array.length){
				i--;
			}
		}

		
		int[] indexArray = new int [2];
		if(i==0){
			indexArray[0] = i;
			indexArray[1] = i+1;
		}else{
			indexArray[0] = i-1;
			indexArray[1] = i;
		}

		return indexArray;
	}
	
	
	
	
	
	
	
	public double aproximateY(double x){
		int[] indices = xSearch(x);
		
		double[] point1=array[indices[0]];
		double[] point2=array[indices[1]];
		double slope = (point2[1]-point1[1])/(point2[0]-point1[0]);
		double aprxY = slope*(x-point1[0])+point1[1];
		return aprxY;
	}
	
	
	
	
	
	
	public void matchPoints(PointArray input){
		double[][] tempArray= new double[input.length()][2];
		for(int i=0; i<input.length(); i++){
			tempArray[i][0]=input.get(i)[0];
			tempArray[i][1]=aproximateY(input.get(i)[0]);
		}
		array=tempArray;
	}
	
	
	
	
	
	public PointArray derivative(){
		
		PointArray derivative = new PointArray(array.length-1);
		
		for(int i = 0; i<(derivative.length()); i++){
			derivative.get(i)[1]= (array[i][1]-array[i+1][1])/(array[i][0]-array[i+1][0]);
			derivative.get(i)[0]= (array[i][0]+array[i+1][0])/2.0;
		}
		
		derivative.matchPoints(this);
		return derivative;
		
	}
	
	
public PointArray secondDerivative(){
		
		PointArray derivative = new PointArray(array.length);
		
		for(int i = 1; i<(derivative.length())-1; i++){
			double dxInverse =2/(array[i+1][0]-array[i-1][0]);
			double d2y=(array[i+1][1]+array[i-1][1]-2.0*array[i][1]);
			derivative.get(i)[1]= d2y*dxInverse*dxInverse;
			derivative.get(i)[0]= array[i][0];
		}
		double aproxSlope =(array[2][1]-array[1][1])/(array[2][0]-array[1][0]) ;
		derivative.get(0)[1]=aproxSlope*(0-array[1][0])+array[1][1];
		derivative.get(0)[0]= array[0][0];
		
		aproxSlope =(array[array.length-2][1]-array[array.length-3][1])/(array[array.length-2][0]-array[array.length-3][0]) ;
		derivative.get(derivative.length()-1)[1]=aproxSlope*(0-array[array.length-3][0])+array[array.length-3][1];
		derivative.get(derivative.length()-1)[0]= array[array.length-1][0];
		
		return derivative;
		
	}





	public double chiSqr(PointArray input){ // bases the chi^2 of of the points of this PointArray
		PointArray tempArray = new PointArray(input.getArrayDeepCopy());
		tempArray.matchPoints(this);
	
		double chiSqr = 0;
	
		for(int i = 0; i<array.length; i++){
			
			double difference = array[i][1]-tempArray.get(i)[1];
			chiSqr+=difference*difference;
		}
		return chiSqr;
	}
	
	public double definiteIntegral(int ind0,int indf){
		double areaSum=0;
		double dx = (array[indf][0]-array[ind0][0])/(double)(indf-ind0+1);
		for(int i=ind0; i<=indf; i++){
			areaSum+=array[i][1]*dx;
		}
		
		return areaSum;
	}
	
	public void renormalize(double width, double norm){ //assumes even spacing between points and renormalizes the function so that the area underneath is equal to the norm variable
		PointArray tempArray = new PointArray(this.getArrayDeepCopy());
		
		double normConst=Math.sqrt(norm/tempArray.multiply(tempArray).definiteIntegral(0,array.length-1));
		for(int i=0; i<array.length; i++){
			array[i][1]*=normConst;
		}
	}
	

}
