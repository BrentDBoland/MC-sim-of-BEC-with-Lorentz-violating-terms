package Thesis;

public interface SimpleQauntumSystemSolver {
	
	abstract double[][] eigenVal(PointArray psiInitial, PointArray psiFinal);
	
	abstract PointArray hamiltonian(PointArray array);
		
	abstract void step();
	
	abstract void renormalize();
	
	abstract double[][] realValues();
	
	abstract double[][] realValuesDeep();
	
	abstract double getAcceptanceRate();
	
	abstract double getChiSqr();
	
	public double getEexp();
	
}
