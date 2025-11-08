
package Thesis;

import org.opensourcephysics.controls.*;
import org.opensourcephysics.frames.*;



public class SimpleQuantumDisplayApp extends AbstractSimulation {
	int N;
	double m;
	double omega;
	double width;
	double norm;
	double stepMax;
	SimpleQauntumSystemSolver solver;
	PlotFrame plotFrame = new PlotFrame("x", "y", "function plot");
	PlotFrame plotFrame2 = new PlotFrame("x", "y", "function2 plot");
  
  

  
  public void initialize() {
	  N = control.getInt("N");
	  m = control.getDouble("mass");
	  omega = control.getDouble("omega");
	  width = control.getDouble("width");
	  norm = control.getDouble("norm");
	  stepMax = control.getDouble("stepMax");
	  solver = new HarmonicOscillator(N, m, omega, width, norm, stepMax);
	  
	  plotFrame.setConnected(true);
	  plotFrame.clearData();
	  plotFrame2.clearData();		
	  
	  PointArray funct = new PointArray(solver.realValuesDeep());
	  for(int n =0; n<funct.length(); n++){
		  plotFrame.append(0, funct.get(n)[0], funct.get(n)[1]);
	  }
	  
	  
	  PointArray tempArray = funct.derivative();							//
	  tempArray.renormalize(width,norm);									//
	  double[][] funct2 = tempArray.getArray();								// appends the 1st derivative of our function to the 2nd frame
	  for(int n =0; n<funct2.length; n++){									//
		  plotFrame2.append(0, funct2[n][0], funct2[n][1]);					//
	  }
	  
	  
	  
	  

	  tempArray=solver.hamiltonian(funct);									//
	  tempArray.renormalize(width,norm);									// appends the 2nd derivative of our function to the 1st frame
	  double[][] funct3 = tempArray.getArray();								//
	  for(int n =0; n<funct3.length; n++){									//
		  plotFrame.append(1, tempArray.get(n)[0], tempArray.get(n)[1]);	//
	  }																		//
	  
	  plotFrame.setMessage("acceptance rate = 0"+"  <E> = "+solver.getEexp());
	  plotFrame.repaint();
	  plotFrame2.repaint();

	  
  }


  public void doStep() {
    solver.step();
    plotFrame.clearData();
	  plotFrame2.clearData();		
	  
	  double[][] funct = solver.realValues();
	  for(int n =0; n<funct.length; n++){
		  plotFrame.append(0, funct[n][0], funct[n][1]);
	  }
	  
	  
	  PointArray tempArray = new PointArray(solver.realValuesDeep());        //
	  tempArray=tempArray.derivative();										//
	  tempArray.renormalize(width,norm);									//
	  double[][] funct2 = tempArray.getArray();								// appends the 1st derivative of our function to the 2nd frame
	  for(int n =0; n<funct2.length; n++){									//
		  plotFrame2.append(0, funct2[n][0], funct2[n][1]);					//
	  }
	  
	  
	  
	  

	  tempArray=solver.hamiltonian(new PointArray(solver.realValuesDeep()));						//
	  tempArray.renormalize(width,norm);									// appends the 2nd derivative of our function to the 1st frame
	  double[][] funct3 = tempArray.getArray();								//
	  for(int n =0; n<funct3.length; n++){									//
		  plotFrame.append(1, funct3[n][0], funct3[n][1]);					//
	  }		
	  
	  
   plotFrame.setMessage("acceptance rate = "+solver.getAcceptanceRate()+"  <E> = "+solver.getEexp());
   plotFrame.repaint();
   plotFrame2.repaint();													
  }


  public void reset() {
    control.setValue("N", 100);
    control.setValue("mass", 3728399991.6);
    control.setValue("omega", 1);
    
    control.setValue("width", 10.0);
    control.setValue("norm", 1.0);
    control.setValue("stepMax", .1);
  }
  

  public static void main(String[] args) {
    SimulationControl.createApp(new SimpleQuantumDisplayApp());
  }
}