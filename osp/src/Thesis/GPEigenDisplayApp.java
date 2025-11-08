
package Thesis;

import java.io.File;

import org.opensourcephysics.controls.*;
import org.opensourcephysics.frames.*;



public class GPEigenDisplayApp extends AbstractSimulation {
	int N;
	int points;
	double m;// Mass of the particles in amu
	double omega;
	double width;
	double norm;
	double cstepMax;//<=1
	double xstepMax;
	double a;
	double error;
	
	double[] AF;
	int[] AFstep;
	double[] AFsSize;
	int currentVar;
	
	final static double a0 = 5.2917721092*Math.pow(10,-11);
	final static double hBar = 6.35077993441*Math.pow(10,-8);
	String loadFileName;
	GPeigen solver = null;
	PlotFrame plotFrame;
	PlotFrame hPlotFrame;
  
  

  
  public void initialize() {
	  
	  
	 
	  m = control.getDouble("mass in amu");
	  omega = control.getDouble("omega");
	  a=a0*control.getDouble("eff. scatt");
	  N = control.getInt("number of eigen states");
	  points = control.getInt("number of points");
	  width = control.getDouble("width 10m^-v");
	  norm = control.getDouble("norm");
	  cstepMax = control.getDouble("c step max");
	  xstepMax = control.getDouble("x step max");
	  error = control.getDouble("percent error");
	  currentVar = control.getInt("currentvar");
	  int n_c=control.getInt("demo: cn=1  n=");
	  
	  if(GPeigen.metro){
		  points=points/10;
	  }
	    
	  loadFileName = control.getString("loadFile");
	  
	  AF=new double[2];
	  AFstep=new int[2];
	  AFstep[currentVar]=0;
	  AFstep[(currentVar+1)%2]=control.getInt("var2");
	  AFsSize=new double[2];
	  AFsSize[0]=hBar*omega/2.0/100.0;;
	  AFsSize[1]=1.0/100.0;
	  
	  

	  
	  for(int i=0; i<2; i++){
		  AF[i]=AFsSize[i]*AFstep[i];
	  }


	  
	  solver = new GPeigen(N, points, m,omega, width, a, AF, norm, cstepMax, xstepMax, error) ;
	  solver.demo(n_c);
	  plotFrame = new PlotFrame("x m*10^-"+width, "psi", "psi plot");
	  hPlotFrame = new PlotFrame("x m*10^-"+width, "Hpsi", "Hpsi plot");

	  if(!loadFileName.equals("")){
		  try {
				DataReader.readFile("C:\\Users\\Public\\Documents\\Thesis\\Data\\GPeigen\\"+loadFileName+".xml");
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		  solver.load();
	  }
	  

	  plotFrame.clearData();
	  hPlotFrame.clearData();	
	  
	  
	  PointArray funct = new PointArray(solver.realValuesDeep());
	  PointArray tempArray=solver.hamiltonian(funct);
	  PointArray V = new PointArray(solver.getVdeep());
	  double[][] hPsi = tempArray.getArrayDeepCopy();	
	  
	  for(int n =0; n<funct.length(); n++){
		  plotFrame.append(0, funct.get(n)[0], funct.get(n)[1]);
		// hPlotFrame.append(0, 1.66053892*6.24150934*Math.pow(10,-9)*V.get(n)[0], V.get(n)[1]);	
	  }
	  
	  
	  
	  
	  

	  									
							
	  for(int n =0; n<hPsi.length; n++){									
		  hPlotFrame.append(1, tempArray.get(n)[0], 1.66053892*6.24150934*Math.pow(10,-9)*tempArray.get(n)[1]);	
				
	  }																		
	  
	  //plotFrame.setMessage("accepted = 0"+"  acceptance rate = 0");
	  hPlotFrame.setMessage("std = "+1.66053892*6.24150934*Math.pow(10,-9)*Math.sqrt(solver.E2exp-solver.Eexp*solver.Eexp)+"    "+solver.E2exp+" "+solver.Eexp*solver.Eexp+"  <E> = "+1.66053892*6.24150934*Math.pow(10,-9)*solver.getEexp()+" hwE="+solver.getEexp()/(hBar*omega),0);
	  plotFrame.repaint();
	  hPlotFrame.repaint();
	  solver.resetAJ();
	  
  }


  public void doStep() {
    solver.step();
    plotFrame.clearData();
    hPlotFrame.clearData();		
	  double[][] funct = solver.realValuesDeep();
	  PointArray tempArray=solver.hamiltonian(new PointArray(solver.realValuesDeep()));	
	  double[][] hPsi = tempArray.getArrayDeepCopy();
	  double[][] V = solver.getVdeep();
	  
	  
	  for(int n =0; n<funct.length; n++){
		  plotFrame.append(0, funct[n][0], funct[n][1]);
		  hPlotFrame.append(0, V[n][0], 1.66053892*6.24150934*Math.pow(10,-9)*V[n][1]);
	  }
				
	  for(int n =0; n<hPsi.length; n++){																	
		  hPlotFrame.append(1, hPsi[n][0], 1.66053892*6.24150934*Math.pow(10,-9)*hPsi[n][1]);													
	  }		
	  
	  
   plotFrame.setMessage("accepted = "+solver.accepted+"rejected = "+solver.rejected+"  acceptance rate = "+solver.getAcceptanceRate(),0);
   hPlotFrame.setMessage("std = "+1.66053892*6.24150934*Math.pow(10,-9)*Math.sqrt(solver.E2exp-solver.Eexp*solver.Eexp)+"    "+solver.E2exp+" "+solver.Eexp*solver.Eexp+"  <E> = "+1.66053892*6.24150934*Math.pow(10,-9)*solver.getEexp()+"   "+solver.getEexp()/(hBar*omega),0);
   plotFrame.repaint();
   hPlotFrame.repaint();
   
   if(solver.getAcceptanceRate()<.5&& (solver.rejected>500)){
	   nextPoten();
   }
  }


  public void reset() {

	    control.setValue("mass in amu", 100);
	    control.setValue("omega", 11*2*Math.PI);
	    control.setValue("eff. scatt", 100);
		control.setValue("number of eigen states", 10);
		control.setValue("number of points", 1000);
	    control.setValue("width 10m^-v", 4);
	    control.setValue("norm", 1.0);
	    control.setValue("percent error", .005);
	    control.setValue("c step max", .001);
	    control.setValue("x step max", .000007);
	    control.setValue("currentvar", 0);
	    control.setValue("var2", 0);
	    control.setValue("demo: cn=1  n=",0);
	    control.setValue("loadFile","");


  }
  
  public void save(String title){
	  String type;
	  if(GPeigen.metro){
		  type="GPeigenMet";
	  }else{
		  type="GPeigen";
	  }
	  try {
		DataWriter.newFile(title, type);
	} catch (Exception e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	}
	  try {
		solver.save();
		DataWriter.addToFile(control);
	} catch (Exception e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	}
  }
  
  public void load(){
	  loadFileName = control.getString("loadFile");
	  loader(loadFileName);
  }
  
  private void loader(String fileName){
	  String path;
	  if(solver.metro){
		  path="C:\\Users\\Public\\Documents\\Thesis\\Data\\GPeigenMet\\"+fileName+".xml";
	  }else{
		  path="C:\\Users\\Public\\Documents\\Thesis\\Data\\GPeigen\\"+fileName+".xml";
	  }
	  File file = new File(path);
	  if(file.isFile()){
		  try {
				DataReader.readFile(path);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			  DataReader.setControlConst(control);
			  control.refresh();
			  
			  solver = new GPeigen(N, points, m,omega, width, a, AF, norm, cstepMax, xstepMax, error) ;
			  solver.load();
	  }
	  
	  
	  
  }
  
  
  
  private void nextPoten(){
	  System.out.println("var="+(solver.E2exp-solver.Eexp*solver.Eexp));
	  save("A"+AFstep[0]+"F"+AFstep[1]);
	  varStep(currentVar);
	  solver.setAF(AF);
	  solver.resetAJ();
	  while(currentVar==2){
		  
	  }
  }
  
  private void varStep(int var){
	  int n=50;
	  if(AFstep[var]==n){
		  if(AFstep[(var+1)%2]==n){
			  
			  AFstep[var]=0;
			  AFstep[(var+1)%2]=1;
			  
			  currentVar++; 
			  
			  
			
			  loader("A"+0+"F"+0);
			  

			  
		  }else{
			  
				  AFstep[var]=0;
				  AFstep[(var+1)%2]=n;
				  
				  loader("A"+AFstep[0]+"F"+AFstep[1]);
				  
				  AFstep[var]=1;
				  
			  }
		  }else{
			  AFstep[var]++;
			  }
	  for(int i=0; i<2; i++){
		  AF[i]=AFsSize[i]*AFstep[i];
	  }
	  System.out.println("A="+AF[0]+" F="+AF[1]);
  }
  
  

  public static void main(String[] args) {
	
	  SimulationControl control = SimulationControl.createApp(new GPEigenDisplayApp());
	  control.addButton("load","load");
  }
}