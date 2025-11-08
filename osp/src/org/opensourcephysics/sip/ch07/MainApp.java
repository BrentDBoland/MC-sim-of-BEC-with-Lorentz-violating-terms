package org.opensourcephysics.sip.ch07;


public class MainApp {
	public static void main(String[] args) {
	    Poisson poisson = new Poisson(20);
	    
	    
	    double unity = poisson.expectation(new Unity(),200);
	    double n = poisson.expectation(new Identity(),200);
	    double n2 = poisson.expectation(new Square(), 200);
	    System.out.println("unity = "+unity+"   n = "+n+"   n2 = "+n2);
	    System.out.println(poisson.probability(1000));
	  }

}
