package Thesis;
import java.io.File;

import javax.xml.parsers.*;

import org.opensourcephysics.controls.SimControl;
import org.w3c.dom.*;

public class DataReader {
    
        static DocumentBuilderFactory factory;
        static DocumentBuilder builder;

        //Here we do the actual parsing
        static Document doc;

        public static void readFile(String path) throws Exception{
        	factory = DocumentBuilderFactory.newInstance();
        	builder = factory.newDocumentBuilder();
        	doc = builder.parse(new File(path));        	
        }
        
        
        
        
        
        public static double[] HarmEigenCoeff(){
        	Element el = doc.getDocumentElement();
            Node c = findChild(el, "cHarmEigen");
            double[] output = null;
            if(c != null){
            	NodeList cMap = c.getChildNodes();
                output = new double[cMap.getLength()];
                for (int i = 0; i < output.length; i++) {
                    output[i] = Double.parseDouble(cMap.item(i).getTextContent());
                }
            }
            
            return output;
        }
        
        
        public static double GPeigenConst(String constant){
        	Element el = doc.getDocumentElement();
            Node c = findChild(el, "GPeigen");
            double output = Double.NaN;
            if(c != null){
            	Node n =findChild(c, constant);
            	if(n != null){
            		output=Double.parseDouble(n.getTextContent());
            		}
            }
            
            return output;
        }
              
        
        
        
        
        private static Node findChild(Node parent, String name) {
        	Node child = null;
        	
        	if(parent.hasChildNodes()) {
        		NodeList nList=parent.getChildNodes();
        		for(int i = 0 ; i<nList.getLength(); i++){
                	if(nList.item(i).getNodeName().equals(name)){
                		child = nList.item(i);
                	}
                }
        		
        	}
        		return child;        		
        	
        }
        
        public static void setControlConst(SimControl control){
        	Element el = doc.getDocumentElement();
        	Node con = findChild(el,"constants");
        	if(con!=null){
        		NodeList conList = con.getChildNodes();
        		for(int i = 0 ; i<conList.getLength(); i++){
                	Node name = findChild(conList.item(i),"name");
                	Node value = findChild(conList.item(i),"value");
                	control.setValue(name.getTextContent(),value.getTextContent());
                }
        	}
        	
        	
        }
        
        
        
        
        
        
        
    }
