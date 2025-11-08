package Thesis;
import java.io.*;
import java.util.Iterator;

import javax.xml.parsers.*;
import javax.xml.transform.*;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.opensourcephysics.controls.SimControl;
import org.w3c.dom.*;

public class DataWriter {
	 static DocumentBuilderFactory factory;
	 static DocumentBuilder builder;
	 static Document doc;
	 static String path;


	 public static void newFile(String title, String type) throws Exception{
		 factory = DocumentBuilderFactory.newInstance();
	    builder = factory.newDocumentBuilder();
	    doc = builder.newDocument();
	    
	    Element el = doc.createElement(type);
        doc.appendChild(el);
	        
	    DOMSource source = new DOMSource(doc);
	    path = "C:\\Users\\Public\\Documents\\Thesis\\Data\\"+type+"\\"+title+".xml";
	    StreamResult result = new StreamResult(new File(path));
	        
	    TransformerFactory transformerFactory = TransformerFactory.newInstance();
	    Transformer transformer = transformerFactory.newTransformer();

	    transformer.transform(source, result);
	 }
	
    
    public static void addToFile(HarmEigen input) throws Exception{
    	factory = DocumentBuilderFactory.newInstance();
    	builder = factory.newDocumentBuilder();
    	doc = builder.parse(new File(path));
    	
    	Element type = doc.getDocumentElement();
        
        Node ch = doc.createElement("cHarmEigen");

        for(int i=0; i<input.coeff.length; i++){
        	Node ci = doc.createElement("c"+i);
        	ci.setTextContent(""+input.coeff[i]);
        	ch.appendChild(ci);
        }
        type.appendChild(ch);
        
        DOMSource source = new DOMSource(doc);
        StreamResult result = new StreamResult(new File(path));
        
        TransformerFactory transformerFactory = TransformerFactory.newInstance();
        Transformer transformer = transformerFactory.newTransformer();

        transformer.transform(source, result);
        
    }
    
    public static void addToFile(SimControl control) throws Exception{
    	factory = DocumentBuilderFactory.newInstance();
    	builder = factory.newDocumentBuilder();
    	doc = builder.parse(new File(path));
    	
    	Element type = doc.getDocumentElement();
        
        Node con = doc.createElement("constants");
        type.appendChild(con);
        
        Iterator<String> propNames = control.getPropertyNames().iterator();

        for(int i=0; propNames.hasNext(); i++){
        	String propName = propNames.next();
        	Node prop = doc.createElement("prop"+i);
        	Node name = doc.createElement("name");
        	Node value = doc.createElement("value");
        	
        	name.setTextContent(propName);
        	value.setTextContent(control.getObject(propName).toString());
        	prop.appendChild(name);
        	prop.appendChild(value);
        	con.appendChild(prop);
        }
        
        
        DOMSource source = new DOMSource(doc);
        StreamResult result = new StreamResult(new File(path));
        
        TransformerFactory transformerFactory = TransformerFactory.newInstance();
        Transformer transformer = transformerFactory.newTransformer();

        transformer.transform(source, result);
        
    }
    
    public static void addToFile(GPeigen input) throws Exception{
    	factory = DocumentBuilderFactory.newInstance();
    	builder = factory.newDocumentBuilder();
    	doc = builder.parse(new File(path));
    	
    	Element type = doc.getDocumentElement();
        
        Node ch = doc.createElement("GPeigen");

        {
        	Node accpt = doc.createElement("accepted");
        	Node rejt = doc.createElement("rejected");
        	Node Eexp = doc.createElement("Eexp");
        	Node E2exp = doc.createElement("E2exp");
        	
        	Node A = doc.createElement("A");
        	Node C = doc.createElement("C");
        	Node F = doc.createElement("F");
        	
        	accpt.setTextContent(""+input.accepted);
        	rejt.setTextContent(""+input.rejected);
        	Eexp.setTextContent(""+input.Eexp);
        	E2exp.setTextContent(""+input.E2exp);
        	
        	A.setTextContent(""+input.A);
        	C.setTextContent("0");
        	F.setTextContent(""+input.F);
        	
        	ch.appendChild(accpt);
        	ch.appendChild(rejt);
        	ch.appendChild(Eexp);
        	ch.appendChild(E2exp);
        	ch.appendChild(A);
        	ch.appendChild(C);
        	ch.appendChild(F);
        	
        }
        type.appendChild(ch);
        
        DOMSource source = new DOMSource(doc);
        StreamResult result = new StreamResult(new File(path));
        
        TransformerFactory transformerFactory = TransformerFactory.newInstance();
        Transformer transformer = transformerFactory.newTransformer();

        transformer.transform(source, result);
        
    }
    
    
    

}
