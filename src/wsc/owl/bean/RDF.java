/**
 *
 */
package wsc.owl.bean;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import wsc.data.pool.NamespaceManager;


@XmlRootElement(name="RDF", namespace="http://www.w3.org/1999/02/22-rdf-syntax-ns#")
public class RDF {

	private Ontology ontology;
	private List<OWLClass> owlClassList = new ArrayList<OWLClass>();
	private List<OWLInst> owlInstList = new ArrayList<OWLInst>();
	@XmlElement
	public Ontology getOntology() {
		return ontology;
	}

	public void setOntology(Ontology ontology) {
		this.ontology = ontology;
	}

	@XmlElement(name="Class", namespace=NamespaceManager.OWL_NAMESPACE)
	public List<OWLClass> getOwlClassList() {
		return owlClassList;
	}

	public void setOwlClassList(List<OWLClass> owlClassList) {
		this.owlClassList = owlClassList;
	}

	@XmlElement(name="Thing", namespace=NamespaceManager.OWL_NAMESPACE)
	public List<OWLInst> getOwlInstList() {
		return owlInstList;
	}

	public void setOwlInstList(List<OWLInst> owlInstList) {
		this.owlInstList = owlInstList;
	}

	/**
	 * 	Unmarshal XML data from the specified InputStream and return the resulting content tree.
	 *
	 * @param owl.file path
	 * @return RDF Object
	 */
	public static RDF parseXML(String filePath) throws JAXBException, FileNotFoundException{
		JAXBContext jc = JAXBContext.newInstance(RDF.class);
		return (RDF)jc.createUnmarshaller().unmarshal(new FileInputStream(new File(filePath)));

	}

}
