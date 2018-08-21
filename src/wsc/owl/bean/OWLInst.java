package wsc.owl.bean;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;

import wsc.data.pool.NamespaceManager;
public class OWLInst {

	private String ID;
	private RDFType rdfType;

	@XmlAttribute(name="ID", namespace=NamespaceManager.RDF_NAMESPACE)
	public String getID() {
		return ID;
	}

	public void setID(String iD) {
		ID = iD;
	}

	@XmlElement(name="type", namespace=NamespaceManager.RDF_NAMESPACE)
	public RDFType getRdfType() {
		return rdfType;
	}

	public void setRdfType(RDFType rdfType) {
		this.rdfType = rdfType;
	}

}
