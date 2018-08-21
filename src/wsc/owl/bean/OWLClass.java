package wsc.owl.bean;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;

import wsc.data.pool.NamespaceManager;
public class OWLClass {

	private String ID;
	private OWLSubClassOf subClassOf;


	@XmlAttribute(name="ID", namespace=NamespaceManager.RDF_NAMESPACE)
	public String getID() {
		return ID;
	}


	public void setID(String iD) {
		ID = iD;
	}


	@XmlElement(name="subClassOf", namespace=NamespaceManager.RDFS_NAMESPACE)
	public OWLSubClassOf getSubClassOf() {
		return subClassOf;
	}


	public void setSubClassOf(OWLSubClassOf subClassOf) {
		this.subClassOf = subClassOf;
	}

}
