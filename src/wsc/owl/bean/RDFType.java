package wsc.owl.bean;

import javax.xml.bind.annotation.XmlAttribute;

import wsc.data.pool.NamespaceManager;

public class RDFType {

	private String resource;

	@XmlAttribute(namespace=NamespaceManager.RDF_NAMESPACE)
	public String getResource() {
		return resource;
	}

	public void setResource(String resource) {
		this.resource = resource;
	}

}
