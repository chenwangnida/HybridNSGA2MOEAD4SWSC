package wsc.owl.bean;

import javax.xml.bind.annotation.XmlAttribute;

public class Ontology {


	private String about;

	@XmlAttribute
	public String getAbout() {
		return about;
	}

	public void setAbout(String about) {
		this.about = about;
	}

}
