package wsc.data.pool;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import org.jgrapht.DirectedGraph;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;
import wsc.graph.ServiceEdge;
import wsc.graph.ServiceInput;
import wsc.graph.ServiceOutput;

public class SWSPool {

	private List<Service> serviceList = new LinkedList<Service>();
	private SemanticsPool semantics;

	// private final Map<String, Service> graphOutputSetMap = new HashMap<String,
	// Service>();

	public List<Service> getServiceList() {
		return serviceList;
	}

	public void setServiceList(List<Service> serviceList) {
		this.serviceList = serviceList;
	}

	public SemanticsPool getSemantics() {
		return semantics;
	}

	public void setSemantics(SemanticsPool semantics) {
		this.semantics = semantics;
	}

	/**
	 * Semantic web service pool initialization
	 *
	 * @param semantics
	 * @param servicefile
	 *            Path
	 *
	 * @return SemanticWebServicePool
	 */
	// public static SWSPool parseXML(SemanticsPool semantics, String
	// serviceFilePath)
	// throws FileNotFoundException, JAXBException {
	// SWSPool swsp = new SWSPool();
	// swsp.semantics = semantics;
	// List<double[]> list = initialQoSfromSLA("qos.xml");
	//
	// Definitions def = Definitions.parseXML(serviceFilePath);
	// for (int i = 0; i < def.getSemExtension().getSemMessageExtList().size(); i +=
	// 2) {
	// swsp.serviceList.add(Service.initialServicefromMECE(def.getSemExtension().getSemMessageExtList().get(i),
	// def.getSemExtension().getSemMessageExtList().get(i + 1)));
	// }
	//
	// // manually add QoS attributes
	// for (int i = 0; i < list.size(); i++) {
	// swsp.qosServiceMap.put(swsp.serviceList.get(i).getServiceID(), list.get(i));
	// swsp.serviceList.get(i).setQos(list.get(i));
	// }
	//
	// System.out.println("No.of Service:" + swsp.serviceList.size());
	// return swsp;
	// }

	public static SWSPool parseWSCServiceFile(SemanticsPool semantics, String fileName) {

		SWSPool swsp = new SWSPool();
		swsp.semantics = semantics;

		List<ServiceInput> inputs = new ArrayList<ServiceInput>();
		List<ServiceOutput> outputs = new ArrayList<ServiceOutput>();
		double[] qos = new double[4];

		try {
			File fXmlFile = new File(fileName);
			DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
			DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
			Document doc = dBuilder.parse(fXmlFile);

			NodeList nList = doc.getElementsByTagName("service");

			for (int i = 0; i < nList.getLength(); i++) {
				org.w3c.dom.Node nNode = nList.item(i);
				Element eElement = (Element) nNode;
				// service name, for example serv904934656
				String name = eElement.getAttribute("name");

				qos[0] = Double.valueOf(eElement.getAttribute("Res"));
				qos[1] = Double.valueOf(eElement.getAttribute("Pri"));
				qos[2] = Double.valueOf(eElement.getAttribute("Ava"));
				qos[3] = Double.valueOf(eElement.getAttribute("Rel"));

				// Get inputs, instance name, for example inst995667695
				org.w3c.dom.Node inputNode = eElement.getElementsByTagName("inputs").item(0);
				NodeList inputNodes = ((Element) inputNode).getElementsByTagName("instance");
				for (int j = 0; j < inputNodes.getLength(); j++) {
					org.w3c.dom.Node in = inputNodes.item(j);
					Element e = (Element) in;
					ServiceInput serviceInput = new ServiceInput(e.getAttribute("name"), false);
					inputs.add(serviceInput);
				}

				// Get outputs instance name, for example inst1348768777
				org.w3c.dom.Node outputNode = eElement.getElementsByTagName("outputs").item(0);
				NodeList outputNodes = ((Element) outputNode).getElementsByTagName("instance");
				for (int j = 0; j < outputNodes.getLength(); j++) {
					org.w3c.dom.Node out = outputNodes.item(j);
					Element e = (Element) out;
					ServiceOutput serviceOutput = new ServiceOutput(e.getAttribute("name"), false);
					outputs.add(serviceOutput);
				}

				Service ws = new Service(name, qos, inputs, outputs);
				swsp.serviceList.add(ws);

				inputs = new ArrayList<ServiceInput>();
				outputs = new ArrayList<ServiceOutput>();
				qos = new double[4];
			}

		} catch (IOException ioe) {
			System.out.println("Service file parsing failed...");
		} catch (ParserConfigurationException e) {
			System.out.println("Service file parsing failed...");
		} catch (SAXException e) {
			System.out.println("Service file parsing failed...");
		}
		return swsp;
	}

	/**
	 * find a single service that can be applied now and update the output list and
	 * delete the service
	 *
	 * @param inputSet
	 */
	public Service findPossibleService(HashSet<String> inputSet) {
		int foundServiceIndex = -1;
		for (int i = 0; i < this.serviceList.size(); i++) {
			Service service = this.serviceList.get(i);

			if (service.searchServiceMatchFromInputSet(this.semantics, inputSet)) {
				foundServiceIndex = i;
				break;
			}
		}
		if (foundServiceIndex == -1) {
			// System.out.println("no matching for inputSet");
			return null;
		}
		Service service = this.serviceList.get(foundServiceIndex);
		this.serviceList.remove(foundServiceIndex);
		for (ServiceOutput output : service.getOutputList()) {

			if (!inputSet.contains(output.getOutput())) {
				inputSet.add(output.getOutput());
			}
		}
		return service;
	}

	/**
	 * find a single service that can be applied now and update the output list and
	 * delete the service
	 *
	 * @param inputSet
	 */
	public List<Service> findPossibleService4Layers(HashSet<String> inputSet) {

		List<Service> services4Layer = new ArrayList<Service>();

		for (int i = 0; i < this.serviceList.size(); i++) {
			Service service = this.serviceList.get(i);

			if (service.searchServiceMatchFromInputSet(this.semantics, inputSet)) {
				services4Layer.add(service);
			}
		}

		this.serviceList.removeAll(services4Layer);

		for (Service service : services4Layer) {
			for (ServiceOutput output : service.getOutputList()) {

				if (!inputSet.contains(output.getOutput())) {
					inputSet.add(output.getOutput());
				}
			}
		}

		return services4Layer;
	}

	/**
	 * find a single service that can be applied now and update the output list and
	 * delete the service
	 *
	 * @param graphOutputList
	 */
	public Service createGraphService(List<ServiceOutput> graphOutputs, List<Service> serviceCandidates,
			SemanticsPool semanticsPool, DirectedGraph<String, ServiceEdge> directedGraph,
			Map<String, Service> graphOutputListMap) {
		int foundServiceIndex = -1;

		for (int i = 0; i < serviceCandidates.size(); i++) {
			Service service = serviceCandidates.get(i);
			if (service.searchServiceGraphMatchFromInputSet(graphOutputs, semanticsPool, service, directedGraph,
					graphOutputListMap)) {
				foundServiceIndex = i;
				break;
			}
		}
		if (foundServiceIndex == -1) {
			System.out.println("no matching for inputSet");
			return null;
		}
		Service service = serviceCandidates.get(foundServiceIndex);
		serviceCandidates.remove(service);
		// add found service outputs to inputSet
		for (ServiceOutput output : service.getOutputList()) {
			// if (!graphOutputList.contains(output.getOutput())) {
			graphOutputs.add(output);
			// graphOutputList.add(output.getOutput());
			// output mapped back to service
			graphOutputListMap.put(output.getOutput(), service);
			// }
		}

		return service;
	}

	public Service createGraphService4Mutation(List<ServiceOutput> graphOutputs, List<Service> serviceCandidates,
			SemanticsPool semanticsPool, DirectedGraph<String, ServiceEdge> directedGraph,
			Map<String, Service> graphOutputListMap, List<String> ioNodeInputs) {
		int foundServiceIndex = -1;

		for (int i = 0; i < serviceCandidates.size(); i++) {
			Service service = serviceCandidates.get(i);
			if (service.searchServiceGraphMatchFromDefinedInputSet(graphOutputs, semanticsPool, service, directedGraph,
					graphOutputListMap, ioNodeInputs)) {
				foundServiceIndex = i;
				break;
			}
		}
		if (foundServiceIndex == -1) {
			System.out.println("no matching for inputSet");
			return null;
		}
		Service service = serviceCandidates.get(foundServiceIndex);
		serviceCandidates.remove(service);
		// add found service outputs to inputSet
		for (ServiceOutput output : service.getOutputList()) {
			// if (!graphOutputList.contains(output.getOutput())) {
			graphOutputs.add(output);
			// graphOutputList.add(output.getOutput());
			// output mapped back to service
			graphOutputListMap.put(output.getOutput(), service);
		}
		// }
		return service;
	}

}
