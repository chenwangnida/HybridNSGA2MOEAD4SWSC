package wsc.graph;

import java.util.ArrayList;
import java.util.List;

import org.jgrapht.graph.DefaultWeightedEdge;


public class ServiceEdge extends DefaultWeightedEdge implements Cloneable {

	String SourceService;
	String TargetService;
	List<ParamterConn> pConnList = new ArrayList<ParamterConn>();

	// average matching type value from source vertice to target vertice
	double avgmt;

	// average semantic distance value value from source vertice to target
	// vertice
	double avgsdt;

	public ServiceEdge(double avgmt, double avgsdt) {
		super();
		this.avgmt = avgmt;
		this.avgsdt = avgsdt;
	}



	@Override
	protected double getWeight() {
		// TODO Auto-generated method stub
		return super.getWeight();
	}


	@Override
	public Object getSource() {
		// TODO Auto-generated method stub
		return super.getSource();
	}

	@Override
	public Object getTarget() {
		// TODO Auto-generated method stub
		return super.getTarget();
	}

	public String getTargetService() {
		return TargetService;
	}

	public void setTargetService(String targetService) {
		TargetService = targetService;
	}

	public String getSourceService() {
		return SourceService;
	}

	public void setSourceService(String sourceService) {
		SourceService = sourceService;
	}

	public List<ParamterConn> getpConnList() {
		return pConnList;
	}

	public void setpConnList(List<ParamterConn> pConnList) {
		this.pConnList = pConnList;
	}

	public double getAvgmt() {
		return avgmt;
	}

	public void setAvgmt(double avgmt) {
		this.avgmt = avgmt;
	}

	public double getAvgsdt() {
		return avgsdt;
	}

	public void setAvgsdt(double avgsdt) {
		this.avgsdt = avgsdt;
	}

	@Override
	public ServiceEdge clone() {
		// TODO Auto-generated method stub
		return (ServiceEdge) super.clone();
	}

	@Override
	public String toString() {
		// return this.avgmt+";"+this.avgsdt + super.toString();
		return this.avgmt + ";" + this.avgsdt + super.toString();
		// return super.toString();
	}

}