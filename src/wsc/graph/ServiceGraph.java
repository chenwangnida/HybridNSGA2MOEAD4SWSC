package wsc.graph;


import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;

public class ServiceGraph extends DefaultDirectedWeightedGraph<String, ServiceEdge>{

	public ServiceGraph(Class<? extends ServiceEdge> edgeClass) {
		super(edgeClass);
		// TODO Auto-generated constructor stub
	}

	private static final long serialVersionUID = -1161512349678026436L;

	/**
	 *
	 */
//	public Map<String, GraphNode> nodeMap = new HashMap<String, GraphNode>();
//	public List<ServiceEdge> edgeList = new ArrayList<ServiceEdge>();

	@Override
	public java.lang.String toString() {
		StringBuilder builder = new StringBuilder();


		builder.append("digraph g {");

		for (ServiceEdge e :this.edgeSet()) {
			String serviceEdgeStr = String.format("%s->%s", this.getEdgeSource(e), this.getEdgeTarget(e));
			builder.append(serviceEdgeStr);
//			builder.append("mt:"+e.getAvgmt()+"dst"+e.getAvgsdt());
			builder.append("; ");
		}
		builder.append("}");
		return builder.toString();
	}
}
