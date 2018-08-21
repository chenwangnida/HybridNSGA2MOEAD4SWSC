package wsc.graph;

public class ParamterConn {
	double matchType;
	double similarity;
	String outputInst;
	String outputrequ;
	String SourceServiceID;
	String TargetServiceID;
	//this parameter is set up to store the target Service information for mutation
	String OriginalTargetServiceId;
	
	boolean isConsidered;
	boolean isSetOriTargetSerId;

	public double getSimilarity() {
		return similarity;
	}

	public void setSimilarity(double similarity) {
		this.similarity = similarity;
	}

	public String getSourceServiceID() {
		return SourceServiceID;
	}

	public String getTargetServiceID() {
		return TargetServiceID;
	}

	public void setTargetServiceID(String targetServiceID) {
		TargetServiceID = targetServiceID;
	}

	public void setSourceServiceID(String sourceServiceID) {
		SourceServiceID = sourceServiceID;
	}

	public double getMatchType() {
		return matchType;
	}

	public String getOutputInst() {
		return outputInst;
	}

	public void setOutputInst(String outputInst) {
		this.outputInst = outputInst;
	}

	public String getOutputrequ() {
		return outputrequ;
	}

	public void setOutputrequ(String outputrequ) {
		this.outputrequ = outputrequ;
	}

	public void setMatchType(double matchType) {
		this.matchType = matchType;
	}

	public boolean isConsidered() {
		return isConsidered;
	}

	public void setConsidered(boolean isConsidered) {
		this.isConsidered = isConsidered;
	}

	public String getOriginalTargetServiceId() {
		return OriginalTargetServiceId;
	}

	public void setOriginalTargetServiceId(String originalTargetServiceId) {
		OriginalTargetServiceId = originalTargetServiceId;
	}

	public boolean isSetOriTargetSerId() {
		return isSetOriTargetSerId;
	}

	public void setSetOriTargetSerId(boolean isSetOriTargetSerId) {
		this.isSetOriTargetSerId = isSetOriTargetSerId;
	}

}
