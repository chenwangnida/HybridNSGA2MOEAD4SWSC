package wsc.graph;

public class ServiceOutput {

	private String output;
	private String serviceId;
	boolean isSatified;

	public String getServiceId() {
		return serviceId;
	}

	public void setServiceId(String serviceId) {
		this.serviceId = serviceId;
	}

	public ServiceOutput(String output, boolean isSatified) {
		super();
		this.output = output;
		this.isSatified = isSatified;
	}

	public ServiceOutput(String output, String serviceId, boolean isSatified) {
		super();
		this.output = output;
		this.serviceId = serviceId;
		this.isSatified = isSatified;
	}

	public String getOutput() {
		return output;
	}

	public void setOutput(String output) {
		this.output = output;
	}

	public boolean isSatified() {
		return isSatified;
	}

	public void setSatified(boolean isSatified) {
		this.isSatified = isSatified;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((output == null) ? 0 : output.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		ServiceOutput other = (ServiceOutput) obj;
		if (output == null) {
			if (other.output != null)
				return false;
		} else if (!output.equals(other.output))
			return false;
		return true;
	}
	
	

//	@Override
//	public int hashCode() {
//		final int prime = 31;
//		int result = 1;
//		result = prime * result + (isSatified ? 1231 : 1237);
//		result = prime * result + ((output == null) ? 0 : output.hashCode());
//		return result;
//	}
//
//	@Override
//	public boolean equals(Object obj) {
//		if (this == obj)
//			return true;
//		if (obj == null)
//			return false;
//		if (getClass() != obj.getClass())
//			return false;
//		ServiceOutput other = (ServiceOutput) obj;
//		if (isSatified != other.isSatified)
//			return false;
//		if (output == null) {
//			if (other.output != null)
//				return false;
//		} else if (!output.equals(other.output))
//			return false;
//		return true;
//	}
}
