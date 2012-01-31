package cnv.var;

import java.util.Hashtable;
import java.util.Vector;

public class IndiPheno {
	private double[] filters;
	private double[] covars;
	private int[] classes;
	private Vector<Hashtable<String,CNVariant[]>> cnvClasses;

	public IndiPheno() {}

	public void setFilters(double[] filters) {
		this.filters = filters;
	}

	public void setCovars(double[] covars) {
		this.covars = covars;
	}

	public void setClasses(int[] classes) {
		this.classes = classes;
	}

	public void setCNVclasses(Vector<Hashtable<String,CNVariant[]>> cnvClasses) {
		this.cnvClasses = cnvClasses;
	}

	public double[] getFilters() {
		return filters;
	}

	public double[] getCovars() {
		return covars;
	}

	public int[] getClasses() {
		return classes;
	}

	public CNVariant[] getCNVs(int cnvClass, int chr) {
		Hashtable<String,CNVariant[]> hash;
		
		if (cnvClass > cnvClasses.size()) {
			System.err.println("Error - specified cnvClass was greater than total number of cnvClasses");
			return null;
		}
		hash = cnvClasses.elementAt(cnvClass);
		if (hash.containsKey(chr+"")) {
			return hash.get(chr+"");
		} else {
			return null;
		}
	}
}
