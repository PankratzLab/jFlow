package org.genvisis.cnv.var;

import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.filesys.CNVariant;

public class IndiPheno {
	private double[] filters;
	private double[] covars;
	private int[] classes;
	private Vector<Hashtable<String,CNVariant[]>> cnvClasses = new Vector<Hashtable<String,CNVariant[]>>();
	private static volatile boolean cnvsLoaded = false;
	
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

	public Vector<Hashtable<String, CNVariant[]>> getCnvClasses() {
		return cnvClasses;
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
	    return getCNVs(cnvClass, chr, true);
	}
	
	public CNVariant[] getCNVs(int cnvClass, int chr, boolean waitIfLoading) {
	    if (!cnvsLoaded) {
	        System.out.println("CNVs Not Loaded Yet!");
	        if (!waitIfLoading) {
	            return null;
	        } else {
	            System.out.println("Waiting for CNVs to be Loaded...");
	            while (!cnvsLoaded) {
	                //try {
	                	try {
							Thread.sleep(100);// was getting  java.lang.IllegalMonitorStateException with wait()
						} catch (InterruptedException ie) {
						}
                     //   wait();
//                    } catch (InterruptedException e) {
//                        e.printStackTrace();
//                    }
	            }
	        }
	    }
	    
        Hashtable<String,CNVariant[]> hash;
	
		if (cnvClass >= cnvClasses.size()) {
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

    public static void setCNVsLoaded() {
        cnvsLoaded = true;
    }
}
