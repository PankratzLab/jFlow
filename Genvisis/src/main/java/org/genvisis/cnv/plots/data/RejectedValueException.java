package org.genvisis.cnv.plots.data;

public class RejectedValueException extends Exception {
	Pipe rejector;

	public RejectedValueException(String rejection, Pipe rejector) {
		super(rejection);
		this.rejector = rejector;
	}
}
