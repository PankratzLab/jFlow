package org.genvisis.cnv.plots.data;


public interface Pipe {
	boolean hasNextPipe();

	Pipe getNextPipe();

	Pipe getPrevPipe();

	void setNextPipe(Pipe p);

	void setPrevPipe(Pipe p);

	String pipeValue(String value) throws RejectedValueException;
}