package common;

import java.util.logging.Logger;

/**
 * Program to get heap statistic of a runtime
 * Author: Rohit Sinha
 * Date: 12/14/13.
 */
public class HeapStats {

	private static final int MB_SIZE = 1024*1024;
	private Runtime runtime;
	private Logger log;

	public HeapStats(){
		//Getting the runtime reference from system
		this.runtime = Runtime.getRuntime();
	}

	public HeapStats(java.util.logging.Logger log) {
		this();	// calling the default constructor
		this.log = log;
	}

	public long getTotalMemory(){
		log.info("Total heap size in this runtime is: " + runtime.totalMemory()/MB_SIZE);
		return (runtime.totalMemory()/MB_SIZE);
	}

	public long getFreeMemory(){
		//log.info("Free heap size in this runtime is: " + runtime.freeMemory()/MB_SIZE);
		return (runtime.freeMemory()/MB_SIZE);
	}

	public long getUsedMemory(){
		log.info("Used heap size in this runtime is: " + (runtime.totalMemory() - runtime.freeMemory())/MB_SIZE);
		return ((runtime.totalMemory() - runtime.freeMemory())/MB_SIZE);
	}

	public long getMaxMemory(){
		log.info("Maximum heap size in this runtime is: " + runtime.maxMemory()/MB_SIZE);
		return (runtime.maxMemory()/MB_SIZE);
	}
}
