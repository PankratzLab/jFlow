package org.genvisis.common;

import java.util.Date;

public class TestLogger {

	private TestLogger() {
		// Prevent instantiation
	}


	// Tests runtime of logging functions surrounding memory usage, not a JUnit test of outcomes
	private static void testMemoryLogRuntimes() {
		int[] useSomeMemory = new int[1000000000];
		int[][] useSomeMoreMemory = new int[2][1000000000];

		Logger log = new Logger();
		long time = new Date().getTime();;
		log.reportTimeElapsed("Load Methods in ", time);

		time = new Date().getTime();
		log.memoryFree();
		log.reportTimeElapsed("memoryFree() in ", time);

		time = new Date().getTime();
		log.memoryMax();
		log.reportTimeElapsed("memoryMax() in ", time);

		time = new Date().getTime();
		log.memoryUsed();
		log.reportTimeElapsed("memoryUsed() in ", time);

		time = new Date().getTime();
		log.memoryTotal();
		log.reportTimeElapsed("memoryTotal() in ", time);

		time = new Date().getTime();
		log.memoryPercentFree();
		log.reportTimeElapsed("memoryPercentFree() in ", time);

		time = new Date().getTime();
		log.memoryPercentTotalFree();
		log.reportTimeElapsed("memoryPercentTotalFree() in ", time);
	}

	public static void main(String... args) {
		testMemoryLogRuntimes();
	}

}
