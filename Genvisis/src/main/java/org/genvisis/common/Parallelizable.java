package org.genvisis.common;

import java.util.Date;

public abstract class Parallelizable implements Runnable {
	
	public abstract void finalAction();

	public static void launch(Parallelizable[] threadSeeds, Logger log) {
		Thread[] threads;
		boolean wait;
		int numThreads;
		long time;
		
        time = new Date().getTime();
		
		numThreads = threadSeeds.length;
		
		threads = new Thread[numThreads];
		for (int i = 0; i<numThreads; i++) {
			threads[i] = new Thread(threadSeeds[i]);
			threads[i].start();
			try {
				Thread.sleep(100L);
			} catch (InterruptedException ex) {}
		}
		
		wait = true;
		while (wait) {
			wait = false;
			try {
				Thread.sleep(1000L);
			} catch (InterruptedException ex) {}
			for (int i = 0; i<numThreads; i++) {
				if (threads[i].isAlive()) {
					wait = true;
				}
			}
		}
		
		threadSeeds[0].finalAction();
		
		log.report(" ...finished in "+ext.getTimeElapsed(time));
	}
	
	public static String[][] splitList(String[] list, int numThreads, boolean collate) {
		String[][] threadSeeds;
		int[] counts;
		int step;
		
		step = (int)Math.ceil((double)list.length/(double)numThreads);
		threadSeeds = new String[numThreads][step];

		if (collate) {
			counts = new int[numThreads];
			for (int i = 0; i<list.length; i++) {
				threadSeeds[i % numThreads][counts[i % numThreads]] = list[i];
				counts[i % numThreads]++;
	        }
		} else {
			for (int i = 0; i<numThreads; i++) {
				for (int j = 0; j<step; j++) {
					if (i*step+j < list.length) {
						threadSeeds[i][j] = list[i*step+j];
					}
                }
            }
		}
		
		for (int i = 0; i<threadSeeds.length; i++) {
			threadSeeds[i] = Array.trimArray(threadSeeds[i]);
        }

		return threadSeeds;
	}
	
}
