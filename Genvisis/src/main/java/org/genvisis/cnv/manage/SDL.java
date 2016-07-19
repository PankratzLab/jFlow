package org.genvisis.cnv.manage;

import java.util.Iterator;
import java.util.concurrent.Callable;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.Producer;

/**
 * similar to {@link MDL} but for sample data
 *
 */
public class SDL implements Iterator<Sample> {
	private Project proj;
	private String[] samples;
	private LOAD_TYPE lType;

	private SDLProducer producer;
	private WorkerTrain<Sample> train;
	private int numThreads;

	public SDL(Project proj, String[] samples, LOAD_TYPE lType, int numThreads) {
		super();
		this.proj = proj;
		this.samples = samples;
		this.lType = lType;
		this.numThreads = numThreads;
		init();
	}

	private void init() {
		this.producer = new SDLProducer(proj, samples, lType);
		this.train = new WorkerTrain<Sample>(producer, numThreads, 2, proj.getLog());
	}

	public enum LOAD_TYPE {
		/**
		 * Loads everything
		 */
		FULL_SAMPLE, /**
		 * Load genos and gc only
		 */
		PARTIAL_GENO_ONLY;
	}

	@Override
	public boolean hasNext() {
		boolean hasNext = train.hasNext();
		if (!hasNext) {
			train.shutdown();
		}
		return hasNext;
	}

	@Override
	public Sample next() {
		return train.next();
	}

	private static class SDLProducer implements Producer<Sample> {
		private Project proj;
		private String[] samples;
		private int index;
		private LOAD_TYPE lType;

		public SDLProducer(Project proj, String[] samples, LOAD_TYPE lType) {
			super();
			this.proj = proj;
			this.samples = samples;
			this.lType = lType;
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < samples.length;
		}

		@Override
		public Callable<Sample> next() {
			SDLWorker worker = new SDLWorker(proj, samples[index], lType);
			index++;
			return worker;
		}

		@Override
		public void shutdown() {

		}

	}

	private static class SDLWorker implements Callable<Sample> {
		private Project proj;
		private String sample;
		private LOAD_TYPE lType;

		public SDLWorker(Project proj, String sample, LOAD_TYPE lType) {
			super();
			this.proj = proj;
			this.sample = sample;
			this.lType = lType;
		}

		@Override
		public Sample call() throws Exception {
			if (sample != null) {
				switch (lType) {
				case FULL_SAMPLE:
					return proj.getFullSampleFromRandomAccessFile(sample);
				case PARTIAL_GENO_ONLY:
					return proj.getPartialSampleFromRandomAccessFile(sample, true, false, false, false, true);
				default:
					throw new IllegalArgumentException("Invalid load type " + lType);
				}
			} else {
				proj.getLog().reportTimeWarning("null sample requested, returning null");
				return null;
			}
		}
	}
}
