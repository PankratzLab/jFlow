package org.genvisis.one.JL;

import java.util.concurrent.Callable;

import org.genvisis.common.Logger;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.Producer;

public class Qtest {

	private static class HiProducer implements Producer<Hi> {

		private int index;

		@Override
		public boolean hasNext() {
			// TODO Auto-generated method stub
			return true;
		}

		@Override
		public Callable<Hi> next() {

			Hi hi = new Hi(index == 0,index);
			index++;
			return hi;
		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub

		}

	}

	private static class Hi implements Callable<Hi> {
		private boolean halt;
		private int index;

		public Hi(boolean halt, int index) {
			super();
			this.halt = halt;
			this.index = index;
		}
		

		private int getIndex() {
			return index;
		}


		@Override
		public Hi call() throws Exception {
			System.out.println("Hi from index " + index);//this will often be printed out of order
//
			return this;
//			if (halt) {
//				while (true) {
//					try {
//						Thread.sleep(100);
//					} catch (InterruptedException ie) {
//					}
//				}
//			} else {
//				System.out.println("Hi from index " + index);
//			}
			// TODO Auto-generated method stub
//			return null;
		}

	}

	public static void main(String[] args) {
		HiProducer producer = new HiProducer();
		WorkerTrain<Hi> train = new WorkerTrain<Hi>(producer, 4, 100, new Logger());
		while (train.hasNext()) {
			System.out.println("SDFD" + train.next().getIndex());//always in order
		}

	}
}
