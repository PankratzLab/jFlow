package org.genvisis.one.JL;

import java.util.concurrent.Callable;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;

public class Qtest {

  private static class HiProducer extends AbstractProducer<Hi> {

    private int index;

    @Override
    public boolean hasNext() {
      // TODO Auto-generated method stub
      return true;
    }

    @Override
    public Callable<Hi> next() {

      Hi hi = new Hi(index == 0, index);
      index++;
      return hi;
    }
  }

  private static class Hi implements Callable<Hi> {

    private final int index;

    public Hi(boolean halt, int index) {
      super();
      this.index = index;
    }

    private int getIndex() {
      return index;
    }

    @Override
    public Hi call() throws Exception {
      System.out.println("Hi from index " + index);// this will often be printed out of order
      //
      return this;
      // if (halt) {
      // while (true) {
      // try {
      // Thread.sleep(100);
      // } catch (InterruptedException ie) {
      // }
      // }
      // } else {
      // System.out.println("Hi from index " + index);
      // }
      // TODO Auto-generated method stub
      // return null;
    }

  }

  public static void main(String[] args) {
    HiProducer producer = new HiProducer();
    try (WorkerTrain<Hi> train = new WorkerTrain<Hi>(producer, 4, 100, new Logger())) {
      while (train.hasNext()) {
        System.out.println("SDFD" + train.next().getIndex());// always in order
      }
    }

  }
}
