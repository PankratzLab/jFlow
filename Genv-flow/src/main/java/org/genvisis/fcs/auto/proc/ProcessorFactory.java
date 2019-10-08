package org.genvisis.fcs.auto.proc;

public interface ProcessorFactory<T extends SampleProcessor> {

  T createProcessor(Object owner, int index);

  void cleanup(Object owner);
}
