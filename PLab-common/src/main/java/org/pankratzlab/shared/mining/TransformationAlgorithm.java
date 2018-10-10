package org.pankratzlab.shared.mining;

import org.pankratzlab.common.Logger;

/**
 * Marker interface for algorithms used in {@link Transformations}
 */
public interface TransformationAlgorithm {

  double[] transform(double[] array, Logger log);
}
