package org.genvisis.cnv.seq;

import org.genvisis.cnv.Resources;
import org.genvisis.cnv.Resources.Genome;
import org.genvisis.seq.GenomeBuild;
import org.genvisis.seq.ReferenceGenome;
import org.pankratzlab.common.Logger;

/**
 * Convenience {@link ReferenceGenome} to construct with a {@link GenomeBuild}'s
 * {@link Genome#getFASTA()} file
 */
public class FastaGenome extends ReferenceGenome {

  public FastaGenome(GenomeBuild build, Logger log) {
    super(Resources.genome(build, log).getFASTA().get(), log);
  }

}
