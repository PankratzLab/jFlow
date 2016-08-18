package org.genvisis.seq.pathway;

import java.util.HashSet;
import java.util.Set;

import org.genvisis.common.Logger;
import org.genvisis.filesys.GeneData;
import org.genvisis.filesys.LocusSet;

public class Pathway extends LocusSet<GeneData> {
  private static final long serialVersionUID = 1L;
  private final String pathwayName;
  private final Set<String> geneNames;
  private final boolean complete;

  public Pathway(String pathwayName, GeneData[] loci, boolean sort, boolean complete, Logger log) {

    super(loci, sort, log);
    this.pathwayName = pathwayName;
    geneNames = new HashSet<String>();
    this.complete = complete;
    for (GeneData element : loci) {
      geneNames.add(element.getGeneName());
    }
  }

  public Pathway(GeneData[] loci, boolean sort, Logger log, String pathwayName,
                 Set<String> geneNames, boolean complete) {
    super(loci, sort, log);
    this.pathwayName = pathwayName;
    this.geneNames = geneNames;
    this.complete = complete;
  }



  public boolean isComplete() {
    return complete;
  }

  public String getPathwayName() {
    return pathwayName;
  }

  public boolean containsGene(String gene) {
    return geneNames.contains(gene);
  }

}
