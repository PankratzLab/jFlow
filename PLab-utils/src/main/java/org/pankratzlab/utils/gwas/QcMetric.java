package org.pankratzlab.utils.gwas;

import org.pankratzlab.common.stats.Maths;
import com.google.common.base.Joiner;

public enum QcMetric {
  CHR("Chromosome Number"),
  MAF("Minor Allele Frequency"),
  CALLRATE("Marker Callrate"),
  HWE("p-value for Hardy-Weinberg Equilibrium"),
  MISHAP_HETERO("p-value for heterozygous flanking mishaps"),
  MISHAP_MIN("p-value for all mishaps"),
  P_MISS("p-value for case-control association with missingess"),
  P_GENDER("p-value for sex association"),
  P_GENDER_MISS("p-value for sex association with missigness");

  private final String key;
  private final String description;

  QcMetric(String description) {
    this(description, null);
  }

  QcMetric(String description, String key) {
    this.description = description;
    this.key = key != null ? key : this.name().toLowerCase();
  }

  public String getKey() {
    return key;
  }

  public String getDescription() {
    return description;
  }

  public String getCLIDescription() {
    String userDescription = getUserDescription();
    return userDescription.substring(0, 1).toLowerCase() + userDescription.substring(1);
  }

  public String getUserDescription() {
    StringBuilder descriptionBuilder = new StringBuilder();
    descriptionBuilder.append("Threshold to reject ").append(getDescription()).append(", using: ");
    descriptionBuilder.append(Joiner.on(", ").join(Maths.OPERATORS));
    descriptionBuilder.append(" (<0 to not filter)");
    return descriptionBuilder.toString();
  }

}
