package org.genvisis.seq.analysis;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.LocusSet.TO_STRING_TYPE;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCFOps.HEADER_COPY_TYPE;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;
import org.genvisis.seq.manage.VCOps;
import org.genvisis.seq.manage.VCOps.ALT_ALLELE_CONTEXT_TYPE;
import org.genvisis.seq.manage.VCOps.VC_SUBSET_TYPE;
import org.genvisis.seq.qc.FilterNGS;
import org.genvisis.seq.qc.FilterNGS.VARIANT_FILTER_BOOLEAN;
import org.genvisis.seq.qc.FilterNGS.VARIANT_FILTER_DOUBLE;
import org.genvisis.seq.qc.FilterNGS.VariantContextFilter;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

public class RegNovo {
  private static final String[] HEADER =
      new String[] {"Sample", "ID", "REF", "ALT", "CHR", "Start", "Stop", "PassesPop", "RefDepth",
                    "AltDepth", "GQ", "Pop_AAC", "AllParents_AAC", "AllOffspring_AAC", "P1RefDepth",
                    "P1AltDepth", "P1GQ", "P1Filter", "P1_AAC", "P2RefDepth", "P2AltDepth", "P2GQ",
                    "P2Filter", "P2_AAC", "GATK_FILTER"};
  // private static final String[] TO_REPORT = new String[] { "SNPEFF_GENE_NAME", "SNPEFF_EFFECT",
  // "SNPEFF_IMPACT", "AAChange.refGene", "SNPEFF_EXON_ID", "esp6500si_all", "g10002014oct_all",
  // "culprit", "snp138" };
  public static final String OFFSPRING = "OFFSPRING";
  public static final String REG_NOVO = "REG_NOVO";

  private final String vcf;
  private final VcfPopulation vpop;
  private final VariantContextFilter filterControl;
  private final VariantContextFilter filterOffspring;
  private final FilterNGS readDepths;
  private final String bams;
  private final String outputDir;
  private final String[][] annos;
  private final Logger log;
  private final Hashtable<String, String> excludeHash;

  public RegNovo(String vcf, String bams, VcfPopulation vpop, VariantContextFilter filterControl,
                 VariantContextFilter filterOffspring, FilterNGS readDepths, String outputDir,
                 Hashtable<String, String> excludeHash, Logger log) {
    super();
    this.vcf = vcf;
    this.vpop = vpop;
    this.filterControl = filterControl;
    this.filterOffspring = filterOffspring;
    this.outputDir = outputDir;
    this.bams = bams;
    new File(outputDir).mkdirs();
    this.readDepths = readDepths;
    annos = VCFOps.getAnnotationKeys(vcf, log);
    this.excludeHash = excludeHash;

    this.log = log;
  }

  private static boolean rare(VariantContext vc, Set<String> controls, String[] annosToGet,
                              Logger log) {
    boolean rare = true;
    if (VCOps.getAAC(VCOps.getSubset(vc, controls), null) > 0) {
      rare = false;
    }
    if (rare) {
      int espIndex = ext.indexOfStr("esp6500si_all", annosToGet);
      int g1000Index = ext.indexOfStr("g10002014oct_all", annosToGet);

      if (espIndex < 0 || g1000Index < 0) {
        log.reportTimeError("esp6500si_all not being tested");
      } else {

        String[] annos = VCOps.getAnnotationsFor(annosToGet, vc, ".");
        if (!annos[espIndex].equals(".") && Double.parseDouble(annos[espIndex]) > 0.01) {
          rare = false;
        }
        if (!annos[g1000Index].equals(".") && Double.parseDouble(annos[g1000Index]) > 0.01) {
          rare = false;
        }
      }
    }

    return rare;
  }

  public void scanForDenovo() {
    String outputVCF = outputDir + ext.rootOf(vcf) + ".denovo.vcf.gz";
    String outputVCFReg = outputDir + ext.rootOf(vcf) + ".denovo.reg.vcf.gz";

    String outputSummary = outputDir + ext.rootOf(vcf) + ".denovoSummary.txt";
    String segSummary = outputDir + "denovo_segment.txt";

    if (!vpop.getSuperPop().containsKey(OFFSPRING)
        || !vpop.getSuperPop().containsKey(VcfPopulation.CONTROL)) {
      log.reportTimeError(vpop.getFileName() + " must contain both " + OFFSPRING + " and "
                          + VcfPopulation.CONTROL);
      return;
    }
    Set<String> offspring = vpop.getSuperPop().get(OFFSPRING);
    Set<String> controls = vpop.getSuperPop().get(VcfPopulation.CONTROL);
    HashSet<String> all = new HashSet<String>();
    all.addAll(controls);
    all.addAll(offspring);

    int total = 0;
    int denovo = 0;
    int regNovoCount = 0;
    if (!Files.exists(segSummary) || !Files.exists(outputVCF)) {
      VCFFileReader reader = new VCFFileReader(new File(vcf), true);
      VariantContextWriter writer =
          VCFOps.initWriter(outputVCF, null, reader.getFileHeader().getSequenceDictionary());
      VCFOps.copyHeader(reader, writer, null, HEADER_COPY_TYPE.FULL_COPY, log);

      VariantContextWriter writerReg =
          VCFOps.initWriter(outputVCFReg, null, reader.getFileHeader().getSequenceDictionary());
      VCFOps.copyHeader(reader, writerReg, null, HEADER_COPY_TYPE.FULL_COPY, log);

      PrintWriter summaryWriter = Files.getAppropriateWriter(outputSummary);
      summaryWriter.println(Array.toStr(HEADER) + "\t" + Array.toStr(annos[0])
                            + "\tOff_Exclude\tOFF_Exclude_Notes\tP1_EXCLUDE\tP1_EXCLUDE_Notes\tP2_EXCLUDE\tP2_EXCLUDE_Notes\t"
                            + REG_NOVO);
      Hashtable<String, Integer> counts = new Hashtable<String, Integer>();
      ArrayList<Segment> segsToReview = new ArrayList<Segment>();
      for (VariantContext vc : reader) {
        total++;
        if (vc.isBiallelic()) {
          VariantContext vcAll = VCOps.getSubset(vc, all, VC_SUBSET_TYPE.SUBSET_STRICT);
          VariantContext vcControls =
              VCOps.getSubset(vcAll, controls, VC_SUBSET_TYPE.SUBSET_STRICT);

          for (String off : offspring) {
            boolean regNovo = true;
            VariantContext vcOffAlt =
                VCOps.getAltAlleleContext(VCOps.getSubset(vcAll, off, VC_SUBSET_TYPE.SUBSET_STRICT),
                                          readDepths, filterOffspring, ALT_ALLELE_CONTEXT_TYPE.ALL,
                                          log);
            if (vcOffAlt.getSampleNames().size() > 0) {
              String[] fam = vpop.getPopulationForInd(off, RETRIEVE_TYPE.SUB);
              for (String element : fam) {

                Set<String> curFam = vpop.getSubPop().get(element);
                String[] excludes = new String[6];
                // boolean hasExclude = false;

                if (excludeHash.containsKey(off)) {
                  excludes[0] = "TRUE";
                  excludes[1] = excludeHash.get(off);

                } else {
                  excludes[0] = "FALSE";
                  excludes[1] = ".";
                }

                int excludeIndex = 2;
                for (String famInd : curFam) {
                  if (!famInd.equals(off)) {
                    if (excludeHash.containsKey(famInd)) {
                      excludes[excludeIndex] = "TRUE";
                      excludeIndex++;
                      excludes[excludeIndex] = excludeHash.get(famInd);
                      excludeIndex++;
                    } else {
                      excludes[excludeIndex] = "FALSE";
                      excludeIndex++;
                      excludes[excludeIndex] = ".";
                      excludeIndex++;
                    }

                  }
                }
                if (regNovo) {
                  // regNovo = !hasExclude;
                }
                if (curFam.size() == 3) {
                  String pString = "";
                  HashSet<String> parents = new HashSet<String>();
                  for (String famInd : curFam) {
                    if (!famInd.equals(off)) {
                      HashSet<String> parent = new HashSet<String>();
                      parent.add(famInd);
                      parents.add(famInd);
                      VariantContext vcFam =
                          VCOps.getSubset(vc, parent, VC_SUBSET_TYPE.SUBSET_STRICT, false);
                      int[] alleleDepths = new int[2];
                      try {
                        alleleDepths = VCOps.getAppropriateAlleleDepths(vcFam, vcFam.getGenotype(0),
                                                                        true, log);
                      } catch (Exception e) {
                        log.reportException(e);
                      }
                      pString += "\t" + alleleDepths[0];
                      pString += "\t" + alleleDepths[1];
                      pString += "\t" + vcFam.getGenotypes().get(0).getGQ() + "";
                      pString += "\t" + filterControl.filter(vcFam).getTestPerformed() + "\t"
                                 + VCOps.getAAC(vcFam, null);
                      if (regNovo) {
                        regNovo = filterControl.filter(vcFam).passed();
                      }
                    }
                  }
                  boolean diffGeno = true;
                  GenotypesContext gcParents = VCOps.getSubset(vc, parents).getGenotypes();
                  Genotype gOff = vcOffAlt.getGenotype(0);
                  for (Genotype parentGeno : gcParents) {
                    if (gOff.sameGenotype(parentGeno)) {
                      diffGeno = false;
                    }
                  }
                  if (diffGeno && VCOps.getAAC(vcOffAlt, null) > VCOps.getAAC(vc, parents)) {
                    String passesAllControls = filterControl.filter(vcControls).getTestPerformed();
                    denovo++;
                    if (filterControl.filter(vcControls).passed()) {
                      if (regNovo) {
                        regNovo = rare(vc, controls, annos[0], log);
                        if (regNovo) {
                          regNovoCount++;
                        }
                      }
                      if (!counts.containsKey(off)) {
                        counts.put(off, 0);
                      }
                      counts.put(off, counts.get(off) + 1);
                    }
                    writer.add(vc);
                    if (regNovo) {
                    }
                    segsToReview.add(VCOps.getSegment(vc));

                    String[] anno = VCOps.getAnnotationsFor(annos[0], vc, ".");
                    int[] alleleDepths =
                        VCOps.getAppropriateAlleleDepths(vcOffAlt, vcOffAlt.getGenotype(0), true,
                                                         log);
                    String sum = "\t" + alleleDepths[0];
                    sum += "\t" + alleleDepths[1];
                    sum += "\t" + vcOffAlt.getGenotypes().get(0).getGQ() + "";
                    sum += "\t" + VCOps.getAAC(VCOps.getAltAlleleContext(vcAll, readDepths,
                                                                         filterControl,
                                                                         ALT_ALLELE_CONTEXT_TYPE.ALL,
                                                                         log),
                                               all);
                    sum += "\t" + VCOps.getAAC(VCOps.getAltAlleleContext(vcControls, readDepths,
                                                                         filterControl,
                                                                         ALT_ALLELE_CONTEXT_TYPE.ALL,
                                                                         log),
                                               null);
                    sum += "\t" + VCOps.getAAC(VCOps.getAltAlleleContext(vcOffAlt, readDepths,
                                                                         filterControl,
                                                                         ALT_ALLELE_CONTEXT_TYPE.ALL,
                                                                         log),
                                               null);
                    summaryWriter.println(off + "\t" + vcOffAlt.getID() + "\t"
                                          + vcOffAlt.getReference() + "\t"
                                          + vcOffAlt.getAlternateAlleles().toString() + "\t"
                                          + vc.getContig() + "\t" + vc.getStart() + "\t"
                                          + vc.getEnd() + "\t" + passesAllControls + sum + pString
                                          + "\t" + vc.getFilters().toString() + "\t"
                                          + Array.toStr(anno) + "\t" + Array.toStr(excludes) + "\t"
                                          + regNovo);
                    summaryWriter.flush();
                    if (regNovo) {
                      writerReg.add(vc);
                    }
                  }
                }
              }
            }
          }
          if (total % 10000 == 0) {
            log.reportTimeInfo("Scanned " + total + " variants and detected" + denovo
                               + " denovos (regNovo=" + regNovoCount + ")");
          }
        }
      }
      reader.close();
      writer.close();
      writerReg.close();
      summaryWriter.close();
      LocusSet<Segment> set =
          new LocusSet<Segment>(segsToReview.toArray(new Segment[segsToReview.size()]), false,
                                log) {
            private static final long serialVersionUID = 1L;
          };
      set.writeRegions(segSummary, TO_STRING_TYPE.REGULAR, false, log);
      String countSummary = ext.addToRoot(outputSummary, ".counts");
      PrintWriter writer2 = Files.getAppropriateWriter(countSummary);
      writer2.println("OFFSPRING\tDENOVO_EXTRA");
      for (String off : counts.keySet()) {
        writer2.println(off + "\t" + counts.get(off));
      }
      writer2.close();
    } else {
      log.reportFileExists(segSummary);
      log.reportTimeWarning("Skipping denovo detection");
    }

    if (bams != null) {
      VCFOps.extractSegments(outputVCF, segSummary, 100, bams,
                             ext.parseDirectoryOfFile(vpop.getFileName()) + "extractedDenovo/",
                             true, true, 2, log);
    }
  }

  public void scanAndReportSelection(Segment[] segs, FilterNGS filterNGS,
                                     VariantContextFilter variantContextFilter, String root) {
    String output = outputDir + root + ".extracted.txt";
    String outputVCF = outputDir + root + ".extracted.vcf";
    Set<String> offspring = vpop.getSuperPop().get(OFFSPRING);
    Set<String> controls = vpop.getSuperPop().get(VcfPopulation.CONTROL);
    Set<String> exclude = vpop.getSuperPop().get(VcfPopulation.EXCLUDE);
    Set<String> all = new HashSet<String>();
    all.addAll(offspring);
    all.addAll(controls);
    all.addAll(exclude);

    Segment[] segsSorted = Segment.sortSegments(segs);
    vpop.getLog().reportTimeInfo("Scanning " + segsSorted.length + " segments");

    String alleleSummary = ext.addToRoot(output, ".alleleSummary");
    String fullFamSummary = ext.addToRoot(output, ".fullFam");
    PrintWriter writerFullFam = Files.getAppropriateWriter(fullFamSummary);
    PrintWriter writerAllele = Files.getAppropriateWriter(alleleSummary);
    writerAllele.println("CHR\tSTART\tSTOP\tID\tChildAlleles\tP1Alleles\tP2Alleles\tREF");

    PrintWriter writer = Files.getAppropriateWriter(output);
    System.out.println(all.size());
    String header =
        "CHR\tSTART\tSTOP\tREF\tALT\tID\tNUM_OFFSPRING\tAAC_OFFSPRING\tHQ_NUMOFFSPRING\tHQ_AAC_OFFSPRING\tNUM_RENTS\tAAC_RENTS\tHQ_NUM_RENTS\tHQ_AAC_RENTS\t";
    header += "NUM_EXCLUDE\tAAC_EXCLUDE\t";
    header += "GATK_FILTER\tGENVISIS_HQ_All\tGENVISIS_HQ_All_ALT";
    header += "\t" + Array.toStr(annos[0]) + "\tINFO";

    writer.println(header);
    writerFullFam.println("OFF\tRENTS\tEXCLUDED\t" + header);

    VCFFileReader reader = new VCFFileReader(new File(vcf), true);
    int count = 0;
    VariantContextWriter writerVC =
        VCFOps.initWriter(outputVCF, null, reader.getFileHeader().getSequenceDictionary());

    VCFOps.copyHeader(reader, writerVC, null, HEADER_COPY_TYPE.FULL_COPY, log);
    ArrayList<Segment> segsToReview = new ArrayList<Segment>();
    for (Segment element : segsSorted) {
      CloseableIterator<VariantContext> iter =
          reader.query(Positions.getChromosomeUCSC(element.getChr(), true), element.getStart(),
                       element.getStop());
      while (iter.hasNext()) {
        VariantContext vc = iter.next();
        count++;
        if (count % 100000 == 0) {
          log.reportTimeInfo("count " + count);
        }
        if (VCOps.isInTheseSegments(vc, segsSorted)) {
          VariantContext vcSub = VCOps.getSubset(vc, all);

          if (!vcSub.isMonomorphicInSamples()) {
            Segment vcSeg = VCOps.getSegment(vc);
            VariantContext vcOff =
                VCOps.getSubset(vc, offspring, VC_SUBSET_TYPE.SUBSET_STRICT, false);
            VariantContext vcRents =
                VCOps.getSubset(vc, controls, VC_SUBSET_TYPE.SUBSET_STRICT, false);
            VariantContext vcExclude =
                VCOps.getSubset(vc, exclude, VC_SUBSET_TYPE.SUBSET_STRICT, false);
            VariantContext vcOffHQ_ALT =
                VCOps.getAltAlleleContext(vcOff, filterNGS, variantContextFilter,
                                          ALT_ALLELE_CONTEXT_TYPE.ALL, log);
            VariantContext vcRentsHQ_ALT =
                VCOps.getAltAlleleContext(vcRents, filterNGS, variantContextFilter,
                                          ALT_ALLELE_CONTEXT_TYPE.ALL, log);

            String out = getSummary(variantContextFilter, vc, vcSub, vcSeg, vcOff, vcRents,
                                    vcExclude, vcOffHQ_ALT, vcRentsHQ_ALT, annos[0], log);
            writerVC.add(vc);
            segsToReview.add(VCOps.getSegment(vcSub));
            writer.println(out + "\t" + vcSub.toStringWithoutGenotypes());

            for (String off : offspring) {
              String[] fam = vpop.getPopulationForInd(off, RETRIEVE_TYPE.SUB);
              for (String element2 : fam) {
                Set<String> curFam = vpop.getSubPop().get(element2);
                if (VCOps.getAAC(vc, curFam) > 0) {
                  if (curFam.size() == 3) {
                    String[] parents = new String[curFam.size() - 1];
                    HashSet<String> parentsHashSetCur = new HashSet<String>();
                    HashSet<String> hashSetExcludeCur = new HashSet<String>();

                    HashSet<String> offHashSetCur = new HashSet<String>();
                    int index = 0;
                    for (String famMember : curFam) {
                      if (!off.equals(famMember)) {
                        parents[index] = famMember;
                        parentsHashSetCur.add(famMember);

                        index++;
                      } else {
                        offHashSetCur.add(famMember);
                      }
                      if (exclude.contains(famMember)) {
                        hashSetExcludeCur.add(famMember);
                      }
                    }
                    VariantContext vcSubCur =
                        VCOps.getSubset(vc, curFam, VC_SUBSET_TYPE.SUBSET_STRICT, false);
                    VariantContext vcOffCur =
                        VCOps.getSubset(vc, offHashSetCur, VC_SUBSET_TYPE.SUBSET_STRICT, false);
                    VariantContext vcRentsCur =
                        VCOps.getSubset(vc, parentsHashSetCur, VC_SUBSET_TYPE.SUBSET_STRICT, false);
                    VariantContext vcExcludeCur =
                        VCOps.getSubset(vc, hashSetExcludeCur, VC_SUBSET_TYPE.SUBSET_STRICT, false);
                    VariantContext vcOffHQ_ALTCur =
                        VCOps.getAltAlleleContext(vcOffCur, filterNGS, variantContextFilter,
                                                  ALT_ALLELE_CONTEXT_TYPE.ALL, log);
                    VariantContext vcRentsHQ_ALTCur =
                        VCOps.getAltAlleleContext(vcRentsCur, filterNGS, variantContextFilter,
                                                  ALT_ALLELE_CONTEXT_TYPE.ALL, log);
                    // private static String getSummary(VariantContextFilter variantContextFilter,
                    // VariantContext vc, VariantContext vcSub, Segment vcSeg, VariantContext vcOff,
                    // VariantContext vcRents, VariantContext vcExclude, VariantContext vcOffHQ_ALT,
                    // VariantContext vcRentsHQ_ALT, Logger log) {

                    String curOut =
                        getSummary(variantContextFilter, vc, vcSubCur, vcSeg, vcOffCur, vcRentsCur,
                                   vcExcludeCur, vcOffHQ_ALTCur, vcRentsHQ_ALTCur, annos[0], log);
                    writerFullFam.println(offHashSetCur.toString() + "\t"
                                          + parentsHashSetCur.toString() + "\t"
                                          + hashSetExcludeCur.toString() + "\t" + curOut);

                    VCOps.Transmission transmission =
                        new VCOps.Transmission(off, parents[0], parents[1]);
                    transmission.parseAlleles(vc);
                    String alleleSummaryText = "";
                    alleleSummaryText += vcSeg.getChr();
                    alleleSummaryText += "\t" + vcSeg.getStart();
                    alleleSummaryText += "\t" + vcSeg.getStop();
                    alleleSummaryText += "\t" + vc.getID();
                    alleleSummaryText += "\t" + transmission.getSummary();
                    writerAllele.println(alleleSummaryText);
                  }
                }
              }
            }
          }
        }
      }
    }
    reader.close();
    writer.close();
    writerVC.close();
    String segFile = ext.rootOf(outputVCF, false) + ".segs";
    LocusSet<Segment> set =
        new LocusSet<Segment>(segsToReview.toArray(new Segment[segsToReview.size()]), false, log) {
          private static final long serialVersionUID = 1L;
        };
    set.writeRegions(segFile, TO_STRING_TYPE.REGULAR, false, log);
    VCFOps.extractSegments(outputVCF, segFile, 100, bams,
                           ext.parseDirectoryOfFile(vpop.getFileName()) + root + "_extracted/",
                           true, true, 2, log);
  }

  private static String getSummary(VariantContextFilter variantContextFilter, VariantContext vc,
                                   VariantContext vcSub, Segment vcSeg, VariantContext vcOff,
                                   VariantContext vcRents, VariantContext vcExclude,
                                   VariantContext vcOffHQ_ALT, VariantContext vcRentsHQ_ALT,
                                   String[] annos, Logger log) {
    String out = "";
    out += vcSeg.getChr();
    out += "\t" + vcSeg.getStart();
    out += "\t" + vcSeg.getStop();
    out += "\t" + vc.getReference().getDisplayString();
    out += "\t" + vc.getAlternateAlleles().toString();
    out += "\t" + vc.getID();
    out += "\t" + VCOps.getAltAlleleContext(vcOff, null, log).getSampleNames().size();
    out += "\t" + VCOps.getAAC(vcOff, null);
    out += "\t" + vcOffHQ_ALT.getSampleNames().size();
    out += "\t" + VCOps.getAAC(vcOffHQ_ALT, null);

    out += "\t" + VCOps.getAltAlleleContext(vcRents, new FilterNGS(), log).getSampleNames().size();
    out += "\t" + VCOps.getAAC(vcRents, null);
    out += "\t" + vcRentsHQ_ALT.getSampleNames().size();
    out += "\t" + VCOps.getAAC(vcRentsHQ_ALT, null);

    out +=
        "\t" + VCOps.getAltAlleleContext(vcExclude, new FilterNGS(), log).getSampleNames().size();
    out += "\t" + VCOps.getAAC(vcExclude, null);

    out += "\t" + vc.getFilters().toString();
    out += "\t" + variantContextFilter.filter(vcSub).getTestPerformed();
    out +=
        "\t" + variantContextFilter.filter(VCOps.getAltAlleleContext(vcSub, new FilterNGS(), log))
                                   .getTestPerformed();
    out += "\t" + Array.toStr(VCOps.getAnnotationsFor(annos, vcSub, "."));
    return out;
  }

  private static VariantContextFilter getQualityFilter(double gqVal, Logger log) {
    VARIANT_FILTER_DOUBLE callRate = VARIANT_FILTER_DOUBLE.CALL_RATE;
    VARIANT_FILTER_DOUBLE gq = VARIANT_FILTER_DOUBLE.GQ;
    gq.setDFilter(gqVal);
    VARIANT_FILTER_DOUBLE dp = VARIANT_FILTER_DOUBLE.DP;
    dp.setDFilter(8);
    // VARIANT_FILTER_DOUBLE vqslod = VARIANT_FILTER_DOUBLE.VQSLOD_LOOSE;
    // vqslod
    VARIANT_FILTER_BOOLEAN fail = VARIANT_FILTER_BOOLEAN.FAILURE_FILTER;
    VARIANT_FILTER_DOUBLE[] qualFilts = new VARIANT_FILTER_DOUBLE[] {callRate, gq};
    VariantContextFilter vContextFilter =
        new VariantContextFilter(qualFilts, new VARIANT_FILTER_BOOLEAN[] {fail}, null, null, log);
    return vContextFilter;
  }

  public static void detectDenovo(String vcf, String vpopFile, String bams, double gqVal,
                                  String excludeFile, Logger log) {
    VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.ANY, log);
    vpop.report();
    Hashtable<String, String> excludeHash = new Hashtable<String, String>();
    if (excludeFile != null) {
      excludeHash =
          HashVec.loadFileToHashString(excludeFile, "Sample", new String[] {"Notes"}, "\t");
      String[] samps = VCFOps.getSamplesInFile(vcf);
      for (String samp : excludeHash.keySet()) {
        if (ext.indexOfStr(samp, samps) < 0) {
          throw new IllegalArgumentException("Error - could not find sample " + samp + " in "
                                             + vcf);
        }
      }
    }
    FilterNGS readDepths = new FilterNGS(0, 0, null);
    readDepths.setAltAlleleDepthFilter(new int[] {5});
    RegNovo regNovo =
        new RegNovo(vcf, bams, vpop, getQualityFilter(gqVal, log), getQualityFilter(gqVal, log),
                    readDepths, ext.parseDirectoryOfFile(vpopFile), excludeHash, log);
    // Segment[] segs = new Segment[] { new Segment("chr9:14079842-14400982"), new
    // Segment("chr17:7571520-7590968"), new Segment("chr8:119933796-119966383") };
    Segment[] segs = new Segment[] {new Segment("chr17:7571520-7590968")};
    regNovo.scanForDenovo();
    regNovo.scanAndReportSelection(segs, readDepths, getQualityFilter(gqVal, log),
                                   "TP53_NFIB_TNFRSF11B");

    //

  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String vcf =
        "D:/data/Project_Tsai_21_25_26_28_spector/joint_genotypes_tsai_21_25_26_28_spector.AgilentCaptureRegions.SNP.recal.INDEL.recal.merge_ARIC.hg19_multianno.eff.gatk.anno_charge.sed1000g.vcf.gz";
    String vpopFile =
        "D:/data/logan/OSv2_seq/RegNovo/noExcludeDenovoV2/OsSamps.vcfPop_noExclude2.txt";
    String excludeFile = "D:/data/logan/OSv2_seq/RegNovo/noExcludeDenovoV2/excludeSamps.txt";
    String bams = null;
    String logfile = null;
    Logger log;
    double gqVal = 40;

    String usage = "\n" + "seq.analysis.RegNovo requires 0-1 arguments\n";
    usage += "   (1) full path to a vcf file (i.e. vcf=" + vcf + " (default))\n" + "";
    usage += "   (2) full path to a vpop file (i.e. vpop=" + vpopFile + " (default))\n" + "";
    usage += "   (3) full path to a directory or file of bams (i.e. bams=" + bams
             + " (no default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("vcf=")) {
        vcf = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("vpop=")) {
        vpopFile = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("bams=")) {
        bams = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      log = new Logger(logfile);
      detectDenovo(vcf, vpopFile, bams, gqVal, excludeFile, log);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
