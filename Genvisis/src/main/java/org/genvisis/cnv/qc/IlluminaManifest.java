package org.genvisis.cnv.qc;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;

public class IlluminaManifest {

  /**
   * @param csv a .csv manifest file <br>
   *          Example Head; <br>
   *          Illumina, Inc. <br>
   *          [Heading] <br>
   *          Descriptor File Name,HumanOmni1-Quad_v1-0-Multi_H.bpm <br>
   *          Assay Format,Infinium HD Super <br>
   *          Date Manufactured,5/2/2011 <br>
   *          Loci Count ,1134514 <br>
   *          [Assay] <br>
   *          IlmnID,Name,IlmnStrand,SNP,AddressA_ID,AlleleA_ProbeSeq,AddressB_ID,AlleleB_ProbeSeq,
   *          GenomeBuild,Chr,MapInfo,Ploidy,Species,Source,SourceVersion,SourceStrand,SourceSeq,
   *          TopGenomicSeq,BeadSetID,Exp_Clusters,Intensity_Only,RefStrand <br>
   * @param array must be {@link ARRAY#ILLUMINA}
   * @param type must be {@link MarkerBlast.FILE_SEQUENCE_TYPE#MANIFEST_FILE}
   * @param output the output file, typically a projects marker positions
   * @param delimiter , typically ","
   * @param log
   * @return a String[] of marker names found in the manifest file
   */
  public static String[] extractMarkerPositionsFromCSVManifest(String csv, ARRAY array,
                                                               MarkerBlast.FILE_SEQUENCE_TYPE type,
                                                               String output, String delimiter,
                                                               Logger log) {

    if (type != MarkerBlast.FILE_SEQUENCE_TYPE.MANIFEST_FILE || array != ARRAY.ILLUMINA) {
      throw new IllegalArgumentException("This method should only be used in preparing marker positions for an "
                                         + ARRAY.ILLUMINA + " array using a "
                                         + MarkerBlast.FILE_SEQUENCE_TYPE.MANIFEST_FILE);
    }
    String[] required = new String[] {"Name", "Chr", "MapInfo"};
    ArrayList<String> markerNames = new ArrayList<>();
    try {
      BufferedReader reader = Files.getAppropriateReader(csv);
      boolean start = false;
      int[] extract = new int[required.length];
      PrintWriter writer = Files.openAppropriateWriter(output);
      while (reader.ready()) {
        String[] line = reader.readLine().trim().split(delimiter);

        if (!start
            && ArrayUtils.countIf(ext.indexFactors(required, line, true, log, false), -1) == 0) {
          start = true;
          extract = ext.indexFactors(required, line, true, log, false);
          writer.println("Name\tChr\tPosition");
        } else if (start) {
          String[] lineMP = null;
          try {
            lineMP = new String[required.length];
            lineMP[0] = line[extract[0]];
            lineMP[1] = line[extract[1]];
            lineMP[2] = line[extract[2]];
            if (lineMP[2].equals("0")) {
              lineMP[1] = "0";
            }
          } catch (ArrayIndexOutOfBoundsException aOfBoundsException) {
            log.reportTimeWarning("Skipping line " + ArrayUtils.toStr(line));
            lineMP = null;
          }
          if (lineMP != null && lineMP[0] != null) {
            markerNames.add(lineMP[0]);
            writer.println(ArrayUtils.toStr(lineMP));
          }
        }
      }

      reader.close();
      writer.close();
      if (!start) {
        throw new IllegalStateException("Could not find required header subset "
                                        + ArrayUtils.toStr(required, ",") + " in " + csv);
      }

    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + csv + "\" not found in current directory");
      return null;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + csv + "\"");
      return null;
    }
    return ArrayUtils.toStringArray(markerNames);
  }

  //--------------------------------------------------------

  /**
   * - From https://github.com/bioinformed/glu-genetics/blob/master/glu/lib/illumina.py String data
   * are encoded as a sequence of one or more length bytes followed by the specified number of data
   * bytes. The lower 7 bits of each length byte encodes the bits that comprise the length of the
   * following byte string. When the most significant bit it set, then an additional length byte
   * follows with 7 additional high bits to be added to the current length. The following string
   * lengths are accommodated by increasing sequences of length bytes: length maximum bytes length
   * ------ -------- 1 127 B 2 16 KB 3 2 MB 4 256 MB 5 32 GB
   * 
   * @param raf
   * @return
   * @throws IOException
   */
  protected static String readString(RandomAccessFile raf) throws IOException {
    int byt = raf.read();
    int len = byt;
    if (len == 0) return "";
    if ((byt & 0x80) > 0) {
      int sh = 7;
      len = byt & 0x7F;
      while ((byt & 0x80) > 0) {
        byt = raf.read();
        len += (byt & 0x7F) << sh;
        sh += 7;
      }
    }
    byte[] b = new byte[len];
    raf.read(b);
    return new String(b);
  }

  /**
   * Reads a 4-byte little-endian long
   * 
   * @param raf
   * @return
   * @throws IOException
   */
  protected static long readLong(RandomAccessFile raf) throws IOException {
    long l = raf.read();
    l += raf.read() << 8;
    l += raf.read() << 16;
    l += raf.read() << 24;
    return l;
  }

  private IlluminaManifest() {}

  public static void parseMarkerPositionsFromBPM(String file, String outFile) throws IOException {
    RandomAccessFile raf = new RandomAccessFile(file, "r");

    int b, p, m;
    b = raf.read();
    p = raf.read();
    m = raf.read();
    if (b != 'B' || p != 'P' || m != 'M') {
      try {
        raf.close();
      } catch (IOException e) {}
      throw new RuntimeException("Invalid file format signature; expected 'BPM', found '"
                                 + ((char) b) + ((char) p) + ((char) m) + "'");
    }

    int version1 = raf.read();
    long version2 = readLong(raf);

    if (version1 != 1 || version2 != 4) {
      throw new RuntimeException("Invalid file format version number; expected (1,4), found ("
                                 + version1 + "," + version2 + ")");
    }

    readString(raf); // manifest name
    readString(raf); // controls

    long count = readLong(raf);
    for (long i = 0; i < count; i++) {
      readLong(raf); // snpEntries
    }
    for (long i = 0; i < count; i++) {
      readString(raf); // names
    }
    List<Integer> normIds = new ArrayList<>();
    for (long i = 0; i < count; i++) {
      normIds.add(raf.read());
    }

    PrintWriter writer = Files.getAppropriateWriter(outFile);

    for (long i = 0; i < count; i++) {
      ManifestEntry me = ManifestEntry.parse(raf, normIds.get((int) i));
      writer.println(me.name + "\t" + (me.mapInfo.equals("0") ? 0 : me.chrom) + "\t" + me.mapInfo);
      me = null;
    }
    writer.close();

    normIds = null;

    raf.close();
  }

  /* Could be turned into an iterable */
  public static IlluminaManifest load(String file) throws IOException {
    RandomAccessFile raf = new RandomAccessFile(file, "r");
    IlluminaManifest im = new IlluminaManifest();
    int b, p, m;
    b = raf.read();
    p = raf.read();
    m = raf.read();
    if (b != 'B' || p != 'P' || m != 'M') {
      try {
        raf.close();
      } catch (IOException e) {}
      throw new RuntimeException("Invalid file format signature; expected 'BPM', found '"
                                 + ((char) b) + ((char) p) + ((char) m) + "'");
    }

    int version1 = raf.read();
    long version2 = readLong(raf);

    if (version1 != 1 || version2 != 4) {
      throw new RuntimeException("Invalid file format version number; expected (1,4), found ("
                                 + version1 + "," + version2 + ")");
    }

    im.manifestName = readString(raf);
    String c = readString(raf);
    im.controls = Arrays.asList(c.split("\n"));

    im.count = readLong(raf);
    im.snpEntries = new ArrayList<>();
    for (long i = 0; i < im.count; i++) {
      im.snpEntries.add(readLong(raf));
    }
    im.names = new ArrayList<>();
    for (long i = 0; i < im.count; i++) {
      im.names.add(readString(raf));
    }
    im.normIds = new ArrayList<>();
    for (long i = 0; i < im.count; i++) {
      im.normIds.add(raf.read());
    }

    im.entries = new ArrayList<>();
    im.nameMap = new HashMap<>();
    for (long i = 0; i < im.count; i++) {
      ManifestEntry me;
      im.entries.add(me = ManifestEntry.parse(raf, im.normIds.get((int) i)));
      im.nameMap.put(me.name, me);
    }

    raf.close();
    return im;
  }

  public void writeCSV(String file) {
    PrintWriter writer = Files.getAppropriateWriter(file);
    writer.println(ArrayUtils.toStr(HEADER, ","));
    for (ManifestEntry me : entries) {
      writer.println(me.formatCSV());
    }
    writer.close();
  }

  public static final String[] HEADER = new String[] {"IlmnID", "Name", "IlmnStrand", "SNP",
                                                      "AssayTypeID", "NormID", "AddressA_ID",
                                                      "AlleleA_ProbeSeq", "AddressB_ID",
                                                      "AlleleB_ProbeSeq", "GenomeBuild", "Chr",
                                                      "MapInfo", "Ploidy", "Species", "Source",
                                                      "SourceVersion", "SourceStrand", "SourceSeq",
                                                      "TopGenomicSeq", "CustomerStrand",
                                                      "GenomicStrand"};

  private String manifestName;
  private List<String> controls;
  private long count;
  private ArrayList<Long> snpEntries;
  private ArrayList<String> names;
  private ArrayList<Integer> normIds;
  private ArrayList<ManifestEntry> entries;
  private Map<String, ManifestEntry> nameMap;

  // TODO add accessor methods as needed 

  static class ManifestEntry {

    long versionId;
    String ilmnId;
    String name;
    long snpNum;
    String designStrand;
    String alleles;
    String chrom;
    String ploidy;
    String species;
    String mapInfo;
    String topGenomicSeq;
    String custStrand;
    long addressA, addressB;
    String alleleProbeA, alleleProbeB;
    String genomeVersion;
    String source;
    String srcVersion;
    String srcStrand;
    String srcSeq;
    int normId;
    int assayTypeId;
    String genomicStrand;
    int normId2;

    private ManifestEntry() {}

    public static ManifestEntry parse(RandomAccessFile raf, int normId) throws IOException {
      ManifestEntry me = new ManifestEntry();
      me.versionId = readLong(raf);
      if (me.versionId != 4 && me.versionId != 7 && me.versionId != 8) {
        throw new RuntimeException("Invalid BPM record version number: " + me.versionId);
      }
      me.ilmnId = readString(raf);
      me.name = readString(raf);
      raf.read();
      raf.read();
      raf.read();
      me.snpNum = readLong(raf);
      raf.read();
      // interning strings cuts memory by ~half
      me.designStrand = readString(raf).intern();
      me.alleles = readString(raf).intern();
      me.chrom = readString(raf).intern();
      me.ploidy = readString(raf).intern();
      me.species = readString(raf).intern();
      me.mapInfo = readString(raf);
      me.topGenomicSeq = readString(raf);
      me.custStrand = readString(raf).intern();
      me.addressA = readLong(raf);
      me.addressB = readLong(raf);
      me.alleleProbeA = readString(raf);
      me.alleleProbeB = readString(raf);
      me.genomeVersion = readString(raf).intern();
      me.source = readString(raf).intern();
      me.srcVersion = readString(raf).intern();
      me.srcStrand = readString(raf).intern();
      me.srcSeq = readString(raf);
      me.normId = normId;

      if (me.versionId == 7 || me.versionId == 8) {
        raf.read();
        raf.read();
        raf.read();
        me.assayTypeId = raf.read();
        readLong(raf);
        readLong(raf);
        readLong(raf);
        readLong(raf);
      } else {
        me.assayTypeId = 0;
      }

      me.genomicStrand = me.versionId == 8 ? readString(raf).intern() : "";

      me.normId2 = me.normId + 100 * me.assayTypeId + 1;

      return me;
    }

    public String formatCSV() {
      StringJoiner line = new StringJoiner(",");
      line.add(ilmnId).add(name).add(designStrand).add(alleles).add(Integer.toString(assayTypeId))
          .add(Integer.toString(normId)).add(Long.toString(addressA)).add(alleleProbeA)
          .add(Long.toString(addressB)).add(alleleProbeB).add(genomeVersion).add(chrom).add(mapInfo)
          .add(ploidy).add(species).add(source).add(srcVersion).add(srcStrand).add(srcSeq)
          .add(topGenomicSeq).add(custStrand).add(genomicStrand);
      return line.toString();
    }

  }

  public static void main(String[] args) throws IOException {
    String dir = "G:\\Manifests\\";
    String[] files = Files.list(dir, ".bpm");
    for (String f : files) {
      System.out.println("Processing " + f);
      IlluminaManifest.load(dir + f).writeCSV(dir + ext.rootOf(f, true) + ".TEST.csv");
    }
    System.out.println("Done!");
  }

  public List<ManifestEntry> getEntries() {
    return Collections.unmodifiableList(entries);
  }

}
