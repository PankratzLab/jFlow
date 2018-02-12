/**
 * 
 */
package org.genvisis.seq.cnv;

import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.concurrent.Callable;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.Producer;
import org.genvisis.common.ext;

/**
 * Prepare aux files (split read and discordant) for lumpy run
 */
public class LumpyPrep implements PairedEndSVAnalysis {

  /**
   * 
   */
  private static final String FAILED = " failed";

  // NICE http://snible.org/java2/uni2java.html
  private static final String EXTRACT_SPLIT_READS_BWA_MEM = "\nimport sys\nimport getopt\nimport string\nfrom optparse import OptionParser\nimport re\n\ndef extractSplitsFromBwaMem(inFile,numSplits,includeDups,minNonOverlap):\n\u0009if inFile == \"stdin\":\n\u0009\u0009data = sys.stdin\n\u0009else:\n\u0009\u0009data = open(inFile, 'r')\n\u0009for line in data:\n\u0009\u0009split = 0\n\u0009\u0009if line[0] == '@':\n\u0009\u0009\u0009print line.strip()\n\u0009\u0009\u0009continue\n\u0009\u0009samList = line.strip().split('\\t')\n\u0009\u0009sam = SAM(samList)\n\u0009\u0009if includeDups==0 and (1024 & sam.flag)==1024:\n\u0009\u0009\u0009continue\n\u0009\u0009for el in sam.tags:\n\u0009\u0009\u0009if \"SA:\" in el:\n\u0009\u0009\u0009\u0009if(len(el.split(\";\")))<=numSplits:\n\u0009\u0009\u0009\u0009\u0009split = 1\n\u0009\u0009\u0009\u0009\u0009mate = el.split(\",\")\n\u0009\u0009\u0009\u0009\u0009mateCigar = mate[3]\n\u0009\u0009\u0009\u0009\u0009mateFlag = int(0)\n\u0009\u0009\u0009\u0009\u0009if mate[2]==\"-\": mateFlag = int(16)\n\u0009\u0009if split:\n\u0009\u0009\u0009read1 = sam.flag & 64 \n\u0009\u0009\u0009if read1 == 64: tag = \"_1\"\n\u0009\u0009\u0009else: tag=\"_2\"\n\u0009\u0009\u0009samList[0] = sam.query + tag \n\u0009\u0009\u0009readCigar = sam.cigar\n\u0009\u0009\u0009readCigarOps = extractCigarOps(readCigar,sam.flag)\n\u0009\u0009\u0009readQueryPos = calcQueryPosFromCigar(readCigarOps)\n\u0009\u0009\u0009mateCigarOps = extractCigarOps(mateCigar,mateFlag)\n\u0009\u0009\u0009mateQueryPos = calcQueryPosFromCigar(mateCigarOps)\n\u0009\u0009\u0009overlap = calcQueryOverlap(readQueryPos.qsPos,readQueryPos.qePos,mateQueryPos.qsPos,mateQueryPos.qePos)\n\u0009\u0009\u0009nonOverlap1 = 1 + readQueryPos.qePos - readQueryPos.qsPos - overlap\n\u0009\u0009\u0009nonOverlap2 = 1 + mateQueryPos.qePos - mateQueryPos.qsPos - overlap\n\u0009\u0009\u0009mno = min(nonOverlap1, nonOverlap2)\n\u0009\u0009\u0009if mno >= minNonOverlap:\n\u0009\u0009\u0009\u0009print \"\\t\".join(samList)\n\n#--------------------------------------------------------------------------------------------------\n# functions\n#--------------------------------------------------------------------------------------------------\n\nclass SAM (object):\n\u0009\"\"\"\n\u0009__very__ basic class for SAM input.\n\u0009\"\"\"\n\u0009def __init__(self, samList = []):\n\u0009\u0009if len(samList) > 0:\n\u0009\u0009\u0009self.query    = samList[0]\n\u0009\u0009\u0009self.flag     = int(samList[1])\n\u0009\u0009\u0009self.ref      = samList[2]\n\u0009\u0009\u0009self.pos      = int(samList[3])\n\u0009\u0009\u0009self.mapq     = int(samList[4])\n\u0009\u0009\u0009self.cigar    = samList[5]\n\u0009\u0009\u0009self.matRef   = samList[6]\n\u0009\u0009\u0009self.matePos  = int(samList[7])\n\u0009\u0009\u0009self.iSize    = int(samList[8])\n\u0009\u0009\u0009self.seq      = samList[9]\n\u0009\u0009\u0009self.qual     = samList[10]\n\u0009\u0009\u0009self.tags     = samList[11:]#tags is a list of each tag:vtype:value sets\n\u0009\u0009\u0009self.valid    = 1\n\u0009\u0009else:\n\u0009\u0009\u0009self.valid = 0\n\u0009\u0009\u0009self.query = 'null'\n\n\u0009def extractTagValue (self, tagID):\n\u0009\u0009for tag in self.tags:\n\u0009\u0009\u0009tagParts = tag.split(':', 2);\n\u0009\u0009\u0009if (tagParts[0] == tagID):\n\u0009\u0009\u0009\u0009if (tagParts[1] == 'i'):\n\u0009\u0009\u0009\u0009\u0009return int(tagParts[2]);\n\u0009\u0009\u0009\u0009elif (tagParts[1] == 'H'):\n\u0009\u0009\u0009\u0009\u0009return int(tagParts[2],16);\n\u0009\u0009\u0009\u0009return tagParts[2];\n\u0009\u0009return None;\n\u0009\n#-----------------------------------------------\ncigarPattern = '([0-9]+[MIDNSHP])'\ncigarSearch = re.compile(cigarPattern)\natomicCigarPattern = '([0-9]+)([MIDNSHP])'\natomicCigarSearch = re.compile(atomicCigarPattern)\n\ndef extractCigarOps(cigar,flag):\n\u0009if (cigar == \"*\"):\n\u0009\u0009cigarOps = []\n\u0009elif (flag & 0x0010):\n\u0009\u0009cigarOpStrings = cigarSearch.findall(cigar)\n\u0009\u0009cigarOps = []\n\u0009\u0009for opString in cigarOpStrings:\n\u0009\u0009\u0009cigarOpList = atomicCigarSearch.findall(opString)\n#\u0009\u0009\u0009print cigarOpList\n\u0009\u0009\u0009# \"struct\" for the op and it's length\n\u0009\u0009\u0009cigar = cigarOp(cigarOpList[0][0], cigarOpList[0][1])\n\u0009\u0009\u0009# add to the list of cigarOps\n\u0009\u0009\u0009cigarOps.append(cigar)\n\u0009\u0009\u0009cigarOps = cigarOps\n\u0009\u0009cigarOps.reverse()\n\u0009\u0009##do in reverse order because negative strand##\n\u0009else:\n\u0009\u0009cigarOpStrings = cigarSearch.findall(cigar)\n\u0009\u0009cigarOps = []\n\u0009\u0009for opString in cigarOpStrings:\n\u0009\u0009\u0009cigarOpList = atomicCigarSearch.findall(opString)\n\u0009\u0009\u0009# \"struct\" for the op and it's length\n\u0009\u0009\u0009cigar = cigarOp(cigarOpList[0][0], cigarOpList[0][1])\n\u0009\u0009\u0009# add to the list of cigarOps\n\u0009\u0009\u0009cigarOps.append(cigar)\n#\u0009\u0009\u0009cigarOps = cigarOps\n\u0009return(cigarOps)\n\ndef calcQueryPosFromCigar(cigarOps):\n\u0009qsPos = 0\n\u0009qePos = 0\n\u0009qLen  = 0\n\u0009# if first op is a H, need to shift start position \n\u0009# the opPosition counter sees if the for loop is looking at the first index of the cigar object    \n\u0009opPosition = 0  \n\u0009for cigar in cigarOps:\n\u0009\u0009if opPosition == 0 and (cigar.op == 'H' or cigar.op == 'S'):\n\u0009\u0009\u0009qsPos += cigar.length\n\u0009\u0009\u0009qePos += cigar.length\n\u0009\u0009\u0009qLen  += cigar.length\n\u0009\u0009elif opPosition > 0 and (cigar.op == 'H' or cigar.op == 'S'):\n\u0009\u0009\u0009qLen  += cigar.length\n\u0009\u0009elif cigar.op == 'M' or cigar.op == 'I':\n\u0009\u0009\u0009qePos += cigar.length\n\u0009\u0009\u0009qLen  += cigar.length\n\u0009\u0009\u0009opPosition += 1\n\u0009d = queryPos(qsPos, qePos, qLen);\n\u0009return d\n\nclass cigarOp (object):\n    \"\"\"\n    sturct to store a discrete CIGAR operations\n    \"\"\"\n    def __init__(self, opLength, op):\n        self.length = int(opLength)\n        self.op     = op\n\nclass queryPos (object):\n    \"\"\"\n    struct to store the start and end positions of query CIGAR operations\n    \"\"\"\n    def __init__(self, qsPos, qePos, qLen):\n        self.qsPos = int(qsPos)\n        self.qePos = int(qePos)\n        self.qLen  = int(qLen)\n\n\ndef calcQueryOverlap(s1,e1,s2,e2):\n\u0009o = 1 + min(e1, e2) - max(s1, s2)\n\u0009return max(0, o)\n\n###############################################\n\nclass Usage(Exception):\n\u0009def __init__(self, msg):\n\u0009\u0009self.msg = msg\u0009\u0009\n\ndef main():\n\u0009\n\u0009usage = \"\"\"%prog -i <file>\n\nextractSplitReads_BwaMem v0.1.0\nAuthor: Ira Hall\u0009\nDescription: Get split-read alignments from bwa-mem in lumpy compatible format. Ignores reads marked as duplicates.\nWorks on read or position sorted SAM input. Tested on bwa mem v0.7.5a-r405. \n\u0009\"\"\"\n\u0009parser = OptionParser(usage)\n\u0009\n\u0009parser.add_option(\"-i\", \"--inFile\", dest=\"inFile\", \n\u0009\u0009help=\"A SAM file or standard input (-i stdin).\",\n\u0009\u0009metavar=\"FILE\")\n\u0009parser.add_option(\"-n\", \"--numSplits\", dest=\"numSplits\", default=2, type = \"int\", \n\u0009\u0009help=\"The maximum number of split-read mappings to allow per read. Reads with more are excluded. Default=2\",\n\u0009\u0009metavar=\"INT\")\n\u0009parser.add_option(\"-d\", \"--includeDups\", dest=\"includeDups\", action=\"store_true\",default=0, \n\u0009\u0009help=\"Include alignments marked as duplicates. Default=False\")\n\u0009parser.add_option(\"-m\", \"--minNonOverlap\", dest=\"minNonOverlap\", default=20, type = \"int\", \n\u0009\u0009help=\"minimum non-overlap between split alignments on the query (default=20)\",\n\u0009\u0009metavar=\"INT\")\n\u0009(opts, args) = parser.parse_args()\n\u0009if opts.inFile is None:\n\u0009\u0009parser.print_help()\n\u0009\u0009print\n\u0009else:\n\u0009\u0009try:\n\u0009\u0009\u0009extractSplitsFromBwaMem(opts.inFile, opts.numSplits, opts.includeDups, opts.minNonOverlap)\n\u0009\u0009except IOError as err:\n\u0009\u0009\u0009sys.stderr.write(\"IOError \" + str(err) + \"\\n\");\n\u0009\u0009\u0009return\nif __name__ == \"__main__\":\n\u0009sys.exit(main()) \n\u0009\n";
  private static final String EXTRACT_SPLIT_READS_BWA_MEM_LICENSE = "#!/usr/bin/env python\n\u0009# The MIT License (MIT)\n\u0009# \n\u0009# Copyright (c) 2011,2012,2013,2014 Ryan M. Layer\n\u0009# Permission is hereby granted, free of charge, to any person obtaining a copy\n\u0009# of this software and associated documentation files (the \"Software\"), to deal\n\u0009# in the Software without restriction, including without limitation the rights\n\u0009# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n\u0009# copies of the Software, and to permit persons to whom the Software is\n\u0009# furnished to do so, subject to the following conditions:\n\u0009# \n\u0009# The above copyright notice and this permission notice shall be included in\n\u0009# all copies or substantial portions of the Software.\n\u0009# \n\u0009# THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n\u0009# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n\u0009# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n\u0009# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n\u0009# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n\u0009# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN\n\u0009# THE SOFTWARE.\n# Version https://github.com/arq5x/lumpy-sv/blob/649309bc079f911a40414b08b47b86fbe6775ee4/scripts/extractSplitReads_BwaMem";

  private String baseBam;
  private String discordantBam;
  private String splitterBam;
  private boolean fail;

  /**
   * @param baseBam input bam file
   * @param discordantBam discordant read bam file, generated with
   *          {@link LumpyPrep#extractDiscordantReads(String, String, Logger)}
   * @param splitterBam split read bam file, generated with
   *          {@link LumpyPrep#extractSplitReads(String, String, Logger)}
   * @param fail
   */
  public LumpyPrep(String baseBam, String discordantBam, String splitterBam, boolean fail) {
    super();
    this.baseBam = baseBam;
    this.discordantBam = discordantBam;
    this.splitterBam = splitterBam;
    this.fail = fail;
  }

  @Override
  public String getBaseBam() {
    return baseBam;
  }

  @Override
  public String getDiscordantBam() {
    return discordantBam;
  }

  @Override
  public String getSplitterBam() {
    return splitterBam;
  }

  @Override
  public boolean isFail() {
    return fail;
  }

  /**
   * Runs prep for multiple bams
   * 
   * @param bams bam files to prep
   * @param outDir where results will be placed
   * @param threads number of threads to use (can only use one per bam file)
   * @param log
   * @return
   */
  public static List<PairedEndSVAnalysis> runPrep(List<String> bams, String outDir, int threads,
                                                  Logger log) {
    List<PairedEndSVAnalysis> preps = new ArrayList<>();
    log.reportTimeWarning("Lumpy preparation assumes that the following are on your system path\nsamtools\npython (2.7)");

    WorkerTrain<LumpyPrep> train = new WorkerTrain<>(new LumpyPrepProducer(bams, outDir, log),
                                                     threads, 10, log);
    while (train.hasNext()) {
      LumpyPrep current = train.next();
      if (!current.isFail()) {
        preps.add(current);
      } else {
        log.reportTimeWarning("Could not prepare bam " + current.baseBam + " for lumpy analysis");
      }
    }
    return preps;
  }

  /**
   * For prepping multiple bams in parallel
   */
  private static class LumpyPrepProducer implements Producer<LumpyPrep> {

    private List<String> bams;
    private String outDir;
    private Logger log;
    private int index;

    /**
     * @param bams bam files to prep
     * @param outDir where results will be placed
     * @param log
     */
    private LumpyPrepProducer(List<String> bams, String outDir, Logger log) {
      super();
      this.bams = bams;
      this.outDir = outDir;
      this.log = log;
      this.index = 0;
    }

    /*
     * (non-Javadoc)
     * @see java.util.Iterator#hasNext()
     */
    @Override
    public boolean hasNext() {
      return index < bams.size();
    }

    /*
     * (non-Javadoc)
     * @see java.util.Iterator#next()
     */
    @Override
    public Callable<LumpyPrep> next() {
      if (index >= bams.size()) {
        throw new NoSuchElementException();
      }
      final String currentBam = bams.get(index);
      index++;
      return () -> run(currentBam, outDir, log);
    }

    /**
     * Runs prep for single bam
     * 
     * @param bam bam file to analyze
     * @param outDir where the results will be stored
     * @param log
     * @return
     */

    private static LumpyPrep run(String bam, String outDir, Logger log) {

      String discordantBam = extractDiscordantReads(bam, outDir, log);
      String splitterBam = extractSplitReads(bam, outDir, log);

      return new LumpyPrep(bam, discordantBam, splitterBam,
                           !Files.exists(discordantBam) || !Files.exists(splitterBam));
    }

    private static String extractSplitReads(String bamFile, String outDir, Logger log) {
      String outFile = outDir + ext.rootOf(bamFile) + ".splitters.unsorted.bam";
      String bat = outDir + ext.rootOf(bamFile) + ".splitters.unsorted.bat";
      String pyScript = outDir + ext.rootOf(bamFile) + ".splitters.unsorted.py";
      Files.write(EXTRACT_SPLIT_READS_BWA_MEM_LICENSE + "\n" + EXTRACT_SPLIT_READS_BWA_MEM,
                  pyScript);
      Files.chmod(pyScript);

      String script = "samtools view -h " + bamFile + " | ";
      script += pyScript + " -i stdin |";
      script += "samtools view -Sb - > " + outFile;

      CmdLine.prepareBatchForCommandLine(new String[] {script}, bat, true, log);
      boolean success = CmdLine.runCommandWithFileChecks(new String[] {bat}, "",
                                                         new String[] {bat, bamFile},
                                                         new String[] {outFile}, true, false, false,
                                                         log);
      if (!success) {
        log.reportError(bat + FAILED);
      }
      return sortBam(outFile, outDir, log);

    }

    private static String extractDiscordantReads(String bamFile, String outDir, Logger log) {
      String outFile = outDir + ext.rootOf(bamFile) + ".discordants.unsorted.bam";
      String bat = outDir + ext.rootOf(bamFile) + ".discordants.unsorted.bat";
      String script = "samtools view -b -F 1294 " + bamFile + " > " + outFile;

      CmdLine.prepareBatchForCommandLine(new String[] {script}, bat, true, log);
      boolean success = CmdLine.runCommandWithFileChecks(new String[] {bat}, "",
                                                         new String[] {bat, bamFile},
                                                         new String[] {outFile}, true, false, false,
                                                         log);
      if (!success) {
        log.reportError(bat + FAILED);
      }
      return sortBam(outFile, outDir, log);

    }

    private static String sortBam(String bamFile, String outDir, Logger log) {
      String outFile = outDir + ext.rootOf(ext.rootOf(bamFile)) + ".sorted.bam";
      String bat = outDir + ext.rootOf(ext.rootOf(bamFile)) + ".sorted.bat";
      String script = "samtools sort " + bamFile + " > " + outFile;

      CmdLine.prepareBatchForCommandLine(new String[] {script}, bat, true, log);
      boolean success = CmdLine.runCommandWithFileChecks(new String[] {bat}, "",
                                                         new String[] {bat, bamFile},
                                                         new String[] {outFile}, true, false, false,
                                                         log);
      if (!success) {
        log.reportError(bat + FAILED);
      }
      return outFile;
    }

    /*
     * (non-Javadoc)
     * @see org.genvisis.common.WorkerTrain.Producer#shutdown()
     */
    @Override
    public void shutdown() {
      // a
    }

  }

}
