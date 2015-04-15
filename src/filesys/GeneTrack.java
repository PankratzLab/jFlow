package filesys;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Vector;

import common.Array;
import common.Files;
import common.HashVec;
import common.Positions;
import common.Sort;

public class GeneTrack implements Serializable {
    public static final long serialVersionUID = 1L;

    private int[][]          starts;
    private GeneData[][]     genes;
    private String geneSetFilename;

    public GeneTrack(String geneSetFilename) {
        Hashtable<String, Vector<GeneData>> hash = new Hashtable<String, Vector<GeneData>>();
        Vector<GeneData> v;
        int chr;
        String[] chrs;
        GeneData[] set = GeneSet.load(geneSetFilename, false).getSet();
        int[] keys, poslar;

        for (int i = 0; i < set.length; i++) {
            chr = set[i].getChr();
            if (hash.containsKey(chr + "")) {
                v = hash.get(chr + "");
            } else {
                hash.put(chr + "", v = new Vector<GeneData>());
            }
            v.add(set[i]);
        }

        starts = new int[26][0];
        genes = new GeneData[26][0];
        chrs = HashVec.getKeys(hash);
        for (int i = 0; i < chrs.length; i++) {
            v = hash.get(chrs[i]);
            chr = Integer.parseInt(chrs[i]);
            poslar = new int[v.size()];
            for (int j = 0; j < v.size(); j++) {
                poslar[j] = v.elementAt(j).getStart();
            }
            keys = Sort.quicksort(poslar);
            starts[chr] = new int[keys.length];
            genes[chr] = new GeneData[keys.length];
            for (int j = 0; j < keys.length; j++) {
                starts[chr][j] = poslar[keys[j]];
                genes[chr][j] = v.elementAt(keys[j]);
            }
        }
    }

    public GeneData[] getBetween(String ucscLocation, int backTrack) {
        int[] position = Positions.parseUCSClocation(ucscLocation);
        return getBetween(position[0], position[1], position[2], backTrack);
    }

    public GeneData[] getBetween(int chr, int start, int stop, int backTrack) {
        Vector<GeneData> v;
        int first, last;
        Segment region;

        // System.out.println("Trying to get genes between chr"+chr+":"+start+"-"+stop);

        v = new Vector<GeneData>();
        region = new Segment((byte) chr, start, stop);
        if (chr < starts.length && starts[chr].length > 0) {
            first = Array.binarySearch(starts[chr], start, false);
            if (first == starts[chr].length) {
                return new GeneData[0];
            }
            try {
                last = Array.binarySearch(starts[chr], stop, first, starts[chr].length, false);
            } catch (Exception e) {
                System.out.println("uh oh");
                last = Array.binarySearch(starts[chr], stop, first, starts[chr].length, false);
            }
            if (last == starts[chr].length) {
                last--;
            }

            if (starts[chr][first] < start) {
                first++;
            }

            if (starts[chr][last] > stop) {
                last--;
            }

            for (int i = Math.max(0, first - backTrack); i < first; i++) {
                if (genes[chr][i].overlaps(region)) {
                    v.add(genes[chr][i]);
                }
            }
            for (int i = first; i <= last; i++) {
                v.add(genes[chr][i]);
            }
        }

        return GeneData.toGeneDataArray(v);
    }

    public int[] lookupPosition(String geneName) {
        for (int i = 1; i < 25; i++) {
            if (genes[i] != null) {
                for (int j = 0; j < genes[i].length; j++) {
                    if (genes[i][j].getGeneName().equalsIgnoreCase(geneName)) {
                        return new int[]
                            { genes[i][j].getChr(), genes[i][j].getStart(), genes[i][j].getStop() };
                    }
                }
            }
        }

        return new int[]
            { -1, -1, -1 };
    }

	/**
	 * @param geneName name of the gene
	 * @return all {@link GeneData} objects associated with this name, it appears there can be multiple
	 */
	public GeneData[] lookupAllGeneDatas(String geneName) {
		ArrayList<GeneData> geneDatas = new ArrayList<GeneData>();
		for (int i = 1; i < 25; i++) {
			if (genes[i] != null) {
				for (int j = 0; j < genes[i].length; j++) {
					if (genes[i][j].getGeneName().equalsIgnoreCase(geneName)) {
						geneDatas.add(genes[i][j]);
					}
				}
			}
		}
		return geneDatas.toArray(new GeneData[geneDatas.size()]);
	}
    
	public GeneData[][] getGenes() {
		return genes;
	}

	public String getGeneSetFilename() {
		return geneSetFilename;
	}

	public void setGeneSetFilename(String geneSetFilename) {
		this.geneSetFilename = geneSetFilename;
	}

	public void serialize(String filename) {
        Files.writeSerial(this, filename);
    }

    public static GeneTrack load(String filename, boolean jar) {
        return (GeneTrack) Files.readSerial(filename, jar, false);
    }

    // public static void main(String[] args) {
    // int numArgs = args.length;
    // String filename = GeneSet.DIRECTORY+GeneSet.REFSEQ_DB;
    //
    // String usage = "\n"+
    // "cnv.filesys.GeneTrack requires 0-1 arguments\n"+
    // "   (1) filename (i.e. file="+filename+" (default))\n"+
    // "";
    //
    // for (int i = 0; i<args.length; i++) {
    // if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
    // System.err.println(usage);
    // System.exit(1);
    // } else if (args[i].startsWith("file=")) {
    // filename = args[i].split("=")[1];
    // numArgs--;
    // }
    // }
    // if (numArgs!=0) {
    // System.err.println(usage);
    // System.exit(1);
    // }
    // try {
    // new GeneTrack(filename).serialize(GeneSet.DIRECTORY+ext.rootOf(GeneSet.REFSEQ_DB)+".gtrack");
    // } catch (Exception e) {
    // e.printStackTrace();
    // }
    // }

}
