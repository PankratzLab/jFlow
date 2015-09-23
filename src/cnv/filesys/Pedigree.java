package cnv.filesys;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import cnv.qc.MendelErrors;
import cnv.qc.MendelErrors.MendelErrorCheck;
import cnv.var.SampleData;
import common.Logger;
import common.ext;
import filesys.FamilyStructure;



public class Pedigree extends FamilyStructure {
    
    private ArrayList<String[]> cached_poPairsIDs = null;
    private ArrayList<int[]> cached_poPairsCompleteOnly = null;
    private HashMap<String, ArrayList<String>> cached_parentToChildrenMap = null;
    private boolean cached_parentMapIsCompleteOnly = false;
    public ArrayList<String[]> cached_all_trios = null;
    public ArrayList<int[]> cached_complete_trios = null;
    public ArrayList<String[]> cached_sib_pairs = null;
    
    public static class PedigreeUtils {
        
        public static void clearCache(Pedigree ped) {
            ped.cached_poPairsIDs = null;
            ped.cached_poPairsCompleteOnly = null;
            ped.cached_parentToChildrenMap = null;
            ped.cached_parentMapIsCompleteOnly = false;
            ped.cached_all_trios = null;
            ped.cached_complete_trios = null;
            ped.cached_sib_pairs  = null;
        }
        
        public static ArrayList<String[]> loadSibs(Pedigree ped, boolean completeOnly, HashSet<String> excludedFIDIIDs, boolean cache) {
            if (ped.cached_all_trios == null) {
                loadCompleteTrios(ped, excludedFIDIIDs, true); // will also create all_trios
            }
            if (ped.cached_all_trios == null) {
                // ERROR
            }
            if (ped.cached_sib_pairs != null) {
                return ped.cached_sib_pairs;
            }
            HashMap<String, ArrayList<String>> parentToChildren = loadParentToChildrenMap(ped, completeOnly, excludedFIDIIDs, cache);
            
            // at this point, only non-excluded IDs are present in all_trios and parentToChildren
            ArrayList<String[]> sibPairs = new ArrayList<String[]>();
            for (String[] trio : ped.cached_all_trios) { 
                ArrayList<String> faChildren = parentToChildren.get(trio[1]);
                if (faChildren == null) {
                    // Error!
                } else if (faChildren.size() == 1 && !trio[0].equals(faChildren.get(0))) {
                    // Error!
                }
                ArrayList<String> moChildren = parentToChildren.get(trio[2]);
                if (moChildren == null) {
                    // Error!
                } else if (moChildren.size() == 1 && !trio[0].equals(faChildren.get(0))) {
                    // Error!
                }
                HashSet<String> unionSet = new HashSet<String>();
                unionSet.addAll(faChildren);
                unionSet.retainAll(moChildren);
                if (unionSet.size() == 0) {
                    continue; // no sibs
                } else {
                    for (String sib : unionSet) {
                        sibPairs.add(new String[]{sib, trio[0]});
                        sibPairs.add(new String[]{trio[0], sib});
                    }
                }
            }
            if (cache) {
                ped.cached_sib_pairs = sibPairs;
            }
            return sibPairs;
        }
        
        public static HashMap<String, ArrayList<String>> loadParentToChildrenMap(Pedigree ped, boolean completeOnly, HashSet<String> excludedFIDIIDs, boolean cache) {
            if (ped.cached_parentToChildrenMap != null) {
                return ped.cached_parentToChildrenMap;
            }
            HashMap<String, ArrayList<String>> parentMap = new HashMap<String, ArrayList<String>>();
            
            for (int i = 0; i < ped.getIDs().length; i++) {
                if (!FamilyStructure.MISSING_ID_STR.equals(ped.getFA(i))
                        && (excludedFIDIIDs == null || !excludedFIDIIDs.contains(ped.getFID(i) + "\t" + ped.getFA(i)))
                        && (!completeOnly || (/*faInd = */ped.getIndexOfFaInIDs(i)) >= 0)) {
                    ArrayList<String> children = parentMap.get(ped.getFID(i) + "\t" + ped.getFA(i));
                    if (children == null) {
                        children = new ArrayList<String>();
                        parentMap.put(ped.getFID(i) + "\t" + ped.getFA(i), children);
                    }
                    children.add(ped.getFID(i) + "\t" + ped.getIID(i));
                }
                if (!FamilyStructure.MISSING_ID_STR.equals(ped.getMO(i))
                        && (excludedFIDIIDs == null || !excludedFIDIIDs.contains(ped.getFID(i) + "\t" + ped.getMO(i)))
                        && (!completeOnly || (/*moInd = */ped.getIndexOfMoInIDs(i)) >= 0)) {
                    ArrayList<String> children = parentMap.get(ped.getFID(i) + "\t" + ped.getMO(i));
                    if (children == null) {
                        children = new ArrayList<String>();
                        parentMap.put(ped.getFID(i) + "\t" + ped.getMO(i), children);
                    }
                    children.add(ped.getFID(i) + "\t" + ped.getIID(i));
                }
            }
            if (cache) {
                ped.cached_parentToChildrenMap = parentMap;
                ped.cached_parentMapIsCompleteOnly = completeOnly;
            }
            return parentMap;
        }
        
        public static ArrayList<String[]> loadPOPairs(Pedigree ped, boolean completeOnly, HashSet<String> excludedFIDIIDs, boolean cache) {
            if (ped.cached_poPairsIDs != null) {
                return ped.cached_poPairsIDs;
            }
            ArrayList<String[]> pairs = new ArrayList<String[]>();
            ArrayList<int[]> completePairs = new ArrayList<int[]>();
            for (int i = 0; i < ped.getIDs().length; i++) {
                int faInd = -1;
                int moInd = -1;
                if (!FamilyStructure.MISSING_ID_STR.equals(ped.getFA(i))
                        && (excludedFIDIIDs == null || !excludedFIDIIDs.contains(ped.getFID(i) + "\t" + ped.getFA(i)))
                        && (!completeOnly || (faInd = ped.getIndexOfFaInIDs(i)) >= 0)) {
                    pairs.add(new String[]{ped.getFID(i) + "\t" + ped.getFA(i), ped.getFID(i) + "\t" + ped.getIID(i)});
                    if (completeOnly) {
                        completePairs.add(new int[]{faInd, i});
                    }
                }
                if (!FamilyStructure.MISSING_ID_STR.equals(ped.getMO(i))
                        && (excludedFIDIIDs == null || !excludedFIDIIDs.contains(ped.getFID(i) + "\t" + ped.getMO(i)))
                        && (!completeOnly || (moInd = ped.getIndexOfMoInIDs(i)) >= 0)) {
                    pairs.add(new String[]{ped.getFID(i) + "\t" + ped.getMO(i), ped.getFID(i) + "\t" + ped.getIID(i)});
                    if (completeOnly) {
                        completePairs.add(new int[]{moInd, i});
                    }
                }
            }
            if (cache) {
                ped.cached_poPairsIDs = pairs;
                if (completeOnly) {
                    ped.cached_poPairsCompleteOnly = completePairs;
                }
            }
            return pairs;
        }
        
        public static ArrayList<int[]> loadCompleteTrios(Pedigree ped, HashSet<String> excludedFIDIIDs, boolean cache) {
            if (ped.cached_complete_trios != null) {
                return ped.cached_complete_trios;
            }
            ArrayList<String[]> allTrios = new ArrayList<String[]>();
            ArrayList<int[]> trios = new ArrayList<int[]>();
            for (int i = 0; i < ped.getIDs().length; i++) {
                int faInd = -1;
                int moInd = -1;
                if (!FamilyStructure.MISSING_ID_STR.equals(ped.getFA(i)) && (excludedFIDIIDs == null || !excludedFIDIIDs.contains(ped.getFID(i) + "\t" + ped.getFA(i)))
                        && !FamilyStructure.MISSING_ID_STR.equals(ped.getMO(i)) && (excludedFIDIIDs == null || !excludedFIDIIDs.contains(ped.getFID(i) + "\t" + ped.getMO(i)))) {
                    if (cache) {
                        allTrios.add(new String[]{ped.getFID(i) + "\t" + ped.getIID(i), ped.getFID(i) + "\t" + ped.getFA(i), ped.getFID(i) + "\t" + ped.getMO(i)}); 
                    }
                    if ((faInd = ped.getIndexOfFaInIDs(i)) >= 0 && (moInd = ped.getIndexOfMoInIDs(i)) >= 0) {
                        trios.add(new int[]{i, faInd, moInd});
                    }
                }
            }
            if (cache) {
                ped.cached_complete_trios  = trios;
                ped.cached_all_trios = allTrios;
            }
            return trios;
        }
        
        public static MendelErrorCheck[] checkMendelErrors(Pedigree pedigree, MarkerData markerData, boolean[] samplesToCheck, String[] sex, ClusterFilterCollection clusterFilters, float gcThreshold) {
            if (pedigree.getProject() == null) {
                System.err.println(ext.getTime() + "]\t Error - cannot run checkMendelErrors without a Project");
                return null;
            }
            MendelErrorCheck[] mendelErrorChecks = new MendelErrorCheck[pedigree.getProject().getSamples().length];
            byte[] genotypes = markerData.getAbGenotypesAfterFilters(clusterFilters, markerData.getMarkerName(), gcThreshold);
            Logger log = pedigree.getProject().getLog();
            if (!pedigree.isProjectOrder()) {
                log.reportTimeError("Pedigree file must be in project order, internal error");
                return null;
            } else {
                for (int i = 0; i < pedigree.getIDs().length; i++) {
                    if (/*this check isn't valid*//*pedigreeEntries[i] != null && */(samplesToCheck == null || samplesToCheck[i]) && pedigree.getIDNAIndex(i) >= 0) {
                        int sampleIndex = pedigree.getIDNAIndex(i);
                        int faDNAIndex = pedigree.getFaDNAIndex(i);
                        int moDNAIndex = pedigree.getMoDNAIndex(i);
                        
                        byte faGenotype = -1;
                        if (faDNAIndex >= 0 && (samplesToCheck == null || samplesToCheck[faDNAIndex])) {
                            faGenotype = genotypes[faDNAIndex];
                        }
                        byte moGenotype = -1;
                        if (moDNAIndex >= 0 && (samplesToCheck == null || samplesToCheck[moDNAIndex])) {
                            moGenotype = genotypes[moDNAIndex];
                        }
                        int sampleSex = -1;
                        try {
                            if (sex != null) {
                                sampleSex = Integer.parseInt(sex[pedigree.getIDNAIndex(i)]);
                            }
                        } catch (NumberFormatException nfe) {
                            
                        }
                        // System.out.println(faGenotype+"\t"+moGenotype);
                        MendelErrors mendelErrors = new MendelErrors(markerData.getChr(), sampleSex, genotypes[sampleIndex], faGenotype, moGenotype);
                        mendelErrorChecks[i] = mendelErrors.checkMendelError();
                    } else {
                        mendelErrorChecks[i] = new MendelErrors(markerData.getChr(), -1, (byte) -1, (byte) -1, (byte) -1).checkMendelError();
                    }
                }
            }
            return mendelErrorChecks;
        }
        
        
    }
    
    private Project project;
    private boolean nullProject;
    private int[][] dnaIndicesInProject;
    public static final int MISSING_DNA_INDEX = -1;
    boolean projectOrder;
    
    public Pedigree(Project proj) {
        this(proj, proj.PEDIGREE_FILENAME.getValue(), true);
    }
    
    /**
     * 
     * @param proj
     * @param projectOrder (Used for ProjectUtils.checkMendelErrors()) Indicates if the order of entries in the project's pedigree file (from the PEDIGREE_FILENAME property) matches the internal order of the project samples.
     */
    public Pedigree(Project proj, boolean projectOrder) {
        this(proj, proj.PEDIGREE_FILENAME.getValue(), projectOrder);
    }
    
    /**
     * 
     * @param proj
     * @param pedigreeFile
     * @param projectOrder (Used for ProjectUtils.checkMendelErrors()) Indicates if the order of entries in the given pedigree file matches the internal order of the project samples.
     */
    public Pedigree(Project proj, String pedigreeFile, boolean projectOrder) {
        super(pedigreeFile, true, proj == null ? new Logger() : proj.getLog());
        this.project = proj;
        this.nullProject = proj == null;
        this.dnaIndicesInProject = new int[this.ids.length][];
        this.projectOrder = projectOrder;
        SampleData sampleData = nullProject ? null : proj.getSampleData(0, false);
        String[] samples = nullProject ? null : proj.getSamples();
        for (int i = 0; i < this.ids.length; i++) {
            int iDNAIndex = proj == null ? MISSING_DNA_INDEX : getSampleIndex(this.dnas[i], sampleData, samples);
            int faDNAIndex = proj == null ? MISSING_DNA_INDEX : getSampleIndex(this.ids[i][FamilyStructure.FA_INDEX], sampleData, samples);
            int moDNAIndex = proj == null ? MISSING_DNA_INDEX : getSampleIndex(this.ids[i][FamilyStructure.MO_INDEX], sampleData, samples);
            this.dnaIndicesInProject[i] = new int[]{iDNAIndex, faDNAIndex, moDNAIndex};
        }
        
    }
    
    private static int getSampleIndex(String sample, SampleData sampleData, String[] projectSamples) {
        int sampleIndex = MISSING_DNA_INDEX;
        if (!sample.equals(FamilyStructure.MISSING_ID_STR) && sampleData != null && projectSamples != null && sampleData.lookup(sample) != null) {
            sampleIndex = ext.indexOfStr(sampleData.lookup(sample)[0], projectSamples);
        }
        return sampleIndex;
    }
    
    public boolean isProjectOrder() {
        return this.projectOrder;
    }
    
    public Project getProject() {
        return this.project;
    }
    
    public int getIDNAIndex(int index) {
        return this.dnaIndicesInProject[index][0];
    }
    public int getFaDNAIndex(int index) {
        return this.dnaIndicesInProject[index][1];
    }
    public int getMoDNAIndex(int index) {
        return this.dnaIndicesInProject[index][2];
    }
    
}