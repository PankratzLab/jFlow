package one.ben.fcs;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.LinkedHashSet;

import one.ben.fcs.AbstractPanel2.AXIS_SCALE;

import org.flowcyt.cfcs.CFCSData;
import org.flowcyt.cfcs.CFCSDataSet;
import org.flowcyt.cfcs.CFCSKeywords;
import org.flowcyt.cfcs.CFCSListModeData;
import org.flowcyt.cfcs.CFCSParameter;
import org.flowcyt.cfcs.CFCSParameters;
import org.flowcyt.cfcs.CFCSSpillover;
import org.flowcyt.cfcs.CFCSSystem;

import java.util.HashMap;

import common.Files;
import common.Matrix;
import common.ext;
import ejml.CommonOps;
import ejml.DenseMatrix64F;

public class FCSDataLoader {
    
    private static final String COMPENSATED_PREPEND = "Comp-";

    public static void main(String[] args) {
        String fcsFilename = "F:\\Flow\\P1-B&C-CD3-APC-Cy7 or CD4-APC-Cy7_ULTRA BRIGHT RAINBOW BEADS_URB_001.fcs";
        try {
            (new FCSDataLoader()).loadData(fcsFilename);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    volatile boolean loaded = false;
    ArrayList<String> paramNamesInOrder;
    ArrayList<AXIS_SCALE> scales;
    LinkedHashSet<String> compensatedNames;
    HashMap<String, Integer> compensatedIndices;
    int eventCount = -1;
    double[][] allData;
    double[][] compensatedData;
    String loadedFile = null;
    boolean eventOrder = true;
    
    enum DATA_SET {
        ALL,
        COMPENSATED,
        UNCOMPENSATED;
    }
    
    public FCSDataLoader() {
        paramNamesInOrder = new ArrayList<String>(); 
        scales = new ArrayList<AXIS_SCALE>();
        compensatedNames = new LinkedHashSet<String>();
        compensatedIndices = new HashMap<String, Integer>();
    }
    
    public void loadData(String fcsFilename) throws IOException {
        if (loaded) {
            throw new RuntimeException("Error - data already loaded!");
        }
        loadedFile = fcsFilename;
//        double[][] spilloverMatrix;
        
        CFCSSystem syst = new CFCSSystem();
        syst.open((new File(fcsFilename)).toURI().toURL());
        System.out.println("Data Sets: " + syst.getCount());
        
        CFCSDataSet dset = syst.getDataSet(0);
        CFCSData data = dset.getData();
        CFCSParameters params = dset.getParameters();
        CFCSKeywords keys = dset.getKeywords();
        CFCSSpillover spillover = keys.getSpillover();
        
        String[] arr = spillover.getParameterNames();
        for (int i = 0, count = arr.length; i < count; i++) {
            compensatedNames.add(arr[i]);
            compensatedIndices.put(arr[i], i);
        }
//        spilloverMatrix = spillover.getSpilloverCoefficients();
        
        for (int i = 0; i < params.getCount(); i++) {
            CFCSParameter param = params.getParameter(i);
            String name = null;
            try { name = param.getShortName(); } catch (Exception e) {}
            if (name == null) {
                try { name = param.getFullName(); } catch (Exception e) {}
            }
            if (name == null) {
                // TODO Error, no name set (included in spillover matrix?) -- will this scenario ever happen?
                name = "P" + i;
            }
            paramNamesInOrder.add(name);
            
            String axisKeywork = "P" + (i + 1) + "DISPLAY";
            AXIS_SCALE scale = AXIS_SCALE.LIN;
            try { scale = AXIS_SCALE.valueOf(keys.getKeyword(axisKeywork).getKeywordValue()); } catch (Exception e) {
                System.err.println("Warning - no axis scale set for parameter " + paramNamesInOrder.get(i) + "; assuming a linear scale.");
            };
            scales.add(scale);
        }
        
        if (data.getType() == CFCSData.LISTMODE) {
            CFCSListModeData listData = ((CFCSListModeData)data); 
            allData = new double[listData.getCount()][];
            for (int i = 0; i < listData.getCount(); i++) {
                double[] newData = new double[params.getCount()];
                listData.getEvent/*AsInTheFile*/(i, newData); // should be getEventAsInTheFile???
                allData[i] = newData;
            }
            eventCount = allData.length;
            compensatedData = compensateSmall(paramNamesInOrder, allData, spillover.getParameterNames(), getInvertedSpilloverMatrix(spillover));
            allData = Matrix.transpose(allData);
            compensatedData = Matrix.transpose(compensatedData);
            System.gc();
            eventOrder = false;
            // now in param-major-order instead of event-major-order
        } else {
            System.err.println("Error - UNSUPPORTED DATA TYPE.");
        }
        loaded = true;
    }
    
    public double[] getData(String columnName) {
        if (columnName.startsWith(COMPENSATED_PREPEND)) {
            return compensatedData[compensatedIndices.get(columnName.substring(COMPENSATED_PREPEND.length()))];
        } else {
            return allData[paramNamesInOrder.indexOf(columnName)];
        }
    }
    
    public ArrayList<String> getAllDisplayableNames(DATA_SET set) {
        ArrayList<String> names = new ArrayList<String>();
        switch (set) {
            case ALL:
                names = new ArrayList<String>(paramNamesInOrder);
                for (int i = 0; i < names.size(); i++) {
                    if (compensatedNames.contains(names.get(i))) {
                        names.add(COMPENSATED_PREPEND + names.get(i));
                    }
                }
                break;
            case COMPENSATED:
                names = new ArrayList<String>(paramNamesInOrder);
                for (int i = 0; i < names.size(); i++) {
                    if (compensatedNames.contains(names.get(i))) {
                        names.set(i, COMPENSATED_PREPEND + names.get(i));
                    }
                }
                break;
            case UNCOMPENSATED:
                names = new ArrayList<String>(paramNamesInOrder);
                break;
            default:
                break;
        }
        return names;
    }

    public void exportData(String outputFilename, DATA_SET set) {
        PrintWriter writer = Files.getAppropriateWriter(outputFilename);
        ArrayList<String> colNames = getAllDisplayableNames(set);
        for (int i = 0, count = colNames.size(); i < count; i++) {
            writer.print(colNames.get(i));
            if (i < count - 1)
                writer.print(",");
        }
        writer.println();
        
        switch (set) {
            case ALL:
                for (int i = 0; i < eventCount; i++) {
                    for (int p = 0; p < allData.length; p++) {
                        writer.print(ext.formDeci(allData[p][i], 10));
                        writer.print(",");
                    }
                    for (int p = allData.length, count = colNames.size(); p < count; p++) {
                        writer.print(ext.formDeci(compensatedData[p - allData.length][i], 10));
                        if (p < count - 1)
                            writer.print(",");
                    }
                    writer.println();
                }
                
                break;
            case COMPENSATED:
                for (int i = 0; i < eventCount; i++) {
                    for (int p = 0, c = 0; p < allData.length; p++) {
                        if (compensatedNames.contains(paramNamesInOrder.get(p))) {
                            writer.print(ext.formDeci(compensatedData[c++][i], 10));
                        } else {
                            writer.print(ext.formDeci(allData[p][i], 10));
                        }
                        if (p < allData.length - 1)
                            writer.print(",");
                    }
                    writer.println();
                }
                
                break;
            case UNCOMPENSATED:
                for (int i = 0; i < eventCount; i++) {
                    for (int p = 0; p < allData.length; p++) {
                        writer.print(ext.formDeci(allData[p][i], 10));
                        if (p < allData.length - 1)
                            writer.print(",");
                    }
                    writer.println();
                }
                
                break;
            default:
                break;
            
        }
        writer.flush();
        writer.close();
    }

    private static DenseMatrix64F getInvertedSpilloverMatrix(CFCSSpillover spillover) {
        double[][] coeffs = spillover.getSpilloverCoefficients();
        DenseMatrix64F spillMatrix = new DenseMatrix64F(coeffs);
        CommonOps.invert(spillMatrix);
        return spillMatrix;
    }
    
    private static double[][] compensateSmall(ArrayList<String> dataColNames, double[][] data, String[] spillColNames, DenseMatrix64F spillMatrix) {
        int[] spillLookup = new int[dataColNames.size()];
        for (int i = 0; i < spillLookup.length; i++) {
            spillLookup[i] = ext.indexOfStr(dataColNames.get(i), spillColNames);
        }
        
        double[][] compensated = new double[data.length][];
        for (int r = 0; r < data.length; r++) {
            compensated[r] = new double[spillColNames.length];
            for (int c = 0, cS = 0; c < data[r].length; c++) {
                if (spillLookup[c] == -1) {
                    continue;
                }
                double sum = 0d;
                for (int c1 = 0; c1 < data[r].length; c1++) {
                    if (spillLookup[c1] == -1) continue;
                    sum += data[r][c1] * spillMatrix.get(spillLookup[c1], spillLookup[c]);
                }
                compensated[r][cS++] = sum;
            }
        }
        
        return compensated;
    }
    
    private static double[][] compensate(ArrayList<String> dataColNames, double[][] data, String[] spillColNames, DenseMatrix64F spillMatrix) {
        int[] spillLookup = new int[dataColNames.size()];
        for (int i = 0; i < spillLookup.length; i++) {
            spillLookup[i] = ext.indexOfStr(dataColNames.get(i), spillColNames);
        }
        
        double[][] compensated = new double[data.length][];
        for (int r = 0; r < data.length; r++) {
            compensated[r] = new double[data[r].length];
            for (int c = 0; c < data[r].length; c++) {
                if (spillLookup[c] == -1) {
                    compensated[r][c] = data[r][c];
                    continue;
                }
                double sum = 0d;
                for (int c1 = 0; c1 < data[r].length; c1++) {
                    if (spillLookup[c1] == -1) continue;
                    sum += data[r][c1] * spillMatrix.get(spillLookup[c1], spillLookup[c]);
                }
                compensated[r][c] = sum;
            }
        }
        
        return compensated;
    }
    
    
    
}
