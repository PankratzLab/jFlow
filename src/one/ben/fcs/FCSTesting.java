package one.ben.fcs;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.flowcyt.cfcs.CFCSData;
import org.flowcyt.cfcs.CFCSDataSet;
import org.flowcyt.cfcs.CFCSGatingParameter;
import org.flowcyt.cfcs.CFCSGatingParameters;
import org.flowcyt.cfcs.CFCSKeyword;
import org.flowcyt.cfcs.CFCSKeywords;
import org.flowcyt.cfcs.CFCSListModeData;
import org.flowcyt.cfcs.CFCSParameter;
import org.flowcyt.cfcs.CFCSParameters;
import org.flowcyt.cfcs.CFCSSpillover;
import org.flowcyt.cfcs.CFCSSystem;

import common.Array;
import common.Files;
import common.ext;
import ejml.CommonOps;
import ejml.DenseMatrix64F;

public class FCSTesting {
    
    public static void main(String[] args) throws IOException, IllegalArgumentException, IllegalAccessException, NoSuchFieldException, SecurityException {
    
        CFCSSystem syst = new CFCSSystem();
        
        
        syst.open((new File("F:\\Flow\\test\\HRS003-48hrs_HRS003-48hrs PANEL 1-T-lymphocytes.fcs")).toURI().toURL());
        System.out.println(CFCSSystem.getVersion());
        
        CFCSDataSet dset = syst.getDataSet(0);
        
        System.out.println("Datatype: " + dset.getDatatype());
        
        CFCSData data = dset.getData();
        CFCSGatingParameters gparams = dset.getGatingParameters();
        CFCSParameters params = dset.getParameters();
        CFCSKeywords keys = dset.getKeywords();
        CFCSSpillover spillover = keys.getSpillover();
        
        System.out.println("Gating: " + gparams.getCount());
        for (int i = 0; i < gparams.getCount(); i++) {
            CFCSGatingParameter gat = gparams.getParameter(i);
            System.out.println(gat.getDetectorType());
            System.out.println(gat.getEmittedPercent());
            System.out.println(gat.getFilter());
            System.out.println(gat.getFullName());
            System.out.println(gat.getLogDecades());
            System.out.println(gat.getLogDecadesAndOffset());
            System.out.println(gat.getOffset());
            System.out.println(gat.getRange());
            System.out.println(gat.getShortName());
            System.out.println(gat.getVoltage());
            System.out.println();
            System.out.println();
        }
        
        System.out.println("Params: " + params.getCount());
        for (int i = 0; i < params.getCount(); i++) {
            CFCSParameter param = params.getParameter(i);
            try { System.out.println(" 1:" + param.getCalibration()); } catch (Exception e) {}
            try { System.out.println(" 2:" + param.getCalibrationFactor()); } catch (Exception e) {}
            try { System.out.println(" 3:" + param.getCalibrationUnit()); } catch (Exception e) {}
            try { System.out.println(" 4:" + param.getDetectorType()); } catch (Exception e) {}
            try { System.out.println(" 5:" + param.getEmittedPercent()); } catch (Exception e) {}
            try { System.out.println(" 6:" + param.getExcitationWavelength()); } catch (Exception e) {}
            try { System.out.println(" 7:" + param.getFieldSize()); } catch (Exception e) {}
            try { System.out.println(" 8:" + param.getFieldSizeString()); } catch (Exception e) {}
            try { System.out.println(" 9:" + param.getFilter()); } catch (Exception e) {}
            try { System.out.println("10:" + param.getFullName()); } catch (Exception e) {}
            try { System.out.println("11:" + param.getGain()); } catch (Exception e) {}
            try { System.out.println("12:" + param.getLaserPower()); } catch (Exception e) {}
            try { System.out.println("13:" + param.getLogDecades()); } catch (Exception e) {}
            try { System.out.println("14:" + param.getLogDecadesAndOffset()); } catch (Exception e) {}
            try { System.out.println("15:" + param.getOffset()); } catch (Exception e) {}
            try { System.out.println("16:" + param.getPreferredDisplay()); } catch (Exception e) {}
            try { System.out.println("17:" + param.getRange()); } catch (Exception e) {}
            try { System.out.println("18:" + param.getShortName()); } catch (Exception e) {}
            try { System.out.println("19:" + param.getSingleExcitationWavelength()); } catch (Exception e) {}
            try { System.out.println("20:" + param.getVoltage()); } catch (Exception e) {}
            try { System.out.println("21:" + Array.toStr(param.getExcitationWavelengthArray(), ", ")); } catch (Exception e) {}
            try { System.out.println("22:" + Array.toStr(param.getPreferredDisplayArgs())); } catch (Exception e) {}
            try { System.out.println("23:" + param.getPreferredDisplayScale().toString()); } catch (Exception e) {}
            
            System.out.println();
         }
        
        System.out.println("Keys: " + keys.getCount());
        System.out.println(keys.getOriginality());
        System.out.println(keys.getPlateId());
        System.out.println(keys.getPlateName());
        System.out.println(keys.getSampleVolume());
        System.out.println(keys.getWellId());
        for (int i = 0; i < keys.getCount(); i++) {
            CFCSKeyword keyw = keys.getKeyword(i);
            System.out.println(keyw.getKeywordName() + " = " + keyw.getKeywordValue());
        }
        
        System.out.println(((CFCSListModeData)data).getCount() + " events");
        
        CFCSListModeData listData = ((CFCSListModeData)data);
        System.out.println("Data datatype: " + listData.getType());
        
        byte[] dataBytes = data.getBytes();
        System.out.println(data.getType());
        System.out.println(dataBytes.length + " data bytes");
        System.out.println(listData.getCount() + " events");
        
        System.out.println();
        System.out.println(spillover.getParameterCount() + " spillover params");
        System.out.println(Array.toStr(spillover.getParameterNames(), ", "));
        System.out.println(spillover.getSpilloverCoefficients().length + " rows of spillover");
        System.out.println(spillover.getSpilloverCoefficients()[0].length + " columns of spillover");
        
        double[][] allData = new double[listData.getCount()][];
        for (int i = 0; i < listData.getCount(); i++) {
            double[] newData = new double[params.getCount()];
            listData.getEvent(i, newData);
            allData[i] = newData;
        }
        
//        Logicle logicle = new Logicle(262144d, 1, 4.5, 2);
//        System.out.println();
//        System.out.println(params.getParameter(5).getShortName());
//        double[] fifthCol = Matrix.extractColumn(allData, 5);
//        double[] scaled = new double[fifthCol.length];
//        for (int i = 0; i < fifthCol.length; i++) {
//            scaled[i] = logicle.scale(fifthCol[i]);
//        }
        
//        double[] gain = new double[params.getCount()];
//        for (int i = 0; i < gain.length; i++) {
//            gain[i] = params.getParameter(i).getGain();
//        }
//        for (int i = 0; i < allData.length; i++) {
//            for (int j = 0; j < allData[i].length; j++) {
//                allData[i][j] *= gain[j];
//            }
//        }
        
//        double[] firstRow = allData[0];
//        double[] newRow = new double[firstRow.length];
//        for (int i = 0; i < firstRow.length; i++) {
//            newRow[i] = params.getParameter(i).scaleToChannel(firstRow[i]);
//        }
        
        
        String[] origParamNames = new String[params.getCount()];
        for (int i = 0; i < origParamNames.length; i++) {
            origParamNames[i] = params.getParameter(i).getShortName();
        }
        
        DenseMatrix64F spillMatrix = getInvertedSpilloverMatrix(spillover);
        double[][] compensated = compensate(origParamNames, allData, spillover.getParameterNames(), spillMatrix);
        
        
        String[] paramNamesInOutputOrder = {"FSC-A","SSC-A","APC-A","APC-Cy7-A","BB515-A","BUV 395-A","BUV 737-A","PE-A","PE-CF594-A","PE-Cy7-A","Time"};
        int[] paramIndices = new int[paramNamesInOutputOrder.length];
        for (int p = 0; p < paramNamesInOutputOrder.length; p++) {
            for (int i = 0, count = params.getCount(); i < count; i++) {
                if (params.getParameter(i).getShortName().equals(paramNamesInOutputOrder[p])) {
                    paramIndices[p] = i;
                    break;
                }
            }
        }
        
        PrintWriter writer = Files.getAppropriateWriter("F:/Flow/test/compExport.csv");
        for (int i = 0, count = params.getCount(); i < count; i++) {
            writer.print(params.getParameter(paramIndices[i]).getShortName());
            if (i < count - 1) writer.print(",");
        }
        writer.println();
        for (double[] dataRow : compensated) {
            for (int i = 0, count = params.getCount(); i < count; i++) {
                writer.print(ext.formDeci(dataRow[paramIndices[i]], 10));
                if (i < count - 1) writer.print(",");
            }
            writer.println();
        }
        writer.flush();
        writer.close();
        System.out.println("test");
    }
    
    private static DenseMatrix64F getInvertedSpilloverMatrix(CFCSSpillover spillover) {
        double[][] coeffs = spillover.getSpilloverCoefficients();
        DenseMatrix64F spillMatrix = new DenseMatrix64F(coeffs);
        CommonOps.invert(spillMatrix);
        return spillMatrix;
    }
    
    private static double[][] compensate(String[] dataColNames, double[][] data, String[] spillColNames, DenseMatrix64F spillMatrix) {
        
        int[] spillLookup = new int[dataColNames.length];
        for (int i = 0; i < spillLookup.length; i++) {
            spillLookup[i] = ext.indexOfStr(dataColNames[i], spillColNames);
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
