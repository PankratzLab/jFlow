package org.genvisis.one.ben.fcs;

import java.io.File;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Arrays;

import org.flowcyt.cfcs.CFCSDataSet;
import org.flowcyt.cfcs.CFCSKeyword;
import org.flowcyt.cfcs.CFCSKeywords;
import org.flowcyt.cfcs.CFCSListModeData;
import org.flowcyt.cfcs.CFCSParameters;
import org.flowcyt.cfcs.CFCSSystem;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Logger;

public class FCSFileDuplicator {

  private FCSFileDuplicator(String destFile) throws MalformedURLException {
    syst = new CFCSSystem();
    syst.create(new File(destFile).toURI().toURL().toString());
  }

  CFCSSystem syst;
  CFCSDataSet dataSet;

  public static boolean createFrom(String srcFile, String destFile, Logger log, String[] gatings) {
    FCSFileDuplicator writer;
    try {
      writer = new FCSFileDuplicator(destFile);
    } catch (MalformedURLException e1) {
      log.reportException(e1);
      return false;
    }
    try {
      writer.copyKeywordsAndData(srcFile, gatings);
    } catch (MalformedURLException e) {
      log.reportException(e);
      return false;
    }

    return true;
  }

  private void copyKeywordsAndData(String srcFile, String[] gatings) throws MalformedURLException {
    CFCSSystem srcSyst = new CFCSSystem();
    File sysFile = new File(srcFile);
    URL fileURL = (sysFile).toURI().toURL();
    srcSyst.open(fileURL);

    CFCSDataSet dset = srcSyst.getDataSet(0);

    dataSet = this.syst.createDataSet(dset);

    CFCSKeywords srcKeys = dset.getKeywords();

    CFCSKeywords myKeys = this.dataSet.getKeywords();
    for (int i = 0, count = srcKeys.getCount(); i < count; i++) {
      CFCSKeyword k = srcKeys.getKeyword(i);
      if (CFCSKeywords.canAdd(k)) {
        myKeys.addKeyword(k);
      }
    }

    CFCSParameters srcParams = dset.getParameters();
    int paramsCount = srcParams.getCount();
    CFCSListModeData srcData = (CFCSListModeData) dset.getData();
    while (!srcData.isLoaded()) {
      Thread.yield();
    }

    CFCSKeyword gatingInfo = new CFCSKeyword();
    gatingInfo.setKeywordName(FCSDataLoader.GATING_KEY);
    gatingInfo.setKeywordSource(1);

    String v = ArrayUtils.toStr(gatings, ",");

    gatingInfo.setKeywordValue(v);
    myKeys.addKeyword(gatingInfo);

    CFCSListModeData myData = (CFCSListModeData) dataSet.getData();

    double[] data = ArrayUtils.doubleArray(paramsCount, Double.NaN);
    for (int i = 0, count = srcData.getCount(); i < count; i++) {
      srcData.getEventAsInTheFile(i, data);
      myData.addEvent(Arrays.copyOf(data, data.length));
    }

    this.syst.close();
    this.syst = null;
  }

}
