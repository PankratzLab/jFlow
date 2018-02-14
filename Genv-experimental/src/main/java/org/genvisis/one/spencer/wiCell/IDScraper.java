package org.genvisis.one.spencer.wiCell;

import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.stream.Collectors;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import com.google.common.collect.Lists;
import com.googlecode.charts4j.collect.Maps;

public class IDScraper {

  private final Collection<String> cellLines;

  /**
   * @param cellLines
   */
  private IDScraper(Collection<String> cellLines) {
    super();
    this.cellLines = cellLines;
  }

  Map<String, String> getCellLineToNWDMap() {
    Map<String, String> cellLineNWDMap = Maps.newHashMap();

    for (String cellLine : cellLines) {
      Document cellLinePage;
      try {
        cellLinePage = Jsoup.connect("https://www.wicell.org/home/stem-cell-lines/catalog-of-stem-cell-lines/"
                                     + cellLine + ".cmsx")
                            .get();
        String data = cellLinePage.data();
        int start = data.indexOf("\"GS") + 1;
        int end = data.indexOf("\"", start);
        String id = data.substring(start, end);
        Document ncbiPage = Jsoup.connect("https://www.ncbi.nlm.nih.gov/biosample/?term=" + id
                                          + "+GeneSTAR+(Genetic+Study+of+Atherosclerosis+Risk)")
                                 .get();
        String nwdID = null;
        for (Element nwdElement : ncbiPage.getElementsContainingOwnText("NWD")) {
          if (nwdID == null) nwdID = nwdElement.text();
          else if (!nwdID.equals(nwdElement.text())) {
            System.err.println("Bad NWD ID");
          }
        }
        if (cellLineNWDMap.putIfAbsent(cellLine, nwdID) != null) {
          System.err.println("Duplicate mapping for " + cellLine);
        }
      } catch (IOException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
        return null;
      }
      try {
        Thread.sleep(1000);
      } catch (InterruptedException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
    }
    return cellLineNWDMap;
  }

  public static void main(String[] args) {
    Map<String, String> cellLineIDMap = new IDScraper(Lists.newArrayList(HashVec.loadFileToStringArray("F:\\WiCellIDs\\CellLines.txt",
                                                                                                       false,
                                                                                                       null,
                                                                                                       false))).getCellLineToNWDMap();
    Files.writeIterable(cellLineIDMap.entrySet().stream()
                                     .map(entry -> entry.getKey() + "\t" + entry.getValue())
                                     .collect(Collectors.toList()),
                        "F:\\WiCellIDs\\idMap.txt");
  }

}
