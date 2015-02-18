package seq.qc;

import java.io.Serializable;

public class FilterNGS
  implements Serializable
{
  private static final long serialVersionUID = 1L;
  private double mappingQualityFilter;
  private double phreadScoreFilter;
  private int[] readDepthFilter;
  
  public FilterNGS(double mappingQualityFilter, double phreadScoreFilter, int[] readDepthFilter)
  {
    this.mappingQualityFilter = mappingQualityFilter;
    this.phreadScoreFilter = phreadScoreFilter;
    this.readDepthFilter = readDepthFilter;
  }
  
  public double getMappingQualityFilter()
  {
    return this.mappingQualityFilter;
  }
  
  public void setMappingQualityFilter(double mappingQualityFilter)
  {
    this.mappingQualityFilter = mappingQualityFilter;
  }
  
  public double getPhreadScoreFilter()
  {
    return this.phreadScoreFilter;
  }
  
  public void setPhreadScoreFilter(double phreadScoreFilter)
  {
    this.phreadScoreFilter = phreadScoreFilter;
  }
  
  public boolean passesPhreadScore(double phreadScore)
  {
    return phreadScore >= this.phreadScoreFilter;
  }
  
  public boolean passesMapQScore(double mapQScore)
  {
    return mapQScore >= this.mappingQualityFilter;
  }
  
  public int[] getReadDepthFilter()
  {
    return this.readDepthFilter;
  }
  
  public void setReadDepthFilter(int[] readDepthFilter)
  {
    this.readDepthFilter = readDepthFilter;
  }
}
