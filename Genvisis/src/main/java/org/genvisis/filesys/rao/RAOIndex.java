package org.genvisis.filesys.rao;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Hashtable;

import org.genvisis.common.Logger;
import org.genvisis.common.SerializedFiles;

/**
 * Serialized index for accessing java objects
 *
 */
public class RAOIndex implements RAObject {
  /**
   * 
   */
  private static final long serialVersionUID = 1L;
  private String indexFileName;
  private Hashtable<String, ArrayList<Long>> index;
  private long maxSize;

  public String getFileName() {
    return indexFileName;
  }

  public RAOIndex(String fileName, Hashtable<String, ArrayList<Long>> index) {
    super();
    indexFileName = fileName;
    this.index = index;
    maxSize = 0;
  }

  public void setFileName(String fileName) {
    indexFileName = fileName;
  }

  public long getMaxSize() {
    return maxSize;
  }

  public void setMaxSize(long maxSize) {
    this.maxSize = maxSize;
  }

  public long[] getPostionsInOrder() {
    long[] pos = new long[index.size()];
    HashSet<Long> all = new HashSet<Long>();
    for (String key : index.keySet()) {
      all.addAll(index.get(key));
    }
    if (all.size() != pos.length) {
      throw new IllegalStateException("All objects could not be detected");
    }
    Arrays.sort(pos);
    return pos;
  }

  public Hashtable<String, ArrayList<Long>> getIndex() {
    return index;
  }

  public void setIndex(Hashtable<String, ArrayList<Long>> index) {
    this.index = index;
  }

  @Override
  public String[] getIndexKeys() {
    return new String[] {indexFileName};
  }

  public void serialize() {
    SerializedFiles.writeSerial(this, indexFileName, true);
  }

  public static RAOIndex load(String indexFileName, Logger log) {
    return (RAOIndex) SerializedFiles.readSerial(indexFileName, false, log, false, true);
  }

  public RAOIndex() {}

  // @Override
  // public void writeExternal(ObjectOutput out) throws IOException {
  // out.writeObject(indexFileName);
  // out.writeObject(index);
  // out.writeLong(maxSize);
  // // TODO Auto-generated method stub
  //
  // }
  //
  // @SuppressWarnings("unchecked")
  // @Override
  // public void readExternal(ObjectInput in) throws IOException, ClassNotFoundException {
  // this.indexFileName = (String) in.readObject();
  // this.index = (Hashtable<String, ArrayList<Long>>) in.readObject();
  // this.maxSize = in.readLong();
  // }

}
