package org.genvisis.filesys.rao;

import java.io.ByteArrayInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.Iterator;
import java.util.zip.GZIPInputStream;

import org.genvisis.common.Logger;

public class RAOReader<T extends RAObject> implements Iterator<T> {

  private final String fullPathToFile;
  private final RAOIndex raoIndex;
  private final String indexFileName;
  private int numLoaded;
  private final FileInputStream fis;
  private long[] allPositions;
  private final Logger log;

  public RAOReader(String fullPathToFile, Logger log) throws FileNotFoundException {
    super();
    this.fullPathToFile = fullPathToFile;
    this.indexFileName = fullPathToFile + RAOWriter.RAO_INDEX_EXT;
    this.raoIndex = RAOIndex.load(indexFileName, log);
    this.log = log;
    this.fis = new FileInputStream(fullPathToFile);
    this.numLoaded = 0;

  }

  public String getFullPathToFile() {
    return fullPathToFile;
  }

  public void close() throws IOException {
    fis.close();
  }

  @Override
  public boolean hasNext() {
    if (allPositions == null) {
      this.allPositions = raoIndex.getPostionsInOrder();
      log.reportTimeInfo("Number to load :" + allPositions.length);
    }
    return numLoaded < allPositions.length;
  }

  @Override
  public T next() {
    T t = loadPosition(allPositions[numLoaded]);
    numLoaded++;
    return t;
  }

  @Override
  public void remove() {
    // TODO Auto-generated method stub

  }

  @SuppressWarnings("unchecked")
  public T loadPosition(long position) {
    // BufferedReader reader = new BufferedReader( new InputStreamReader(new ZipInputStream(new
    // FileInputStream(fullPathToFile))));
    // while(reader.ready()){
    // System.out.println(reader.read());
    // }

    try {
      System.out.println("Position " + position);
      fis.getChannel().position(position);
      byte[] data = new byte[(int) raoIndex.getMaxSize()];
      fis.read(data);
      ByteArrayInputStream bis = new ByteArrayInputStream(data);
      ObjectInputStream iis = new ObjectInputStream(new GZIPInputStream(bis));
      T t = (T) iis.readObject();
      iis.close();
      return t;
    } catch (IOException e) {
      e.printStackTrace();
    } catch (ClassNotFoundException e) {
      e.printStackTrace();
    }

    return null;
    // inputStream.

  }

  public static class RAOReaderQuery<T extends RAObject> implements Iterator<T> {

    @Override
    public boolean hasNext() {
      // TODO Auto-generated method stub
      return false;
    }

    @Override
    public T next() {
      // TODO Auto-generated method stub
      return null;
    }

    @Override
    public void remove() {
      // TODO Auto-generated method stub

    }

  }

}
