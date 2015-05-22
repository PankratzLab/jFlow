package filesys.rao;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.util.Iterator;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipInputStream;

import common.Array;
import common.Files;
import common.Logger;

public class RAOReader<T extends RAObject> {

	private String fullPathToFile;
	private RAOIndex raoIndex;
	private String indexFileName;
	private FileInputStream fis;
	private Logger log;

	public RAOReader(String fullPathToFile, String indexFileName, Logger log) {
		super();
		this.fullPathToFile = fullPathToFile;
		this.raoIndex = RAOIndex.load(indexFileName, log);
		this.indexFileName = indexFileName;
		this.log = log;
		try {
			this.fis = new FileInputStream(fullPathToFile);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	@SuppressWarnings("unchecked")
	public T loadPosition(long position) {
		// BufferedReader reader = new BufferedReader( new InputStreamReader(new ZipInputStream(new FileInputStream(fullPathToFile))));
		// while(reader.ready()){
		// System.out.println(reader.read());
		// }

		try {
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
