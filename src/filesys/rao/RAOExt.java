package filesys.rao;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.zip.Deflater;
import java.util.zip.GZIPOutputStream;

public class RAOExt {

	public static ByteArrayOutputStream convertToByte(RAObject raObject) throws IOException {
		ByteArrayOutputStream bos = new ByteArrayOutputStream(0);
		ObjectOutputStream oos = new ObjectOutputStream(bos);
		
		oos.writeObject(raObject);
		oos.flush();
		oos.close();
		return bos;
	}

	public static ByteArrayOutputStream convertAndCompress(RAObject raObject) throws IOException {
		ByteArrayOutputStream bos = new ByteArrayOutputStream(0);
		GZIPOutputStream gzipOutput = new GZIPOutputStream(bos){
 
            {
                def.setLevel(Deflater.BEST_COMPRESSION);
            }
        };
 
		gzipOutput.write(convertToByte(raObject).toByteArray());
		gzipOutput.close();
		return bos;
	}
}
//
// ObjectInputStream ois = new ObjectInputStream(new GZIPInputStream(new ByteArrayInputStream(bos.toByteArray())));
// try {
// CNVariant cnv = (CNVariant)ois.readObject();
// cnv.getUCSClocation();
// System.out.println(cnv.getUCSClocation());
// } catch (ClassNotFoundException e) {
// // TODO Auto-generated catch block
// e.printStackTrace();
// }
// System.exit(1);